# SPDX-License-Identifier: GPL-3.0-only
"""zeros25k_rerun.py — in-house recreation and re-run on the 25,000 exact zeros of ζ.   Run disowned:  see experiments/zeros25k/README.md

Stage A  RECREATE   every γ_k of the SemPriCon table re-found from scratch at 32 digits (mpmath siegelz root, secant, tol 1e-28), then
                    · sign change of Z(t) across γ_k ± 1e-9 (a simple real root, not a touch)       :CALCULATED
                    · strictly increasing, min gap                                                  :CALCULATED
                    · table error  |γ_k(table) − γ_k(32 digits)|  — the table's real precision      :CALCULATED
Stage B  COUNT      Backlund exact count N(T) = θ(T)/π + 1 + S(T) at the mid-gap above every 50th zero (and above the 25,000th):
                    S(T) = arg ζ along 2+iT → ½+iT / π.  N(T) must equal k, an integer to 1e-6.  A missed or doubled zero shows as k ± 1.   :CALCULATED
                    (completeness: this is independent of how the table was made)
Stage C  RE-RUN     the wave and prime-resolution tests with all 25,000 tones instead of 3,066 (ψ_K jump across [n−½, n+½] against Λ(n)):
                    · reconstruction correlation of W_T by K = 10, 100, 1000, 3066, 25000 tones
                    · resolved fraction of primes p against p·2π/T, and the crossover p* against T/2π (573 → 3,492 was an estimate)
                    · per-integer 'ping': classify every n ≤ NMAX as prime power / not from the zero sum alone, K = 3066 and 25000
Resumable: stage A and B checkpoint to experiments/zeros25k/*.pkl; rerun the same command to continue.
"""
from __future__ import annotations
import os, sys, json, time, math, pickle, hashlib
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "zeros25k"); os.makedirs(OUT, exist_ok=True)
TABLE = os.path.join(ROOT, "..", "VAPMIP", "SemPriCon", "SpecTools", "zeros", "zeros_exact_25000.tsv")
HP = os.path.join(ROOT, "data", "zeros_25000_hp.tsv")
LOG = os.path.join(OUT, "run.log")
DPS = 32
WORKERS = int(os.environ.get("ZEROS_WORKERS", "4"))


def log(*a):
    s = time.strftime("%H:%M:%S ") + " ".join(str(x) for x in a)
    with open(LOG, "a") as f: f.write(s + "\n")
    st = dict(msg=s, t=time.time())
    with open(os.path.join(OUT, "status.json"), "w") as f: json.dump(st, f)


def load_table():
    rows = [l.split("\t") for l in open(TABLE).read().split("\n")[1:] if l.strip()]
    return [r[1] for r in rows]


# ---------------------------------------------------------------- stage A
def _refine(args):
    k, s = args
    mp.mp.dps = DPS
    x0 = mp.mpf(s)
    try:
        r = mp.findroot(mp.siegelz, (x0, x0 + mp.mpf("1e-9")), solver="secant", tol=mp.mpf("1e-28"), maxsteps=60)
        e = mp.mpf("1e-9")
        sign_ok = bool(mp.siegelz(r - e) * mp.siegelz(r + e) < 0)
        return k, mp.nstr(r, 32), float(r - x0), sign_ok, None
    except Exception as ex:
        return k, None, None, False, repr(ex)


def stage_a():
    from multiprocessing import Pool
    tab = load_table(); N = len(tab)
    ck = os.path.join(OUT, "stageA.pkl")
    done = pickle.load(open(ck, "rb")) if os.path.exists(ck) else {}
    todo = [(k, tab[k]) for k in range(N) if k not in done]
    log("stage A: todo", len(todo), "of", N, "workers", WORKERS, "dps", DPS)
    t0 = time.time()
    with Pool(WORKERS, maxtasksperchild=200) as pool:
        for b in range(0, len(todo), 400):
            blk = todo[b:b + 400]
            for k, r, d, ok, err in pool.map(_refine, blk, chunksize=8):
                done[k] = (r, d, ok, err)
            pickle.dump(done, open(ck + ".tmp", "wb")); os.replace(ck + ".tmp", ck)
            log("stage A", len(done), "/", N, "elapsed", round(time.time() - t0))
    bad = [k for k, v in done.items() if v[0] is None or not v[2]]
    mp.mp.dps = DPS
    g = [mp.mpf(done[k][0]) if done[k][0] else None for k in range(N)]
    inc = all(g[i] is not None and g[i + 1] is not None and g[i + 1] > g[i] for i in range(N - 1))
    gaps = [float(g[i + 1] - g[i]) for i in range(N - 1) if g[i] is not None and g[i + 1] is not None]
    dif = np.array([abs(done[k][1]) for k in range(N) if done[k][1] is not None])
    with open(HP, "w") as f:
        f.write("k\tgamma_32digits  (recreated in-house; mpmath siegelz root, dps 32; see experiments/zeros25k_rerun.py)\n")
        for k in range(N): f.write(f"{k + 1}\t{done[k][0]}\n")
    sha = hashlib.sha256(open(HP, "rb").read()).hexdigest()
    res = dict(N=N, failures=len(bad), first_failures=bad[:10], strictly_increasing=bool(inc), min_gap=min(gaps), max_gap=max(gaps),
               table_err_max=float(dif.max()), table_err_median=float(np.median(dif)), table_err_at=int(np.argmax(dif)) + 1,
               table_err_gt_1e_15=int((dif > 1e-15).sum()), table_err_gt_1e_12=int((dif > 1e-12).sum()),
               gamma_first=done[0][0], gamma_last=done[N - 1][0], hp_sha256=sha)
    log("stage A result", json.dumps(res))
    return res


# ---------------------------------------------------------------- stage B
def _count(args):
    k, a, b = args                                   # expect N(T) = k at T between zero k and k+1 (1-based k)
    mp.mp.dps = 25
    T = (mp.mpf(a) + mp.mpf(b)) / 2
    th = mp.siegeltheta(T)
    z_prev = mp.zeta(mp.mpc(2, T)); tot = mp.arg(z_prev)
    sig = mp.mpf(2); h = mp.mpf("0.05")
    def step(s0, s1, z0, depth=0):
        z1 = mp.zeta(mp.mpc(s1, T))
        d = mp.arg(z1 / z0)
        if abs(d) > 1.0 and depth < 12:
            sm = (s0 + s1) / 2
            da, za = step(s0, sm, z0, depth + 1)
            db, zb = step(sm, s1, za, depth + 1)
            return da + db, zb
        return d, z1
    while sig > mp.mpf("0.5") + mp.mpf("1e-12"):
        nxt = max(mp.mpf("0.5"), sig - h)
        d, z_prev = step(sig, nxt, z_prev)
        tot += d; sig = nxt
    S = tot / mp.pi
    N = th / mp.pi + 1 + S
    return k, float(T), float(N), float(abs(N - mp.nint(N)))


def stage_b():
    from multiprocessing import Pool
    tab = load_table(); N = len(tab)
    ck = os.path.join(OUT, "stageB.pkl")
    done = pickle.load(open(ck, "rb")) if os.path.exists(ck) else {}
    mp.mp.dps = DPS
    gap_last = mp.zetazero(N + 1).imag                 # independent: next zero above the table, from mpmath
    log("zero", N + 1, "from mpmath.zetazero:", mp.nstr(gap_last, 20))
    ks = list(range(50, N, 50)) + [N]
    jobs = []
    for k in ks:
        if k in done: continue
        a = tab[k - 1]; b = tab[k] if k < N else mp.nstr(gap_last, 25)
        jobs.append((k, a, b))
    log("stage B: heights todo", len(jobs), "of", len(ks))
    t0 = time.time()
    with Pool(max(1, WORKERS - 1), maxtasksperchild=20) as pool:
        for r in pool.imap_unordered(_count, jobs, chunksize=2):
            done[r[0]] = r[1:]
            if len(done) % 25 == 0:
                pickle.dump(done, open(ck + ".tmp", "wb")); os.replace(ck + ".tmp", ck)
                log("stage B", len(done), "/", len(ks), "elapsed", round(time.time() - t0))
    pickle.dump(done, open(ck, "wb"))
    off = {k: (done[k][1] - k) for k in ks}
    wrong = {k: v for k, v in off.items() if abs(v) > 0.01}
    frac = max(done[k][2] for k in ks)
    res = dict(heights=len(ks), count_wrong=len(wrong), wrong_first=dict(list(wrong.items())[:10]), max_integer_defect=frac,
               max_abs_N_minus_k=max(abs(v) for v in off.values()), final_N=done[N][1], next_zero_mpmath=str(mp.nstr(gap_last, 25)))
    log("stage B result", json.dumps(res))
    return res


# ---------------------------------------------------------------- stage C
def lam_table(X):
    lam = np.zeros(X + 1)
    sv = np.ones(X + 1, bool); sv[:2] = False
    for p in range(2, int(X ** .5) + 1):
        if sv[p]: sv[p * p::p] = False
    ps = np.nonzero(sv)[0]
    for p in ps:
        q = int(p)
        while q <= X: lam[q] = math.log(p); q *= int(p)
    return lam, ps


def jump_recovered(ns, g, chunk=64):
    """ψ_K(n+½) − ψ_K(n−½) from the zero sum, vectorised over integers n:  1 − 2 Re Σ [(n+½)^ρ − (n−½)^ρ]/ρ  (the −½ln(1−x⁻²) tail is O(n⁻²), kept)."""
    ns = np.asarray(ns, float); rho = 0.5 + 1j * g
    out = np.empty(len(ns))
    for a in range(0, len(ns), chunk):
        n = ns[a:a + chunk][:, None]
        hi = np.exp(rho[None, :] * np.log(n + .5)); lo = np.exp(rho[None, :] * np.log(n - .5))
        out[a:a + chunk] = 1 - 2 * np.real((hi - lo) / rho[None, :]).sum(axis=1) \
            - 0.5 * (np.log(1 - (n[:, 0] + .5) ** -2.0) - np.log(1 - (n[:, 0] - .5) ** -2.0))
    return out


def heartbeat_chunked(u, g, chunk=100):
    rho = 0.5 + 1j * g; out = np.empty(len(u))
    for a in range(0, len(u), chunk):
        out[a:a + chunk] = -2 * np.real(np.exp(1j * np.outer(u[a:a + chunk], g)) / rho).sum(axis=1)
    return out


def stage_c(NMAX=int(os.environ.get("ZEROS_NMAX", "8000"))):
    g_all = np.array([float(l.split("\t")[1]) for l in open(HP).read().split("\n")[1:] if l.strip()])
    N = len(g_all); T = g_all[-1]; K3 = 3066
    res = dict(T=T, T_over_2pi=T / (2 * math.pi), T3066_over_2pi=g_all[K3 - 1] / (2 * math.pi))
    # C1 reconstruction of W_T = (ψ(x) − x)/√x (u = ln x between jumps) from K tones
    X = 9_000_000
    lam, ps = lam_table(X); psi = np.cumsum(lam)
    u = np.linspace(6, 15.9, 6000)                                # 6000 pts (30000 in the atlas): memory-lean
    x = np.exp(u); wt = (psi[np.floor(x).astype(int)] - x) / np.sqrt(x)
    rec = {}
    for K in (10, 100, 1000, K3, N):
        t0 = time.time()
        h = heartbeat_chunked(u, g_all[:K]) - math.log(2 * math.pi) * np.exp(-u / 2)
        rec[K] = dict(corr=float(np.corrcoef(wt, h)[0, 1]), rms=float(np.sqrt(np.mean((wt - h) ** 2))), secs=round(time.time() - t0, 1))
        log("stage C1 K", K, rec[K])
    res["reconstruction"] = rec
    # C2 + C3 per-integer ping, n = 2..NMAX
    ns = np.arange(2, NMAX + 1)
    true = lam[ns]
    ping = {}
    for K in (K3, N):
        t0 = time.time(); rj = jump_recovered(ns, g_all[:K])
        pp = true > 0                                         # prime powers (true jump ln p ≥ ln 2)
        thr = 0.5 * math.log(2)                               # halfway between 0 and the smallest jump ln 2
        pred = rj > thr
        err_pp = float(np.mean(~pred[pp])); err_np = float(np.mean(pred[~pp]))
        primes = np.isin(ns, ps)
        fr = rj[primes] / np.log(ns[primes])
        ping[K] = dict(secs=round(time.time() - t0, 1), miss_rate_primepowers=err_pp, false_alarm_rate_others=err_np,
                       rms_err_composite=float(np.sqrt(np.mean(rj[~pp] ** 2))), rms_err_prime=float(np.sqrt(np.mean((rj[pp] - true[pp]) ** 2))))
        # blocks by n: first n at which the miss rate over a window of 100 primes exceeds 0.5
        pn = ns[primes]; hit = (fr > 0.5).astype(float)
        w = 50; roll = np.convolve(hit, np.ones(w) / w, mode="valid")
        below = np.nonzero(roll < 0.5)[0]
        ping[K]["p_star_resolved_half"] = int(pn[below[0] + w // 2]) if len(below) else None
        ping[K]["resolved_by_block"] = {f"{a}-{b}": float(np.mean(hit[(pn >= a) & (pn < b)])) for a, b in
                                        ((2, 100), (100, 300), (300, 600), (600, 1000), (1000, 2000), (2000, 3500), (3500, 5000), (5000, 8001))
                                        if np.any((pn >= a) & (pn < b))}
        np.save(os.path.join(OUT, f"jump_K{K}.npy"), rj)
        log("stage C2 K", K, json.dumps(ping[K]))
    res["ping"] = ping; res["NMAX"] = NMAX
    json.dump(res, open(os.path.join(OUT, "stageC.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    stages = sys.argv[1:] or ["A", "B", "C"]
    allres = {}
    p = os.path.join(OUT, "results.json")
    if os.path.exists(p): allres = json.load(open(p))
    log("START stages", stages, "pid", os.getpid())
    for s in stages:
        allres[s] = dict(A=stage_a, B=stage_b, C=stage_c)[s]()
        json.dump(allres, open(p, "w"), indent=1)
    log("DONE")
