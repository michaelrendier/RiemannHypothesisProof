# SPDX-License-Identifier: GPL-3.0-only
"""nodal_extract_run.py — the nodal-line extractor on (R1) the prime wave built from K zero-tones and (R2) ξ on the strip.   :CALCULATED
Writes experiments/nodal_extract/{r1.json, r2.json, xi_grid.npz}.   Run: python3 experiments/nodal_extract_run.py [R1] [R2]"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, ROOT)
OUT = os.path.join(HERE, "nodal_extract"); os.makedirs(OUT, exist_ok=True)
from rhp import nodal_extract as NX
G = np.array([float(l.split("\t")[1]) for l in open(os.path.join(ROOT, "data", "zeros_25000_hp.tsv")).read().split("\n")[1:] if l.strip()])
def log(*a):
    with open(os.path.join(OUT, "run.log"), "a") as f: f.write(time.strftime("%H:%M:%S ") + " ".join(str(x) for x in a) + "\n")

def r1():
    psi, _ = NX.psi_table(100300); res = {}
    for lo, hi in ((1000, 1500), (3000, 3500), (10000, 10500), (30000, 30500), (100000, 100200)):
        ex, dr, kd = NX.exact_prime_nodes(lo, hi, psi)
        u0, u1 = math.log(lo), math.log(hi); u = np.arange(u0, u1, 1e-5); xm = math.sqrt(lo * hi)
        row = dict(window=[lo, hi], exact_nodes=int(len(ex)), exact_jump=int(sum(k == "jump" for k in kd)), exact_slide=int(sum(k == "slide" for k in kd)), per_K={})
        exu = np.log(ex)
        for K in (10, 100, 1000, 3066, 25000):
            t0 = time.time(); w = NX.zero_tone_wave(u, G[:K]); pos, d = NX.sign_nodes(u, w); x = np.exp(pos)
            tol = 0.25 / xm
            f_e, f_f = NX.match_nodes(exu, pos, tol)
            chance = float(np.mean([NX.match_nodes(np.log(ex + sh), pos, tol)[0] for sh in (1.3, 2.7, 4.1, 5.9, 7.7)]))
            row["per_K"][K] = dict(nodes=int(len(pos)), rice_predicted=NX.rice_rate(G[:K]) * (u1 - u0), exact_matched=f_e, found_matched=f_f, chance=chance, secs=round(time.time() - t0, 1))
            log("R1", lo, hi, K, row["per_K"][K])
        res[f"{lo}-{hi}"] = row
    json.dump(res, open(os.path.join(OUT, "r1.json"), "w"), indent=1)

def _col(args):
    sig, ts = args
    from rhp import nodal as N
    from mpmath import mpc
    return [complex(N.xi(mpc(sig, t))) for t in ts]

def r2(tmax=80.0, dt=0.025, ds=0.01, smax=1.5):
    from multiprocessing import Pool
    sig = np.round(np.arange(0.5, smax + 1e-9, ds), 10); ts = np.arange(0.5, tmax + 1e-9, dt)
    path = os.path.join(OUT, "xi_grid.npz")
    if os.path.exists(path): Zh = np.load(path)["Z"]
    else:
        t0 = time.time()
        with Pool(4) as pool: cols = pool.map(_col, [(float(s), ts) for s in sig], chunksize=2)
        Zh = np.array(cols).T                                        # rows t, columns σ ≥ ½
        np.savez_compressed(path, Z=Zh, sig=sig, ts=ts); log("R2 grid", Zh.shape, round(time.time() - t0), "s")
    # Re ξ = 0 on both halves (mirror ξ(σ+it) = conj ξ(1−σ+it)), Im ξ = 0 on σ > ½ only (Im ≡ 0 on the equator)
    sigfull = np.concatenate([1 - sig[:0:-1], sig]); Zfull = np.concatenate([np.conj(Zh[:, :0:-1]), Zh], axis=1)
    Rl, Rc = NX.zero_set_2d(sigfull, ts, Zfull.real)
    sig_up = sig[1:]; Il, Ic = NX.zero_set_2d(sig_up, ts, Zh[:, 1:].imag)
    pts_off = NX.seg_intersections(Rl, Il)
    # rungs: where each Re ξ = 0 line meets the equator = zeros on the line
    eq = np.array([[0.5, ts[0]], [0.5, ts[-1]]]); rung_pts = NX.seg_intersections(Rl, [eq])
    from mpmath import zetazero
    known = G[G < tmax]
    found = np.sort(rung_pts[:, 1]); 
    err = [float(np.min(np.abs(known - f))) for f in found]
    npos, nneg, lab = NX.sign_domains(Zfull.real)
    ends = {"back_to_equator_side": 0, "reaches_sigma_max": 0, "reaches_t_boundary": 0, "closed": 0}
    for L, c in zip(Rl, Rc):
        if c: ends["closed"] += 1; continue
        for p in (L[0], L[-1]):
            pass
        e = [L[0], L[-1]]
        kinds = set()
        for p in e:
            if abs(p[0] - sigfull[0]) < 2 * ds or abs(p[0] - sigfull[-1]) < 2 * ds: kinds.add("reaches_sigma_max")
            elif abs(p[1] - ts[0]) < 2 * dt or abs(p[1] - ts[-1]) < 2 * dt: kinds.add("reaches_t_boundary")
            else: kinds.add("back_to_equator_side")
        for k in kinds: ends[k] += 1
    res = dict(t_max=tmax, grid=[len(ts), len(sigfull)], zeros_known_below_tmax=int(len(known)), rungs_found=int(len(found)),
               rung_max_err_vs_table=float(max(err)) if err else None, rung_median_err=float(np.median(err)) if err else None,
               unmatched_known=int(sum(1 for k in known if np.min(np.abs(found - k)) > 0.05)) if len(found) else int(len(known)),
               re_lines=len(Rl), re_closed=int(sum(Rc)), re_length=float(sum(NX.polyline_length(l) for l in Rl)),
               im_lines_offline=len(Il), im_closed=int(sum(Ic)), offline_crossings=int(len(pts_off)), offline_points=pts_off.tolist()[:10],
               re_domains_pos=npos, re_domains_neg=nneg, line_ends=ends,
               im_offline_ends=[[float(L[0][0]), float(L[0][1]), float(L[-1][0]), float(L[-1][1])] for L in Il][:40])
    json.dump(res, open(os.path.join(OUT, "r2.json"), "w"), indent=1); log("R2", json.dumps({k: v for k, v in res.items() if k != "im_offline_ends"}))

if __name__ == "__main__":
    for st in (sys.argv[1:] or ["R1", "R2"]): {"R1": r1, "R2": r2}[st]()
    log("DONE")
