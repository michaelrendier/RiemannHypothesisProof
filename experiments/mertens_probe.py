"""The Möbius / Mertens wave from the 25,000 zeros (rhp.mertens_wave).  Writes mertens_probe.json."""
import json, sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import mertens_wave as W, nodal_extract as NE

t0 = time.time(); R = {}; rng = np.random.default_rng(5)
g, zp = W.load_zp(); c = W.residues(g, zp); K = len(g); print("zeros", K, flush=True)

# 1  the slope table: M tone = ψ tone × (−1/ζ'(ρ))
az = np.abs(zp); gap = np.minimum(np.r_[np.inf, np.diff(g)], np.r_[np.diff(g), np.inf])
R["slope"] = dict(mean_abs_zp=float(az.mean()), min_abs_zp=float(az.min()), argmin_k=int(az.argmin() + 1), max_amplitude_ratio=float((1 / az).max()),
                  corr_log_abs_zp_log_gap=float(np.corrcoef(np.log(az), np.log(gap))[0, 1]),
                  partial_variance={str(k): float(2 * np.sum(1 / np.abs((0.5 + 1j * g[:k]) * zp[:k]) ** 2)) for k in (100, 1000, 10000, K)})
print("slope:", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in R["slope"].items() if k != "partial_variance"}, R["slope"]["partial_variance"], flush=True)

# 2  exact waves to 10^7
X = 10_000_000; mu = W.mobius_table(X); M = np.cumsum(mu.astype(np.int64)); n = np.arange(X + 1)
from rhp.nodal_extract import psi_table
psi, _ = psi_table(X)
decs = [(10 ** a, 10 ** (a + 1)) for a in range(2, 7)]
def sgn_changes(a):
    s = np.sign(a); s = s[s != 0]; return int((s[1:] != s[:-1]).sum())
R["exact"] = {"max_M_over_sqrt": float(np.max(np.abs(M[2:]) / np.sqrt(n[2:]))), "decades": {}}
for lo, hi in decs:
    sl = slice(lo, hi); wm = M[sl] / np.sqrt(n[sl]); wp = (psi[sl] - n[sl]) / np.sqrt(n[sl])
    a = psi[lo - 1:hi - 1] - n[lo:hi]; b = psi[lo:hi] - n[lo:hi]
    psi_nodes = int(((a * b) < 0).sum() + ((n[lo:hi] < psi[lo:hi]) & (psi[lo:hi] < n[lo:hi] + 1)).sum())
    R["exact"]["decades"][f"{lo:.0e}-{hi:.0e}"] = dict(rms_M=float(np.sqrt(np.mean(wm ** 2))), rms_psi=float(np.sqrt(np.mean(wp ** 2))), corr_M_psi=float(np.corrcoef(wm, wp)[0, 1]),
                                                          M_sign_changes=sgn_changes(M[sl]), psi_crossings=psi_nodes)
    print(f"  {lo:.0e}-{hi:.0e}: rms W_M {R['exact']['decades'][f'{lo:.0e}-{hi:.0e}']['rms_M']:.3f}  rms W_ψ {R['exact']['decades'][f'{lo:.0e}-{hi:.0e}']['rms_psi']:.3f}  corr {R['exact']['decades'][f'{lo:.0e}-{hi:.0e}']['corr_M_psi']:+.3f}  M sign changes {R['exact']['decades'][f'{lo:.0e}-{hi:.0e}']['M_sign_changes']}  ψ crossings {psi_nodes}", flush=True)
print("  max |M|/√x to 1e7:", round(R["exact"]["max_M_over_sqrt"], 4), flush=True)

# 3  reconstruction error of M from K zero-tones
R["recon"] = {}
for lo, hi in ((100, 300), (1000, 3000), (3000, 10000), (10000, 20000)):
    x = rng.uniform(lo, hi, 300); x = np.floor(x) + 0.5; ex = M[np.floor(x).astype(int)]
    R["recon"][f"{lo}-{hi}"] = {Kk: float(np.sqrt(np.mean(((W.m_reconstruct(x, g, c, Kk) - ex) / np.sqrt(x)) ** 2))) for Kk in (300, 1000, 3000, 10000, K)}
    print("  recon rms error/√x", f"{lo}-{hi}", {k: round(v, 3) for k, v in R["recon"][f"{lo}-{hi}"].items()}, flush=True)

# 4  μ ping
R["ping"] = {}
for lo, hi in ((2, 200), (200, 1000), (1000, 3000), (3000, 10000), (10000, 20000)):
    ns = np.sort(rng.choice(np.arange(lo, hi), min(500, hi - lo), replace=False)); row = {}
    for Kk in (300, 1000, 3000, 10000, K):
        cl = W.classify(W.mu_ping(ns, g, c, Kk)); tr = mu[ns].astype(int)
        nz = tr != 0; row[Kk] = dict(acc3=float(np.mean(cl == tr)), squarefree_acc=float(np.mean((cl != 0) == nz)), sign_acc=float(np.mean(cl[nz] == tr[nz])))
    R["ping"][f"{lo}-{hi}"] = row; print("  μ-ping acc (3-class)", f"{lo}-{hi}", {k: round(v["acc3"], 3) for k, v in row.items()}, flush=True)

# 5  the convolution identity: Σ_{n ≤ x} M_K(x/n) = 1
R["conv"] = {str(x): {Kk: W.convolution_residual(x, g, c, Kk) for Kk in (1000, 3000, K)} for x in (50, 200, 1000, 3000)}
print("  convolution residual Σ M_K(x/n) − 1:", {x: {k: round(v, 4) for k, v in d.items()} for x, d in R["conv"].items()}, flush=True)
# 6  two questions from the same tones: 'is n a hole' (ψ ping: prime power) and 'what is n's factor content' (μ ping)
from rhp import wave
_, lam = psi_table(20000); R["two_questions"] = {}
for lo, hi in ((2, 200), (200, 1000), (1000, 3000), (3000, 10000), (10000, 20000)):
    ns = np.sort(rng.choice(np.arange(lo, hi), min(400, hi - lo), replace=False)); row = {}
    isp = lam[ns] > 0; nz = mu[ns] != 0
    for Kk in (3000, K):
        hole = np.array([wave.jump_fraction(int(m), g, Kk) for m in ns]) > 0.5
        cl = W.classify(W.mu_ping(ns, g, c, Kk)); tr = mu[ns].astype(int)
        row[Kk] = dict(hole_balanced_acc=float(0.5 * (np.mean(hole[isp]) + np.mean(~hole[~isp]))), mu_acc3=float(np.mean(cl == tr)),
                       composite_squarefree_acc=float(np.mean((cl != 0)[nz & ~isp] == True)) if (nz & ~isp).any() else None)
    R["two_questions"][f"{lo}-{hi}"] = row; print("  hole vs factor-content", f"{lo}-{hi}", {k: (round(v["hole_balanced_acc"], 3), round(v["mu_acc3"], 3)) for k, v in row.items()}, flush=True)
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1, default=float)
print(f"done {time.time()-t0:.0f}s")
