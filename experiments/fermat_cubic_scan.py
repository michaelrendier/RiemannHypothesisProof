"""Zero-phase scan on the Fermat cubic: zeros of L(E,s) from Fermat-side counts alone against the zeros from the functional equation (rhp.fermat_cubic).  Writes fermat_cubic_scan.json."""
import json, math, sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import fermat_cubic as FC
from rhp.sedenion_slots import primes_upto
t0 = time.time(); R = {}
X = 2_000_000; a = FC.ap_table(X); P = primes_upto(X); an = FC.an_table(120, a)
# independent zeros from the approximate functional equation
TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0
zs, ts, vals = FC.find_zeros(an, 0.2, TMAX, 0.05, nterms=110, dps=50); R["afe_zeros"] = zs.tolist()
print(f"AFE zeros of L(E,1+it), t in (0.2, {TMAX}): {len(zs)}  first: {np.round(zs[:10], 6).tolist()}   [{time.time()-t0:.0f}s]", flush=True)
# the scan from a_p only
ns, ls = FC.lam_von_mangoldt(P, a, X); grid = np.arange(0.05, TMAX + 1, 0.01); A = FC.scan(grid, ns, ls)
lx = math.log(X)
from rhp import zerophase as ZP
pk_t, pk_h = ZP.peaks(A, grid, 0.35)
match = [(float(z), float(pk_t[np.argmin(np.abs(pk_t - z))]), float(abs(pk_t[np.argmin(np.abs(pk_t - z))] - z))) for z in zs]
hz = np.interp(zs, grid, A); mid = 0.5 * (zs[1:] + zs[:-1]); hm = np.interp(mid, grid, A)
R["scan"] = dict(x=X, ln_x=lx, n_prime_powers=int(len(ns)), peaks_above_0p35=int(len(pk_t)), zeros=len(zs), matched_within_0p03=int(sum(m[2] < 0.03 for m in match)), max_offset=max(m[2] for m in match),
                 height_mean=float(hz.mean()), height_min=float(hz.min()), height_max=float(hz.max()), midpoint_median=float(np.median(hm)), midpoint_max=float(hm.max()), median_all=float(np.median(A)),
                 unmatched_peaks=[float(p) for p in pk_t if np.min(np.abs(zs - p)) > 0.05])
print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in R["scan"].items() if k != "unmatched_peaks"}, "unmatched peaks:", np.round(R["scan"]["unmatched_peaks"], 3).tolist(), flush=True)
for z, p, d in match[:12]: print(f"   zero {z:9.5f}   scan peak {p:9.5f}   offset {d:.4f}")
# control: the same scan with shuffled a_p signs (destroys the arithmetic, keeps the sizes)
rng = np.random.default_rng(0); sg = rng.choice([-1.0, 1.0], len(ls)); Ac = FC.scan(grid, ns, ls * sg); R["control_shuffled_signs"] = dict(max=float(Ac.max()), median=float(np.median(Ac)), peaks_above_0p35=int(len(ZP.peaks(Ac, grid, 0.35)[0])))
print("control (random signs):", R["control_shuffled_signs"], flush=True)
np.savez_compressed(pathlib.Path(__file__).with_name("fermat_cubic_scan.npz"), grid=grid, A=A, zs=zs)
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1); print(f"done {time.time()-t0:.0f}s")
