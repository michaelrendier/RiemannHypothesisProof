"""Zero-phase windings against Archimedean (polygonal) windings (rhp.zerophase).  Writes zerophase_probe.json."""
import json, sys, pathlib, math, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import zerophase as Z

t0 = time.time(); R = {}; X = 50_000
w, lx = Z.weights(X); print("ln x =", round(lx, 2), flush=True)
g = np.array([float(l.split()[-1]) for l in open("data/zeros_25000_hp.tsv") if l[0].isdigit()])

# 1  the lifts: bearing of each zero-tone under the two lifts
u = np.linspace(0.5, 6, 4000); R["bearing"] = {}
for k in (1, 2, 3, 10, 100):
    gk = g[k - 1]; bs = Z.bearing(Z.lift_stereo(Z.zero_spiral(u, gk))); bl = Z.bearing(Z.lift_linear(u, gk * u, 6.0))
    R["bearing"][k] = dict(gamma=float(gk), pitch_deg=Z.pitch_deg(gk), stereo_mean_deg=float(np.degrees(bs.mean())), stereo_std_deg=float(np.degrees(bs.std())),
                           theory_deg=math.degrees(math.atan(2 * gk)), linear_std_deg=float(np.degrees(bl.std())))
    print(f"zero {k}: γ = {gk:.3f}  pitch {R['bearing'][k]['pitch_deg']:.3f}°  stereo bearing {R['bearing'][k]['stereo_mean_deg']:.4f}° ± {R['bearing'][k]['stereo_std_deg']:.1e} (theory {R['bearing'][k]['theory_deg']:.4f}°)  linear lift ± {R['bearing'][k]['linear_std_deg']:.2f}°", flush=True)

# 2  zero-phase scan, t in [10, 100]
ts = np.arange(10.0, 100.0, 0.02); Az = Z.scan_zero(ts, w, lx)
zin = g[(g > 10) & (g < 100)]; pk_t, pk_h = Z.peaks(Az, ts, 0.3)
match = [float(np.min(np.abs(pk_t - z))) for z in zin]
mid = 0.5 * (zin[1:] + zin[:-1]); Amid = np.interp(mid, ts, Az)
R["zero_scan"] = dict(x=X, zeros_in_range=len(zin), peaks_above_0p3=len(pk_t), peaks_within_0p05_of_a_zero=int(np.sum(np.array(match) < 0.05)), max_offset=float(max(match)),
                      peak_height_mean=float(np.interp(zin, ts, Az).mean()), peak_height_min=float(np.interp(zin, ts, Az).min()), midpoint_median=float(np.median(Amid)), midpoint_max=float(Amid.max()), median_all=float(np.median(Az)))
print("zero scan:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in R["zero_scan"].items()}, flush=True)

# 3  Archimedean scan, h in [0.02, 8]
hs = np.arange(0.02, 8.0, 0.001); Aa = Z.scan_arch(hs, w, lx)
ph, pa = Z.peaks(Aa, hs, 0.0); top = np.argsort(pa)[::-1][:12]
R["arch_scan"] = dict(x=X, max=float(Aa.max()), median=float(np.median(Aa)), argmax=float(hs[Aa.argmax()]), top_peaks=[(float(ph[i]), float(pa[i])) for i in top])
R["polygons"] = {}
for m in (3, 4, 5, 6, 8):
    h = Z.polygon_h(m); vals = {j: float(np.interp(j * h, hs, Aa)) for j in (1, 2, 3, 4) if j * h < 8}
    R["polygons"][m] = dict(h=h, A_at_multiples=vals)
print("arch scan: max", round(R["arch_scan"]["max"], 3), "at h =", round(R["arch_scan"]["argmax"], 3), " median", round(R["arch_scan"]["median"], 3), flush=True)
print("   top peaks (h, A):", [(round(a, 3), round(b, 3)) for a, b in R["arch_scan"]["top_peaks"][:8]], flush=True)
for m, d in R["polygons"].items(): print(f"   {m}-gon h = {d['h']:.4f}: A at j·h =", {j: round(v, 3) for j, v in d["A_at_multiples"].items()}, flush=True)
np.savez_compressed(pathlib.Path(__file__).with_name("zerophase_scans.npz"), ts=ts, Az=Az, hs=hs, Aa=Aa)
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1, default=float); print(f"done {time.time()-t0:.0f}s")
