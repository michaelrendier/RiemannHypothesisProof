"""Flow tests (rhp.flow): F2 islands vs rigid points (spacing statistics), F3 the Fermat finite-tone wave (nodes and strip).  Writes flow_tests.json."""
import json, math, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import flow as FL, nodal_extract as NX
from rhp.sedenion_slots import primes_upto
rng = np.random.default_rng(1); R = {}
# F2: islands (primes) against rigid points (zeros)
g = np.array([float(l.split()[-1]) for l in open("data/zeros_25000_hp.tsv") if l[0].isdigit()]); xz = FL.unfold_zeros(g)
P = primes_upto(300_000_000).astype(float); R["F2"] = {}
for name, x in (("zeros 1..25000", xz), ("zeros 20001..25000", xz[20000:])):
    v, f01, m = FL.spacing_stats(x); R["F2"][name] = dict(spacing_var=v, frac_below_0p1=f01, nvar_L5=FL.number_variance(x, 5.0, rng=rng), nvar_L10=FL.number_variance(x, 10.0, rng=rng))
for lo, hi in ((1e6, 1e7), (1e7, 1e8), (1e8, 3e8)):
    p = P[(P >= lo) & (P < hi)]; from scipy.special import expi; xs = expi(np.log(p))                      # unfold by the smooth counting function li(p): mean spacing 1
    v, f01, m = FL.spacing_stats(xs); R["F2"][f"primes {lo:.0e}-{hi:.0e}"] = dict(spacing_var=v, frac_below_0p1=f01, nvar_L5=FL.number_variance(xs, 5.0, rng=rng), nvar_L10=FL.number_variance(xs, 10.0, rng=rng))
poisson = np.cumsum(rng.exponential(1.0, 25000)); v, f01, m = FL.spacing_stats(poisson); R["F2"]["Poisson control"] = dict(spacing_var=v, frac_below_0p1=f01, nvar_L5=FL.number_variance(poisson, 5.0, rng=rng), nvar_L10=FL.number_variance(poisson, 10.0, rng=rng))
for k, v in R["F2"].items(): print(f"{k:22s} spacing var {v['spacing_var']:.3f}  P(s<0.1) {v['frac_below_0p1']:.3f}  Σ²(5) {v['nvar_L5']:.3f}  Σ²(10) {v['nvar_L10']:.3f}", flush=True)
print("   references: Poisson var 1, Σ²(L) = L; GUE var ≈ 0.18, Σ²(5) ≈", round((math.log(2 * math.pi * 5) + 0.5772 + 1) / math.pi ** 2, 3), "Σ²(10) ≈", round((math.log(2 * math.pi * 10) + 0.5772 + 1) / math.pi ** 2, 3))
# F3: Fermat finite-tone wave
R["F3"] = {}
for n, p in ((3, 7), (5, 31), (7, 29), (4, 17)):
    th = FL.frobenius_angles(p, n); g_ = (n - 1) * (n - 2) // 2; L = math.log(p); per = 2 * math.pi / L
    key = np.round(np.mod(th, 2 * math.pi), 7); uniq, cnt = np.unique(key, return_counts=True)
    # 1-D: nodes of W_F(u), Rice
    u = np.arange(0, 4000, 0.01); w = FL.fermat_wave(u, p, th); pos, dr = NX.sign_nodes(u, w); gam = th / L
    rice = math.sqrt(np.sum(gam ** 2) / len(gam)) / math.pi * (u[-1] - u[0])
    # integer k only (what the arithmetic sees)
    k = np.arange(1, 200001, dtype=float); wk = np.array([FL.fermat_wave(k[i:i + 5000] * L, p, th) for i in range(0, len(k), 5000)]).ravel()
    # 2-D strip over 2 periods
    sig = np.linspace(0.0, 1.0, 241); ts = np.linspace(0.0, 2 * per, 1201); S, Tt = np.meshgrid(sig, ts); Z = FL.fermat_F(S + 1j * Tt, p, th)
    Rl, Rc = NX.zero_set_2d(sig, ts, Z.real); Il, Ic = NX.zero_set_2d(sig, ts, Z.imag)
    pts = NX.seg_intersections(Rl, Il); on_line = pts[np.abs(pts[:, 0] - 0.5) < 2e-3] if len(pts) else pts
    zeros_expected = np.sort(np.mod(th, 2 * math.pi) / L); zeros_expected = np.concatenate([zeros_expected, zeros_expected + per])
    ang = []
    for Lp in Rl:
        a0, a1 = Lp[:-1], Lp[1:]; s0 = a0[:, 0] - 0.5; s1 = a1[:, 0] - 0.5
        for j in np.nonzero((s0 * s1 <= 0) & (s0 != s1))[0]:
            v = a1[j] - a0[j]; a = math.degrees(math.atan2(v[1], v[0])); a = a - 180 if a > 90 else a + 180 if a < -90 else a; ang.append(a)
    R["F3"][f"n={n}, p={p}"] = dict(tones=len(th), distinct_angles=int(len(uniq)), multiplicities=sorted(cnt.tolist(), reverse=True), sqrt_2g=math.sqrt(2 * g_), rms_wave=float(np.sqrt(np.mean(w ** 2))), sqrt_g=math.sqrt(g_), nodes_u=int(len(pos)), rice_pred=float(rice), integer_k_sign_changes=int(len(NX.sign_nodes(k, wk)[0])),
                                    strip_Re_lines=len(Rl), strip_Im_lines=len(Il), crossings_on_line=int(len(on_line)), zeros_expected_in_2_periods=len(zeros_expected),
                                    rung_angle_median_abs=float(np.median(np.abs(ang))) if ang else None, rung_angle_max_abs=float(np.max(np.abs(ang))) if ang else None, rung_crossings=len(ang))
    r = R["F3"][f"n={n}, p={p}"]; print(f"n={n} p={p}: 2g={r['tones']} distinct {r['distinct_angles']} mult {r['multiplicities']}  RMS {r['rms_wave']:.3f} (√g {r['sqrt_g']:.3f})  nodes {r['nodes_u']} (Rice {r['rice_pred']:.0f})  strip: Re lines {r['strip_Re_lines']}, crossings on line {r['crossings_on_line']} (expected {r['zeros_expected_in_2_periods']}); rung angle median {r['rung_angle_median_abs']} max {r['rung_angle_max_abs']}", flush=True)
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
