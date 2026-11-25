"""The orange-slice wedge of the zeta pathway on the value sphere, and a 'through the sphere' (interior) chart.  Writes wedge_and_interior.json.
Expectations written before the run:
  W1  the wedge is not empty: it is a depletion at LARGE |zeta| on the Re < 0 side.  Near the origin the path crosses in every direction (it goes through 0 at every zero), so the slice is cut from the
      large-|zeta| end of the sphere (the north cap in the pole view), starting at the pole and widening toward azimuth pi.
  W2  each lobe (the path between consecutive zeros) is a petal along a direction; petals pointing to -Re are shorter: their peak |zeta| is much smaller than for petals pointing to +Re, and most of them
      never leave the southern hemisphere (|zeta| < 1): they 'turn around' before the equator.
  W3  the n = 1 term alone gives the circle r = 2 cos(phi) (never Re < 0); the observed frontier max|zeta|(phi) follows a cos(phi) law plus a floor from the n >= 2 terms.
  I1  in the interior chart (u = zeta/(1+|zeta|), w = tanh(Z'/(2Z)) as the third coordinate, Aizawa's x, y, z analogue) a large share of the path lies inside the unit sphere and a large share outside:
      the path goes through the sphere."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.zeta_uns_path import maths as M
R = {}; NZ = 25000
g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()])
ts = M.sample_by_phase(10.0, NZ, 16, jitter=True); z = M.zeta_line(ts); a = np.abs(z); ph = np.angle(z); dt = np.gradient(ts); w = dt / dt.sum()
# W1: 2-D density on the sphere: azimuth (72) x latitude (36, by Zs = (m^2-1)/(m^2+1)), phase-uniform (what is drawn) and t-weighted
Zs = (a ** 2 - 1) / (a ** 2 + 1); Hd, ae, le = np.histogram2d(ph, Zs, bins=[72, 36], range=[[-math.pi, math.pi], [-1, 1]]); Hd = Hd / Hd.mean()
R['W1'] = dict(empty_cells=int((Hd == 0).sum()), cells=int(Hd.size), cells_below_0p1_of_mean=int((Hd < 0.1).sum()),
               north_cap_Re_neg_share_of_points_with_abs_gt_3=float(np.mean(z[a > 3].real < 0)), north_cap_Re_neg_share_abs_gt_1=float(np.mean(z[a > 1].real < 0)), south_Re_neg_share_abs_lt_1=float(np.mean(z[a < 1].real < 0)),
               abs_gt_3_points=int((a > 3).sum()), abs_gt_6_points=int((a > 6).sum()), abs_gt_10_points=int((a > 10).sum()))
for thr in (1, 2, 3, 6): R['W1'][f'share_of_abs_gt_{thr}_that_has_Re_lt_0'] = float(np.mean(z[a > thr].real < 0))
lat_band = [(-1, -0.5), (-0.5, 0), (0, 0.5), (0.5, 1)]
left = (np.abs(ph) > 2.4); right = (np.abs(ph) < 0.74); bands = {}
for lo, hi in lat_band: m = (Zs >= lo) & (Zs < hi); bands[f'Zs[{lo},{hi})'] = dict(left_share=float(np.sum(w[m & left]) / np.sum(w[m])), right_share=float(np.sum(w[m & right]) / np.sum(w[m])), ratio_left_over_right=float(np.sum(w[m & left]) / max(np.sum(w[m & right]), 1e-12)))
R['W1']['by_latitude_band_(left=|arg|>138deg, right=|arg|<42deg)'] = bands
# W2: lobes
pk = []; 
for i in range(len(g) - 1):
    m = (ts > g[i]) & (ts < g[i + 1])
    if m.any(): j = np.argmax(a[m]); pk.append((float(ph[m][j]), float(a[m][j]), float(z[m][j].real)))
pk = np.array(pk); phi, rad = pk[:, 0], pk[:, 1]; sec = lambda lo, hi: (np.abs(phi) >= lo) & (np.abs(phi) < hi)
rows = {}
for name, (lo, hi) in {'toward +Re (|arg|<30deg)': (0, math.radians(30)), '30-60deg': (math.radians(30), math.radians(60)), '60-90deg': (math.radians(60), math.radians(90)), '90-120deg': (math.radians(90), math.radians(120)), '120-150deg': (math.radians(120), math.radians(150)), 'toward -Re (|arg|>150deg)': (math.radians(150), 4)}.items():
    m = sec(lo, hi); rows[name] = dict(lobes=int(m.sum())) if m.sum() == 0 else dict(lobes=int(m.sum()), median_peak=float(np.median(rad[m])), p95_peak=float(np.percentile(rad[m], 95)), max_peak=float(rad[m].max()), share_peak_below_1=float(np.mean(rad[m] < 1)))
R['W2'] = dict(lobes=len(pk), by_direction=rows, lobes_with_peak_gt_1_share=float(np.mean(rad > 1)))
# W3: frontier max|zeta| per azimuth bin vs cos law
nb = 36; edges = np.linspace(-math.pi, math.pi, nb + 1); fr = np.array([a[(ph >= edges[k]) & (ph < edges[k + 1])].max() for k in range(nb)]); ctr = 0.5 * (edges[1:] + edges[:-1])
A = np.stack([np.cos(ctr), np.ones(nb)], 1); coef, *_ = np.linalg.lstsq(A, fr, rcond=None); fit = A @ coef
R['W3'] = dict(frontier_max_abs_by_azimuth_bin_deg={int(round(math.degrees(c))): round(float(v), 2) for c, v in zip(ctr, fr)}, fit_A_cos_plus_B=[float(coef[0]), float(coef[1])], R2=float(1 - np.var(fr - fit) / np.var(fr)), frontier_at_pi_over_at_0=float(fr[0] / fr[nb // 2]) if False else float(np.max(fr[:2]) / np.max(fr[nb // 2 - 1:nb // 2 + 1])))
# I1: interior chart
drv = M.driver(ts) if False else None
h = 1e-4; Zv = M.Z_rs(ts); dZ = (M.Z_rs(ts + h) - M.Z_rs(ts - h)) / (2 * h); drv = dZ / Zv; wv = np.tanh(0.5 * drv)
u = z / (1 + a); r3 = np.sqrt(np.abs(u) ** 2 + wv ** 2)
R['I1'] = dict(share_inside_unit_sphere_phase_weighted=float(np.mean(r3 < 1)), share_inside_t_weighted=float(np.sum(w[r3 < 1])), r3_median=float(np.median(r3)), r3_p05=float(np.percentile(r3, 5)), r3_p95=float(np.percentile(r3, 95)),
               share_w_positive=float(np.mean(wv > 0)), sign_flips_of_w_per_zero=float(np.sum(np.diff(np.sign(wv)) != 0) / NZ), driver_abs_median=float(np.median(np.abs(drv))))
# wedge in the interior chart: occupancy of 24 x 24 x 24
Hh, _ = np.histogramdd(np.stack([u.real, wv, u.imag], 1), bins=24, range=[[-1, 1]] * 3); R['I1']['occupied_cells_of_24^3_inside_unit_ball'] = int((Hh > 0).sum()); cx = (np.arange(24) + 0.5) / 12 - 1; X, Y, Zg = np.meshgrid(cx, cx, cx, indexing='ij'); ball = X ** 2 + Y ** 2 + Zg ** 2 < 1
R['I1']['ball_cells'] = int(ball.sum()); R['I1']['ball_cells_visited'] = int(((Hh > 0) & ball).sum()); ang3 = np.arctan2(Zg, X); left3 = ball & (np.abs(ang3) > 2.4); R['I1']['left_wedge_ball_cells'] = int(left3.sum()); R['I1']['left_wedge_ball_cells_visited'] = int(((Hh > 0) & left3).sum())
json.dump(R, open(root / 'experiments/wedge_and_interior.json', 'w'), indent=1, default=float)
print(json.dumps({k: R[k] for k in ('W1',)}, indent=0)[:1800]); print(json.dumps(R['W2'], indent=0)); print(json.dumps(R['W3'], indent=0)); print(json.dumps(R['I1'], indent=0))
