"""Hat and curvature tests: the zeros' pair correlation against the project's hat (witches_hat), the Aizawa z-well against the hat shapes, the Higgs numbers, and the lobes as circles (inversion / geodesic curvature).
Writes hat_curvature_tests.json.  Expectations written before the run:
  A  the empirical pair correlation R2(u) of the 25,000 unfolded zeros follows the sine kernel 1 - (sin pi u/pi u)^2 (Montgomery, Odlyzko): RMSE small, hole exponent 2 (R2 ~ u^2 at small u), first maximum near the fermion1d brim.
  B  the Aizawa z-well (tilted quartic) matches the hat family: RMSE to the quartic n = 4 shape is small after normalising hilltop -> brim, and larger to the fermion1d (sine-kernel) hat; the tilt sits near the saddle-node of the well (c = 0.6 against (2/3) a^(3/2) = 0.617).
  C  Higgs: lam = m_h^2/(2 v^2) = 0.129, v = 246.22 GeV, and d* = W(1)/ln 10 = 0.24631 is within 0.04% of v in TeV (a coincidence to flag, not a claim).
  D  zeta = S + e^(-2 i theta) S-bar: a lobe is approximately a circle through 0 centred on S(t_mid): the fit |zeta|^2 = 2 Re(zeta c-bar) has median R^2 above 0.8, c is close to S, and on the sphere the lobes have nearly constant geodesic curvature."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.witches_hat import maths as WH
from ValaQuenta.modules.lagrangian_sm import maths as LS
from ValaQuenta.modules.zeta_uns_path import maths as ZM
from ValaQuenta.modules.aizawa_attractor import maths as AZ
R = {}
g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()])
# A
u, R2 = WH.zeta_pair_correlation(g, umax=3.0, bins=60); sk = WH.pair_correlation('fermion1d', u)
m = u < 0.35; sl = np.polyfit(np.log(u[m]), np.log(np.maximum(R2[m], 1e-6)), 1)[0]; i1 = int(np.argmax(R2[u < 2.0])); f_, meta = WH.hat_shape('fermion1d')
R['A'] = dict(rmse_R2_vs_sine_kernel=float(np.sqrt(np.mean((R2 - sk) ** 2))), hole_log_slope_u_lt_0p35=float(sl), expected_slope=2.0, first_max_u_empirical=float(u[i1]), first_max_u_fermion1d_brim=float(meta['brim_u']),
              R2_at_u=[[round(float(a), 3), round(float(b), 3), round(float(c_), 3)] for a, b, c_ in list(zip(u, R2, sk))[::6]])
# B
a, b, c, e = (AZ.STANDARD[k] for k in 'abce'); Vz = lambda z: z ** 4 / 12 - a * z ** 2 / 2 - c * z; fps = sorted([float(f['point'][2]) for f in AZ.fixed_points()]); zh, zm = fps[1], fps[2]
x = np.linspace(0, 1.3, 131); z = zh + x * (zm - zh); Vn = (Vz(z) - Vz(zh)) / (Vz(zh) - Vz(zm)); V4 = WH.hat_shape('quartic', 4)[0](x)
Vf, mf = WH.hat_shape('fermion1d'); rm4 = float(np.sqrt(np.mean((Vn - V4) ** 2))); rmf = float(np.sqrt(np.mean((Vn - Vf(x)) ** 2)))
from scipy.optimize import brentq
c_crit = 2.0 / 3.0 * a ** 1.5
R['B'] = dict(hilltop_z=zh, brim_z=zm, rmse_to_quartic_n4=rm4, rmse_to_fermion1d=rmf, tilt_c=c, critical_tilt_saddle_node=c_crit, c_over_critical=c / c_crit, linear_sigma_model=dict(lam=1 / 3.0, mu2=a, H=c, sigma0=zm, m_sigma2=zm ** 2 - a, m_pi2_over=c / zm))
# C
hs = LS.higgs_sector(); ds = 0.5671432904097838 / math.log(10)
R['C'] = dict(m_h=hs['m_h'], v=hs['v'], lam=hs['lam'], mu2=hs['mu2'], m_h_over_v=hs['m_h'] / hs['v'], d_star=ds, v_TeV=hs['v'] / 1000.0, d_star_over_v=ds / (hs['v'] / 1000.0), percent_apart=100 * abs(ds / (hs['v'] / 1000.0) - 1))
# D
ts = ZM.sample_by_phase(10.0, 25000, 16, jitter=True); zz = ZM.zeta_line(ts); r2s = []; ratios = []; angd = []; cvk = []; cnt = 0
for i in range(300, 24990, 5):
    mk = (ts > g[i]) & (ts < g[i + 1])
    if mk.sum() < 10: continue
    w = zz[mk]; X = np.stack([w.real, w.imag], 1); y = np.abs(w) ** 2; cc, *_ = np.linalg.lstsq(2 * X, y, rcond=None); pred = 2 * X @ cc; ss = np.sum((y - y.mean()) ** 2); r2 = 1 - np.sum((y - pred) ** 2) / ss
    S = ZM.spring_sum(0.5 * (g[i] + g[i + 1])); c_ = complex(cc[0], cc[1]); r2s.append(r2); ratios.append(abs(c_) / abs(S)); angd.append(abs(np.angle(c_ / S))); cnt += 1
    P = ZM.sphere_point(w); tt = ts[mk]; dP = np.gradient(P, tt, axis=0); ddP = np.gradient(dP, tt, axis=0); kg = np.einsum('ij,ij->i', np.cross(P, dP), ddP) / np.maximum(np.linalg.norm(dP, axis=1) ** 3, 1e-9); kk = kg[2:-2]; cvk.append(float(np.std(kk) / (abs(np.median(kk)) + 1e-9)) if len(kk) > 3 else np.nan)
r2s = np.array(r2s); ratios = np.array(ratios); angd = np.array(angd); cvk = np.array(cvk)
R['D'] = dict(lobes=cnt, circle_fit_R2_median=float(np.median(r2s)), R2_p10_p90=[float(np.percentile(r2s, 10)), float(np.percentile(r2s, 90))], share_R2_gt_0p8=float(np.mean(r2s > 0.8)), share_R2_gt_0p95=float(np.mean(r2s > 0.95)),
              centre_over_S_abs_median=float(np.median(ratios)), centre_over_S_abs_IQR=[float(np.percentile(ratios, 25)), float(np.percentile(ratios, 75))], centre_vs_S_angle_median_rad=float(np.median(angd)), centre_vs_S_angle_share_lt_0p5=float(np.mean(angd < 0.5)),
              geodesic_curvature_cv_median=float(np.nanmedian(cvk)))
json.dump(R, open(root / 'experiments/hat_curvature_tests.json', 'w'), indent=1, default=float); print(json.dumps({k: v for k, v in R.items()}, indent=0, default=float)[:4500])
# small data file for the page (committed): pair correlation vs the sine-kernel hat, and the Aizawa z-well vs the hat shapes
hat = dict(u=[round(float(v), 4) for v in u], r2=[round(float(v), 4) for v in R2], sk=[round(float(v), 4) for v in sk],
           az=dict(x=[round(float(v), 4) for v in x], vn=[round(float(v), 4) for v in Vn], v4=[round(float(v), 4) for v in V4], vf=[round(float(v), 4) for v in Vf(x)]))
open(root / 'visualizations/javascript/hat_data.js', 'w').write('// generated by experiments/hat_curvature_tests.py: pair correlation of the 25,000 zeros against 1 - sinc^2, and the Aizawa z-well against the project\'s hat shapes. Not hand-edited.\nwindow.HAT = ' + json.dumps(hat) + ';\n')
