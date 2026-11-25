"""Curvature and torsion of both paths; the rate-plane Mach number M = |Re lambda| / |Im lambda|; the speed ceilings.  Writes curvature_and_mach.json.
Expectations written before the run:
  K1  zeta in the value plane: exact curvature kappa = Im(conj(a) b)/|a|^3 with a = Z' - i theta' Z, b = Z'' - 2 i theta' Z' - (theta'^2 + i theta'') Z; at a zero kappa = 2 theta'/|Z'| exactly, so the circle radius is |Z'|/(2 theta') (compare with the circle fits).
  K2  Aizawa: space-curve curvature and torsion finite and bounded (a smooth attractor), in both the Cartesian and the UNS torus charts.
  M1  Aizawa is subsonic (M <= max|z-b|/d = 0.33).  zeta: in the circle model M = |tan(theta - theta_0)| per lobe, so about half the phase of every lobe is supersonic and M diverges at the zeros.
  M2  lowering Aizawa's rotation rate d until max|z - b| / d reaches 1 (d ~ 1.1) changes its character (reported either way: Lyapunov exponent, radius range, ring thickness).
  S1  speed through the origin |Z'(gamma)| has no fixed ceiling: its mean square grows like (ln T)^3 (Gonek); the largest value seen grows slowly; the rigorous-style bound is far above it."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.zeta_uns_path import maths as ZM
from ValaQuenta.modules.aizawa_attractor import maths as AZ
R = {}; g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()])
ts = ZM.sample_by_phase(10.0, 25000, 16, jitter=True); h = 1e-3
Z0 = ZM.Z_rs(ts); Zp = (ZM.Z_rs(ts + h) - ZM.Z_rs(ts - h)) / (2 * h); Zpp = (ZM.Z_rs(ts + h) - 2 * Z0 + ZM.Z_rs(ts - h)) / h ** 2
thp = ZM.theta_prime(ts); thpp = 1 / (2 * ts) + 1 / (24 * ts ** 3) + 7 / (480 * ts ** 5)
a = Zp - 1j * thp * Z0; b = Zpp - 2j * thp * Zp - (thp ** 2 + 1j * thpp) * Z0; kap = np.imag(np.conj(a) * b) / np.abs(a) ** 3; speed = np.abs(a)
dt = np.gradient(ts); w = dt / dt.sum()
# K1: curvature at zeros and the circle radius
zpg = np.array([(float(ZM.Z_rs(x + 1e-4)[0]) - float(ZM.Z_rs(x - 1e-4)[0])) / 2e-4 for x in g]); thg = ZM.theta_prime(g); rad_pred = np.abs(zpg) / (2 * thg)
R['K1'] = dict(median_abs_kappa_path=float(np.median(np.abs(kap))), p90_abs_kappa=float(np.percentile(np.abs(kap), 90)), share_kappa_positive=float(np.mean(kap > 0)), circle_radius_at_zero_Zp_over_2thetap_median=float(np.median(rad_pred)), radius_p5_p95=[float(np.percentile(rad_pred, 5)), float(np.percentile(rad_pred, 95))])
# compare with the lobe circle fits (|c|)
rr = []; cc = []
for i in range(300, 24990, 25):
    mk = (ts > g[i]) & (ts < g[i + 1])
    if mk.sum() < 10: continue
    z = ZM.zeta_line(ts[mk]); X = np.stack([z.real, z.imag], 1); yv = np.abs(z) ** 2; c_, *_ = np.linalg.lstsq(2 * X, yv, rcond=None); rr.append(np.hypot(*c_)); cc.append(rad_pred[i])
R['K1']['circle_fit_radius_vs_Zp_over_2thetap_corr'] = float(np.corrcoef(rr, cc)[0, 1]); R['K1']['median_ratio_fit_over_pred'] = float(np.median(np.array(rr) / np.array(cc)))
# K2: Aizawa space-curve curvature and torsion
o = AZ.orbit(200000, 3000, 0.0045); hh = 0.0045
def kt(P):
    v = np.gradient(P, hh, axis=0); ac = np.gradient(v, hh, axis=0); jk = np.gradient(ac, hh, axis=0); cr = np.cross(v, ac); k = np.linalg.norm(cr, axis=1) / np.linalg.norm(v, axis=1) ** 3; t = np.einsum('ij,ij->i', cr, jk) / np.maximum(np.linalg.norm(cr, axis=1) ** 2, 1e-12); return k[5:-5], t[5:-5]
k1, t1 = kt(o); r = np.hypot(o[:, 0], o[:, 1]); P2 = np.stack([(np.log(r) + 3.0) * 0.55 * np.cos(np.arctan2(o[:, 1], o[:, 0])), (o[:, 2] - 0.7) * 0.9, (np.log(r) + 3.0) * 0.55 * np.sin(np.arctan2(o[:, 1], o[:, 0]))], 1); k2, t2 = kt(P2)
q = lambda v: [float(np.percentile(v, 5)), float(np.median(v)), float(np.percentile(v, 95)), float(np.max(np.abs(v)))]
R['K2'] = dict(cartesian_kappa_p5_med_p95_max=q(k1), cartesian_torsion_p5_med_p95_maxabs=q(t1), uns_kappa_p5_med_p95_max=q(k2), uns_torsion_p5_med_p95_maxabs=q(t2), speed_cartesian_max=float(np.linalg.norm(np.gradient(o, hh, axis=0), axis=1).max()))
# M1: Mach number
rate = Zp / Z0; M = np.abs(rate) / thp; R['M1_zeta'] = dict(share_supersonic_phase_uniform=float(np.mean(M > 1)), share_supersonic_t_weighted=float(np.sum(w[M > 1])), M_median=float(np.median(M)), M_p90=float(np.percentile(M, 90)), share_M_gt_10=float(np.mean(M > 10)))
# tan model per lobe: rate/theta' = -tan(theta - theta_peak)
res = []
for i in range(300, 24990, 40):
    mk = (ts > g[i]) & (ts < g[i + 1])
    if mk.sum() < 12: continue
    zz = np.abs(Z0[mk]); j = int(np.argmax(zz)); thv = ZM.theta(ts[mk]); psi = thv - thv[j]; pred = -np.tan(psi); obs = rate[mk] / thp[mk]; sel = np.abs(pred) < 4; res.append(np.corrcoef(pred[sel], obs[sel])[0, 1] if sel.sum() > 4 else np.nan)
res = np.array(res); R['M1_zeta']['tan_model_corr_median'] = float(np.nanmedian(res)); R['M1_zeta']['tan_model_corr_p10'] = float(np.nanpercentile(res, 10))
ra = o[:, 2] - AZ.STANDARD['b']; R['M1_aizawa'] = dict(max_abs_rate=float(np.abs(ra).max()), d=AZ.STANDARD['d'], M_max=float(np.abs(ra).max() / AZ.STANDARD['d']), share_supersonic=float(np.mean(np.abs(ra) > AZ.STANDARD['d'])))
# M2: scan d
rows = []
for dv in (0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0):
    pp = dict(AZ.STANDARD, d=dv); oo = AZ.orbit(40000, 4000, 0.0045, pp); zz_ = oo[:, 2] - pp['b']; rr_ = np.hypot(oo[:, 0], oo[:, 1]); lam = AZ.lyapunov_max(steps=30000, p=pp)
    rows.append(dict(d=dv, M_max=float(np.abs(zz_).max() / dv), share_M_gt_1=float(np.mean(np.abs(zz_) > dv)), lyapunov=float(lam), r_min=float(rr_.min()), r_max=float(rr_.max()), ring_thickness_ln_r=float(np.log(rr_.max() / rr_.min()))))
R['M2_scan_d'] = rows
# S1: ceilings
az = np.abs(zpg); rows = []
for lo, hi in ((0, 100), (100, 300), (300, 1000), (1000, 3000), (3000, 10000), (10000, 25000)):
    T_ = g[hi - 1]; ms = float(np.mean(az[lo:hi] ** 2)); gon = (1 / 12.0) * math.log(g[(lo + hi) // 2] / (2 * math.pi)) ** 3
    rows.append(dict(zeros=[lo, hi], max_speed_through_zero=float(az[lo:hi].max()), mean_square=ms, rms=math.sqrt(ms), gonek_leading=gon, ratio_ms_over_gonek=ms / gon))
R['S1'] = rows
# path speed ceiling |dzeta/dt| per decade of t
R['S1_path_speed'] = [dict(t_range=[lo, hi], max_speed=float(speed[(ts >= lo) & (ts < hi)].max())) for lo, hi in ((10, 100), (100, 1000), (1000, 10000), (10000, 22000))]
json.dump(R, open(root / 'experiments/curvature_and_mach.json', 'w'), indent=1, default=float); print(json.dumps(R, indent=0, default=float)[:6000])
