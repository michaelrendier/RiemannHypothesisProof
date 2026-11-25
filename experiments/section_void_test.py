"""Test: is the empty region of the phase-portrait SECTION the shape the paths cannot reach?  Writes section_void_test.json.
Phase portrait coordinates: u = ln|radius|, a = asinh(d ln|radius| / dt).  Section: arg = 0 crossings (Aizawa: phi = 0; zeta: Im zeta crosses 0 upward with Re > 0).
Expectations written before the run:
  T1  zeta: the exact rate at the section is large (median ~ 34); the page's tanh(Z'/2Z) encoding saturates there (clipped line): an artefact.
  T2  zeta's section points lie near the hyperbola rate * |zeta| = Z' (the speed component) with Z' of order |zeta'| ~ 3: the section traces a curve, not the whole portrait.
  T3  void test: the cells the section leaves empty are NOT empty for the full path (the section is a different slice, locked near zeros): a large share of full-path points fall in the section's void.
  T4  Aizawa: the section is an unbiased sample of the (u, z) distribution (phi independent): KS p not small; the hole inside the ring is empty for the section AND the full orbit.
  T5  the forbidden region of the full zeta portrait is bounded by the hyperbola |rate| e^u = |Z'| <= B'(t) = 2 sum_{n<=N} n^(-1/2) (theta' + ln n): the largest observed |rate| e^u stays below it."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.zeta_uns_path import maths as ZM
from ValaQuenta.modules.aizawa_attractor import maths as AZ
from scipy.stats import ks_2samp
R = {}
ts = ZM.sample_by_phase(10.0, 25000, 16, jitter=True); z = ZM.zeta_line(ts); h = 1e-4; Zv = ZM.Z_rs(ts); rate = (ZM.Z_rs(ts + h) - ZM.Z_rs(ts - h)) / (2 * h) / Zv; u = np.log(np.abs(z) + 1e-300); a = np.arcsinh(rate); dt = np.gradient(ts)
sec = np.array([i for i in range(1, len(ts)) if z[i - 1].imag < 0 <= z[i].imag and z[i].real > 0]); R['T1'] = dict(section_points=len(sec), rate_median=float(np.median(rate[sec])), share_rate_gt_14p5=float(np.mean(np.abs(rate[sec]) > 14.5)), share_full_path_gt_14p5=float(np.mean(np.abs(rate) > 14.5)))
Zp = rate * np.abs(z); R['T2'] = dict(section_Zprime_median=float(np.median(np.abs(Zp[sec]))), section_Zprime_p5_p95=[float(np.percentile(np.abs(Zp[sec]), 5)), float(np.percentile(np.abs(Zp[sec]), 95))], full_Zprime_median=float(np.median(np.abs(Zp))), full_Zprime_p95=float(np.percentile(np.abs(Zp), 95)),
              section_rate_sign_positive=float(np.mean(rate[sec] > 0)), corr_u_a_section=float(np.corrcoef(u[sec], a[sec])[0, 1]), corr_u_a_full=float(np.corrcoef(u, a)[0, 1]))
# grids
U0, U1, A0, A1, G = -6.0, 3.3, -7.0, 7.0, 60; cell = lambda uu, aa: (np.clip(((uu - U0) / (U1 - U0) * G).astype(int), 0, G - 1), np.clip(((aa - A0) / (A1 - A0) * G).astype(int), 0, G - 1))
cu, ca = cell(u, a); Hf = np.zeros((G, G)); np.add.at(Hf, (cu, ca), dt / dt.sum()); su, sa = cell(u[sec], a[sec]); Hs = np.zeros((G, G)); np.add.at(Hs, (su, sa), 1)
void_s = Hs == 0; full_in_void = float(Hf[void_s].sum()); R['T3'] = dict(grid=G, cells=int(G * G), section_empty_cells=int(void_s.sum()), full_path_share_in_section_empty_cells=full_in_void, full_empty_cells=int((Hf == 0).sum()), section_occupied_cells=int((Hs > 0).sum()), full_occupied_cells=int((Hf > 0).sum()))
# T5: hyperbola bound
def Bp(t): N = int(math.sqrt(t / (2 * math.pi))); n = np.arange(1, N + 1); return 2 * np.sum(n ** -0.5 * (float(ZM.theta_prime(t)) + np.log(n)))
blk = [(0, 1000), (1000, 5000), (5000, 25000)]; rows = []
for lo, hi in blk:
    m = (np.arange(len(ts)) >= np.searchsorted(ts, ZM.zeros_mp(1)[0] if False else (ts[0] if lo == 0 else ts[min(len(ts) - 1, lo * 16)]))) & (np.arange(len(ts)) < min(len(ts), hi * 16)); K = np.abs(Zp[m]); tt = ts[m]; B = np.array([Bp(t) for t in tt[::2000]])
    rows.append(dict(zeros=[lo, hi], max_abs_Zprime=float(K.max()), bound_Bp_at_end=float(Bp(tt[-1])), max_over_bound=float(K.max() / Bp(tt[-1]))))
R['T5'] = rows
# Aizawa
o = AZ.orbit(600000, 3000, 0.0045); x, y, zz = o[:, 0], o[:, 1], o[:, 2]; r = np.hypot(x, y); ua = np.log(r); ra = zz - AZ.STANDARD['b']; aa = np.arcsinh(ra)
secA = np.array([i for i in range(1, len(y)) if y[i - 1] < 0 <= y[i] and x[i] > 0]); R['T4'] = dict(section_points=len(secA), KS_u_p=float(ks_2samp(ua[secA], ua[::50])[1]), KS_rate_p=float(ks_2samp(ra[secA], ra[::50])[1]))
cu_, ca_ = cell(ua, aa); HA = np.zeros((G, G)); np.add.at(HA, (cu_, ca_), 1); su_, sa_ = cell(ua[secA], aa[secA]); HS = np.zeros((G, G)); np.add.at(HS, (su_, sa_), 1)
# hole inside the ring: cells enclosed by the band: flood from the corner over empty cells, remaining empty cells are the hole
from scipy.ndimage import label
emp = HA == 0; lab, nlab = label(emp); outside = set(np.unique(lab[[0, -1], :].ravel())) | set(np.unique(lab[:, [0, -1]].ravel())); outside.discard(0)
hole = np.isin(lab, [k for k in range(1, nlab + 1) if k not in outside]); R['T4'].update(dict(full_occupied_cells=int((HA > 0).sum()), hole_cells=int(hole.sum()), section_points_in_hole=int(HS[hole].sum()), full_points_in_hole=int(HA[hole].sum()), u_range_full=[float(ua.min()), float(ua.max())], rate_range_full=[float(ra.min()), float(ra.max())]))
R['_meta'] = dict(note='phase portrait cells: u in [-6,3.3], asinh(rate) in [-7,7], 60 x 60')
json.dump(R, open(root / 'experiments/section_void_test.json', 'w'), indent=1, default=float); print(json.dumps(R, indent=0, default=float))
