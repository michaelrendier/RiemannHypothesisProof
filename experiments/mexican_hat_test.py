"""Does the Aizawa / zeta portrait contain a Mexican-hat potential?  Writes mexican_hat_test.json.
Expectations written before the run:
  H1  dz/dt = -dV/dz with V(r, z) = z^4/12 - a z^2/2 - c z + r^2 z + (e/2) r^2 z^2 exactly (so z moves in a tilted quartic double well, a hat cross-section), but dr/dt = (z - b) r is NOT -dV/dr:
      the system is not a gradient flow.
  H2  the tilted double well has minima at the roots of c + a z - z^3/3 where z^2 > a (z = 1.943 and -1.105) and a hump at -0.838; the orbit's z range sits inside the right-hand well and the hump region.
  H3  Aizawa's position density in the (x, y) plane is a RING with a HOLE (sombrero-like topology: V_eff = -ln rho has a ring valley and an infinite central cap); zeta's value-plane density has no hole: it rises as 1/|zeta| toward the
      origin (V_eff = ln|zeta|, the 2-D point-mass potential in form), which is kinematic: any smooth curve crossing a point at nonzero speed does it, so it is not evidence of gravity."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.aizawa_attractor import maths as AZ
from ValaQuenta.modules.zeta_uns_path import maths as ZM
P = AZ.STANDARD; a, b, c, d, e, f = (P[k] for k in 'abcdef'); R = {}
V = lambda r, z: z ** 4 / 12 - a * z ** 2 / 2 - c * z + r ** 2 * z + e / 2 * r ** 2 * z ** 2
rng = np.random.default_rng(3); res_z = []; res_r = []
for _ in range(200):
    r, z = rng.uniform(0.1, 1.5), rng.uniform(-1, 2); h = 1e-6; dVdz = (V(r, z + h) - V(r, z - h)) / (2 * h); dVdr = (V(r + h, z) - V(r - h, z)) / (2 * h)
    zdot = c + a * z - z ** 3 / 3 - r * r * (1 + e * z); rdot = (z - b) * r; res_z.append(abs(zdot + dVdz)); res_r.append(abs(rdot + dVdr))
R['H1'] = dict(max_abs_residual_zdot_eq_minus_dV_dz=float(max(res_z)), median_abs_residual_rdot_eq_minus_dV_dr=float(np.median(res_r)), note='f = 0 part; the f z x^3 term adds a phi-dependent force')
Vz = lambda z: z ** 4 / 12 - a * z ** 2 / 2 - c * z; fps = AZ.fixed_points(); wells = []
for fp in fps: z0 = float(fp['point'][2]); wells.append(dict(z=z0, V=float(Vz(z0)), V_second_derivative=float(z0 ** 2 - a), kind='minimum' if z0 ** 2 - a > 0 else 'maximum (hump)'))
o = AZ.orbit(300000, 3000, 0.0045); x, y, zz = o[:, 0], o[:, 1], o[:, 2]; r = np.hypot(x, y)
R['H2'] = dict(wells=wells, orbit_z_range=[float(zz.min()), float(zz.max())], V_at_orbit_z_min=float(Vz(zz.min())), V_at_orbit_z_max=float(Vz(zz.max())), vacuum_like_z0_for_c_zero=float(math.sqrt(3 * a)))
# H3: radial density in the plane
edges = np.linspace(0, 1.6, 65); cen = 0.5 * (edges[1:] + edges[:-1]); H, _ = np.histogram(r, bins=edges); dens = H / (2 * math.pi * cen * np.diff(edges)); dens /= dens.max(); ipk = int(np.argmax(dens))
R['H3_aizawa'] = dict(hole_radius=float(r.min()), ring_peak_radius=float(cen[ipk]), density_at_peak=1.0, density_at_0p5_of_peak_inner_edge=float(cen[np.argmax(dens > 0.5)]), V_eff_min_at_r=float(cen[ipk]), r_max=float(r.max()))
ts = ZM.sample_by_phase(10.0, 25000, 16, jitter=True); z = ZM.zeta_line(ts); ar = np.abs(z); dt = np.gradient(ts); w = dt / dt.sum()
e2 = np.logspace(-2.3, 0.5, 30); Hz, _ = np.histogram(ar, bins=e2, weights=w); cz = np.sqrt(e2[1:] * e2[:-1]); rho = Hz / (math.pi * (e2[1:] ** 2 - e2[:-1] ** 2)); m = cz < 0.3; sl = np.polyfit(np.log(cz[m]), np.log(rho[m]), 1)
R['H3_zeta'] = dict(slope_log_density_vs_log_r_for_r_lt_0p3=float(sl[0]), expected_for_a_curve_through_the_origin=-1.0, density_at_r_0p02=float(rho[np.argmin(np.abs(cz - 0.02))]), density_at_r_0p3=float(rho[np.argmin(np.abs(cz - 0.3))]), hole=False,
                    note='rho ~ 1/r follows from crossing 0 at nonzero speed (time within eps is proportional to eps); V_eff = -ln rho = ln r + const has the form of the 2-D point-mass potential only formally')
# perturbation: ring broadening by f
p0 = dict(P, f=0.0); s = np.array([0.1, 0.0, 0.0])
for _ in range(60000): s = AZ.rk4_step(s, 0.0045, p0)
q = np.empty((60000, 3))
for i in range(60000): s = AZ.rk4_step(s, 0.0045, p0); q[i] = s
rq = np.hypot(q[:, 0], q[:, 1]); R['perturbation_f'] = dict(f0_r_range=[float(rq.min()), float(rq.max())], f0_ring_thickness_ln_r=float(np.log(rq.max() / rq.min())), f01_r_range=[float(r.min()), float(r.max())], f01_ring_thickness_ln_r=float(np.log(r.max() / r.min())))
json.dump(R, open(root / 'experiments/mexican_hat_test.json', 'w'), indent=1, default=float); print(json.dumps(R, indent=0, default=float))
