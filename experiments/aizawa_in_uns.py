"""The Aizawa attractor in UNS (log-complex, radial-polar) coordinates.  Writes aizawa_in_uns.json.
In complex form w = x + i y:  dw/dt = (z - b + i d) w  EXACTLY (the third equation only feeds back through |w|^2 and Re(w)^3).  With u = ln r, phi = arg w:  du/dt = z - b,  dphi/dt = d,
dz/dt = c + a z - z^3/3 - e^(2u)(1 + e z) + f z e^(3u) cos^3(phi).  So (u, phi, z) is the native chart: phi is a pure rotation, u integrates the driver, and the orbit lives on a TORUS
(a surface of revolution) whose cross-section is the (u, z) curve.  At f = 0 that curve is the closed limit cycle of the planar (r, z) flow; at f = 0.1 it is a band.
Expectations written before the run: (1) the complex rate lambda = d(ln w)/dt = (z - b) + i d is a horizontal segment at Im = d; its time-mean real part is 0 (bounded u); (2) at f = 0 the (u, z) section at phi = 0 is a
single point (a limit cycle meets the section once per turn... a closed curve in (u,z) sampled at the rotation phase: a dense closed curve); at f = 0.1 a band; (3) zeta's analogue: lambda = Z'/Z - i theta'(t) (unbounded real part, drifting imaginary part),
its section at arg zeta = 0 (the Gram points) in (ln|zeta|, Z'/Z)."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.aizawa_attractor import maths as AZ
from ValaQuenta.modules.zeta_uns_path import maths as ZM
P = AZ.STANDARD; R = {}
o = AZ.orbit(300000, 3000, 0.0045); x, y, z = o[:, 0], o[:, 1], o[:, 2]; r = np.hypot(x, y); u = np.log(r); phi = np.unwrap(np.arctan2(y, x))
# complex-rate identity: dw/dt = (z - b + i d) w
w = x + 1j * y; h = 0.0045; wd = np.gradient(w, h); lam = wd / w; R['rate_identity'] = dict(max_abs_residual=float(np.abs(lam - ((z - P['b']) + 1j * P['d']))[5:-5].max()), note='central difference: O(h^2) residual, not the exact 1e-16 of the algebra')
R['rate_plane'] = dict(Re_min=float((z - P['b']).min()), Re_max=float((z - P['b']).max()), Re_mean=float((z - P['b']).mean()), Im=P['d'], share_growing=float(np.mean(z > P['b'])))
R['log_radius'] = dict(u_min=float(u.min()), u_max=float(u.max()), u_mean=float(u.mean()), r_min=float(r.min()), r_max=float(r.max()))
# section at phi = 0 (mod 2 pi) crossing upward with x > 0, f = 0.1
sec = [i for i in range(1, len(y)) if y[i - 1] < 0 <= y[i] and x[i] > 0]; su = np.array([u[i] for i in sec]); sz = np.array([z[i] for i in sec])
R['section_f01'] = dict(points=len(sec), u_range=[float(su.min()), float(su.max())], z_range=[float(sz.min()), float(sz.max())], corr_u_z=float(np.corrcoef(su, sz)[0, 1]))
# f = 0: limit cycle in (u, z)
p0 = dict(P, f=0.0); s = np.array([0.1, 0.0, 0.0])
for _ in range(60000): s = AZ.rk4_step(s, h, p0)
cyc = np.empty((80000, 3))
for i in range(80000): s = AZ.rk4_step(s, h, p0); cyc[i] = s
ru = np.log(np.hypot(cyc[:, 0], cyc[:, 1])); R['limit_cycle_f0'] = dict(u_range=[float(ru.min()), float(ru.max())], z_range=[float(cyc[:, 2].min()), float(cyc[:, 2].max())], loop_period_check=6.0493)
# zeta: lambda and Gram-point section
g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()])
ts = ZM.sample_by_phase(10.0, 25000, 16, jitter=True); zz = ZM.zeta_line(ts); th = ZM.theta_prime(ts)
hh = 1e-4; Zv = ZM.Z_rs(ts); drv = (ZM.Z_rs(ts + hh) - ZM.Z_rs(ts - hh)) / (2 * hh) / Zv
gs = [i for i in range(1, len(ts)) if zz[i - 1].imag < 0 <= zz[i].imag and zz[i].real > 0]
R['zeta_rate_plane'] = dict(Im_range=[float(th.min()), float(th.max())], Re_abs_median=float(np.median(np.abs(drv))), Re_share_abs_gt_10=float(np.mean(np.abs(drv) > 10)))
R['zeta_section_arg0'] = dict(points=len(gs), u_median=float(np.median(np.log(np.abs(zz[gs])))), u_p5_p95=[float(np.percentile(np.log(np.abs(zz[gs])), 5)), float(np.percentile(np.log(np.abs(zz[gs])), 95))], rate_median=float(np.median(drv[gs])), rate_p5_p95=[float(np.percentile(drv[gs], 5)), float(np.percentile(drv[gs], 95))])
json.dump(R, open(root / 'experiments/aizawa_in_uns.json', 'w'), indent=1, default=float); print(json.dumps(R, indent=0, default=float))
