"""Renders rhp.strange_path: the Aizawa RK4 orbit (the JavaScript's geometry) and the zeta critical-line path on the Riemann sphere.  Writes strange_path.json and two figures."""
import json, sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from rhp import strange_path as P
R = {}
def draw(pos, line, name, title, lim=None, elev=22, azim=-60):
    X = pos.reshape(-1, 3); fig = plt.figure(figsize=(7, 7), facecolor='black'); ax = fig.add_subplot(111, projection='3d', facecolor='black')
    ax.scatter(X[:, 0], X[:, 2] if False else X[:, 1], X[:, 2], c=line, cmap='rainbow', s=0.4, alpha=0.7, linewidths=0)       # slots are (x, z, y): plot (x, y=slot3, z=slot2) upright
    ax.set_axis_off(); ax.view_init(elev=elev, azim=azim); ax.set_title(title, color='white', fontsize=9)
    m = lim or np.abs(X).max(); ax.set_xlim(-m, m); ax.set_ylim(-m, m); ax.set_zlim(-m, m)
    fig.savefig(f'figures/{name}.png', dpi=130, facecolor='black', bbox_inches='tight'); plt.close(fig)
t0 = time.time(); pos, line = P.generate_geometry(P.aizawa(), points=30000, burn_in=2000, scale=1.0, h=0.01); R['aizawa'] = dict(points=30000, secs=round(time.time() - t0, 2), extent=float(np.abs(pos).max()))
X = pos.reshape(-1, 3); draw(np.stack([X[:, 0], X[:, 2], X[:, 1]], 1).ravel(), line, 'aizawa_rk4', 'Aizawa attractor, RK4 (the JavaScript geometry)', lim=1.6)
t0 = time.time(); pos, line, ts, w = P.zeta_sphere_path(points=9000, t0=0.0, t1=120.0)
absw = np.abs(w); mins = [i for i in range(1, len(w) - 1) if absw[i] < absw[i - 1] and absw[i] <= absw[i + 1] and absw[i] < 0.15]
from mpmath import zetazero
zeros = [float(zetazero(k).imag) for k in range(1, 30)]
R['zeta'] = dict(points=9000, t=[0.0, 120.0], secs=round(time.time() - t0, 1), pole_visits=len(mins), zeros_below_120=len([z for z in zeros if z < 120]), visit_ts=[round(float(ts[i]), 2) for i in mins[:6]], first_zeros=[round(z, 2) for z in zeros[:6]])
X = pos.reshape(-1, 3); draw(np.stack([X[:, 0], X[:, 2], X[:, 1]], 1).ravel(), line, 'zeta_sphere_path', 'zeta(1/2 + it), t in [0,120], on the Riemann sphere', lim=1.05)
json.dump(R, open('experiments/strange_path.json', 'w'), indent=1); print(json.dumps(R, indent=1))
