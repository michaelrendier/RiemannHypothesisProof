"""Both objects seen from the pole (the view in Cody's riemann-path.png): the Aizawa orbit down its axis, and zeta(1/2+it) on the Riemann sphere down the Z axis."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from rhp import strange_path as P
pa, la = P.generate_geometry(P.aizawa(), points=120000, burn_in=3000, h=0.0045); A = pa.reshape(-1, 3)          # slots x, z, y
pz, lz, ts, w = P.zeta_sphere_path(points=9000, t0=0.0, t1=120.0); Z = pz.reshape(-1, 3)
fig, ax = plt.subplots(1, 3, figsize=(15, 5), facecolor='black')
for a in ax: a.set_facecolor('black'); a.set_aspect('equal'); a.axis('off')
ax[0].scatter(A[:, 0], A[:, 2], c=la, cmap='rainbow', s=0.05, alpha=0.5); ax[0].set_title('Aizawa orbit, looking down the z axis (the JS up axis)', color='white', fontsize=9)
ax[1].scatter(Z[:, 0], Z[:, 2], c=lz, cmap='rainbow', s=0.4, alpha=0.7); ax[1].set_title('zeta(1/2+it), t in [0,120], Riemann sphere seen from the pole', color='white', fontsize=9)
ax[2].scatter(A[:, 0], A[:, 1], c=la, cmap='rainbow', s=0.05, alpha=0.5); ax[2].set_title('Aizawa orbit, side view (x vs z)', color='white', fontsize=9)
fig.savefig('figures/poles_view.png', dpi=110, facecolor='black', bbox_inches='tight'); print('ok')
