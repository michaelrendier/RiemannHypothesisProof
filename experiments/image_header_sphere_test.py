"""Run the numbers of the repo-header image (Gemini 'The Riemann Sphere') against the maths.
SPDX-License-Identifier: GPL-3.0-only.  Run: ../ValaQuenta/.venv/bin/python experiments/image_header_sphere_test.py
Tests: (A) the image's pole/axis labels on the Cayley chart; (B) hemisphere shares of the VALUE sphere
(zeta(1/2+it) through stereographic projection: equator |z|=1, S pole = zero, N pole = infinity) and the
DC term (n=1); (C) the twist/snap count: theta(gamma_n)/pi against n (a snap per pi of twist), scatter = primes."""
import json, sys
import numpy as np
from mpmath import mp, mpf, mpc, zeta, siegeltheta, zetazero
mp.dps = 20
out = {}
# (A) labels
w = lambda s: (s - 1) / s
A = {"s=1 (pole of zeta)": "w=0 -> north pole", "s=0": "w=inf -> south pole", "s=inf": "w=1 -> ON THE EQUATOR (not a pole)"}
triv = []
for k in range(1, 6):
    ww = complex(w(-2 * k)); lat = np.degrees(np.arcsin((abs(ww)**2 - 1) / (abs(ww)**2 + 1)))  # |w|>1 -> northern in this formula; sign flipped below
    triv.append({"s": -2 * k, "w": ww, "|w|": abs(ww), "lat_deg_south": -abs(lat) if abs(ww) > 1 else abs(lat)})
nt = [abs(complex(w(mpc(0.5, float(zetazero(n).imag))))) for n in range(1, 11)]
A["trivial zeros s=-2k"] = triv
A["|w| at first 10 nontrivial zeros (should be 1)"] = nt
out["A_labels"] = A
# (B) value sphere
T = 400.0; t = np.arange(0.5, T, 0.02)
z = np.array([complex(zeta(mpc(0.5, x))) for x in t])
def share(c): return float(np.mean(c))
mod = np.abs(z)
B = {"samples": len(t), "mean zeta": complex(z.mean()), "share Re>0": share(z.real > 0), "share |zeta|>1 (north of value-equator)": share(mod > 1),
     "mean ln|zeta| (Selberg: 0)": float(np.mean(np.log(mod))),
     "share Re(zeta-1)>0 (n=1 DC removed)": share((z - 1).real > 0), "mean(zeta-1)": complex((z - 1).mean())}
# Riemann-Siegel pieces
th = np.array([float(siegeltheta(mpf(x))) for x in t])
N = np.floor(np.sqrt(t / (2 * np.pi))).astype(int)
S = np.array([sum(n ** (-0.5 - 1j * x) for n in range(1, Nn + 1)) for x, Nn in zip(t, N)])
mirror = np.exp(-2j * th) * np.conj(S)
B["mean S (n<=N)"] = complex(S.mean()); B["mean mirror term e^{-2i theta} conj(S)"] = complex(mirror.mean())
B["share Re S>0"] = share(S.real > 0); B["share Re mirror>0"] = share(mirror.real > 0)
Z = np.real(np.exp(1j * th) * z)
B["mean Z (AC?)"] = float(Z.mean()); B["mean |Z| (rectified)"] = float(np.abs(Z).mean()); B["Z>0 share"] = share(Z > 0)
out["B_value_sphere"] = B
# (C) twist per snap
zs = np.array([float(x.split("\t")[1]) for x in open("data/zeros_25000_hp.tsv").read().splitlines()[1:]])
ths = np.array([float(siegeltheta(mpf(g))) for g in zs])
x = ths / np.pi - np.arange(0, len(zs))      # twist in units of pi, minus index
C = {"mean theta(gamma_n)/pi - (n-1)": float(x.mean()), "std": float(x.std()), "first100 std": float(x[:100].std()), "last5000 std": float(x[-5000:].std()),
     "Selberg std of S(T) at T=gamma_25000": float(np.sqrt(np.log(np.log(zs[-1])) / 2) / np.pi),
     "mean spacing in theta (units of pi)": float(np.diff(ths).mean() / np.pi), "std of spacing (units of pi)": float(np.diff(ths).std() / np.pi)}
out["C_twist_per_snap"] = C
json.dump(out, open("experiments/image_header_sphere_test.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
