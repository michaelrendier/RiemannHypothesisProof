# SPDX-License-Identifier: GPL-3.0-only
"""fermat_side_probe.py — first probe of the Fermat side with exact, established objects (Riemann-side tests' analogues).   :CALCULATED
 (1) L(s, χ₋₄) = the n = 2 Fermat character: zeros on ½ + it (Λ(s) = (π/4)^{-(s+1)/2} Γ((s+1)/2) L(s) is real there).
 (2) Dedekind zeta of ℚ(i) = ζ(s) L(s, χ₋₄) = Σ r₂(n)/4 n^{-s}: lines r₂(n)/4 at u = ln n; Λ_K(n) = Λ(n)(1 + χ₋₄(n)): primes p ≡ 1 mod 4 double, p ≡ 3 mod 4 vanish at p and appear at p².
 (3) Fermat curves xⁿ + yⁿ = zⁿ over F_p: Weil's bound |N − (p+1)| ≤ 2g√p, g = (n−1)(n−2)/2 (the Riemann hypothesis for the curve, PROVEN): ratio ≤ 1.
 (4) sums of two squares (the n = 2 'solution lattice') thin like 0.7642/√ln x (Landau–Ramanujan), primes like 1/ln x.
Writes experiments/fermat_side_probe.json."""
import json, math, os, sys
from collections import Counter
import numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rhp import nodal_extract as NX
mp.mp.dps = 20; chi = [0, 1, 0, -1]; out = {}
Lam = lambda t: (mp.pi / 4) ** (-(mp.mpc(0.5, t) + 1) / 2) * mp.gamma((mp.mpc(0.5, t) + 1) / 2) * mp.dirichlet(mp.mpc(0.5, t), chi)
out["lambda_imag_over_real_t7_t20"] = [float(abs(Lam(t).imag / Lam(t).real)) for t in (7, 20)]
ts = np.arange(0.5, 100, 0.05); f = np.array([float(Lam(t).real) for t in ts]); idx = np.nonzero(f[:-1] * f[1:] < 0)[0]
zs = [float(mp.findroot(lambda x: Lam(x).real, (ts[i], ts[i + 1]), solver="anderson", tol=1e-15)) for i in idx]
out["L_chi_m4_zeros_below_100"] = zs
r2q = lambda n: sum(chi[d % 4] for d in range(1, n + 1) if n % d == 0)
out["r2_over_4_n1_32"] = [r2q(n) for n in range(1, 33)]
psi, lam = NX.psi_table(60); LamK = lambda n: float(lam[n] * (1 + chi[n % 4]))
out["LambdaK"] = {str(n): LamK(n) for n in (2, 3, 4, 5, 7, 9, 13, 49)}
def count(n, p):
    c = Counter(pow(x, n, p) for x in range(p)); aff = sum(c[a] * c[(1 - a) % p] for a in c)
    inf = sum(c[a] * c[(-a) % p] for a in c if a != 0) // (p - 1)
    return aff + inf
wb = {}
for n in (3, 4, 5, 7):
    g = (n - 1) * (n - 2) // 2; ratios = []
    for p in range(5, 400):
        if all(p % d for d in range(2, int(p ** .5) + 1)) and p % n: ratios.append(abs(count(n, p) - (p + 1)) / (2 * g * math.sqrt(p)))
    wb[n] = dict(genus=g, max_ratio=max(ratios), primes_tested=len(ratios))
out["weil_bound"] = wb
X = 10 ** 6; s2 = np.zeros(X + 1, bool); a = np.arange(0, int(X ** .5) + 1); A = (a * a)[:, None] + (a * a)[None, :]; s2[A[A <= X]] = True
out["two_squares_density"] = {str(x): dict(density=float(s2[1:x + 1].sum() / x), landau_ramanujan=0.7642 / math.sqrt(math.log(x))) for x in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6)}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fermat_side_probe.json"), "w"), indent=1)
print(len(zs), "zeros;", {n: round(v["max_ratio"], 3) for n, v in wb.items()})
