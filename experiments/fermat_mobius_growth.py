"""The inverse (Möbius) side of the Fermat finite-tone wave: coefficients c_k of 1/P(T), P(T) = ∏(1 − ω_j T), normalised by p^{k/2}.  Simple roots → bounded; a root of order m → growth ~ k^{m−1}.  Writes fermat_mobius_growth.json."""
import json, math, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import flow as FL
R = {}
for n, p in ((3, 7), (4, 17), (5, 31), (7, 29)):
    th = FL.frobenius_angles(p, n); poly = np.poly1d([1.0 + 0j])
    for t in th: poly = poly * np.poly1d([-np.exp(1j * t), 1.0])          # (1 − e^{iθ} T), unit-modulus roots (T → T/√p)
    a = poly.coeffs[::-1].real; deg = len(a) - 1; c = [1.0]
    for k in range(1, 161): c.append(-sum(a[i] * c[k - i] for i in range(1, min(k, deg) + 1)))
    c = np.array(c); key = np.round(np.mod(th, 2 * math.pi), 7); _, cnt = np.unique(key, return_counts=True); m = int(cnt.max())
    env = [float(np.max(np.abs(c[lo:hi]))) for lo, hi in ((10, 20), (20, 40), (40, 80), (80, 160))]
    slope = [math.log(env[i + 1] / env[i]) / math.log(2) for i in range(3)]
    R[f"n={n}, p={p}"] = dict(max_multiplicity=m, envelope_max_abs_c=env, growth_exponent_between_doublings=slope, expected_exponent=m - 1)
    print(f"n={n} p={p}: max root multiplicity {m}; max|c_k|/p^(k/2) over k in [10,20),[20,40),[40,80),[80,160): {[round(e,2) for e in env]}; local exponents {[round(s,2) for s in slope]} (expected {m-1})")
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
