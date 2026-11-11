# SPDX-License-Identifier: GPL-3.0-only
"""nshape_shape.py — The Shape (rhp/nshape.py): the Ω-onion, the two factoring conics, the tower of division-algebra lattices.   :CALCULATED
Writes experiments/nshape_shape.json.   Run: python3 experiments/nshape_shape.py"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
from rhp import nshape as S
out = {}; t0 = time.time()

# ---- A. the Ω-onion to 10^7
X = 10 ** 7; om = S.big_omega(X); xs = [10 ** k for k in range(3, 8)]; lc = S.layer_counts(om, xs, kmax=9)
rows = {}
for x in xs:
    c = lc[x]; tot = x - 1
    rows[str(x)] = dict(counts=c, fractions=[v / tot for v in c], mode=int(np.argmax(c)) + 1, prime_fraction=c[0] / tot, one_over_ln=1 / math.log(x),
                        landau_ratio=[c[k - 1] / S.landau(x, k) for k in (1, 2, 3, 4)], mean_omega=float(om[2:x + 1].mean()), lnln=math.log(math.log(x)))
out["onion"] = rows; print("A done", round(time.time() - t0)); sys.stdout.flush()

# ---- B. the two conics
NB = 10 ** 6; P = S.primes_upto(NB); isp = np.zeros(NB + 1, bool); isp[P] = True
n = np.arange(NB + 1); odd = (n % 2 == 1); oddcomp = odd & ~isp & (n > 1)
h = S.hyperbola_points(NB)
out["hyperbola"] = dict(N_max=NB, lattice_reach_equals_odd_composites=bool((h[odd & (n > 1)] == oddcomp[odd & (n > 1)]).all()), odd_composites=int(oddcomp.sum()))
NS = 30000; hs = S.hyperbola_search(NS)
out["hyperbola"]["direct_search_equals_odd_composites_to"] = NS; out["hyperbola"]["direct_search_ok"] = bool((hs[odd[:NS + 1] & (n[:NS + 1] > 1)] == oddcomp[:NS + 1][odd[:NS + 1] & (n[:NS + 1] > 1)]).all())
cnt = S.circle_primitive_reps(NB)
# theory: odd N, primitive reps x>y>=1 coprime = 2^(w1-1) if every prime factor ≡ 1 mod 4 else 0
w1 = np.zeros(NB + 1, np.int8); bad = np.zeros(NB + 1, bool)
for p in P[1:]:
    p = int(p); w1[p::p] += (1 if p % 4 == 1 else 0)
    if p % 4 == 3: bad[p::p] = True
pred = np.where(odd & ~bad & (n > 1), 2 ** np.maximum(w1.astype(np.int64) - 1, 0), 0); pred[w1 == 0] = 0
mask = odd & (n > 1)
out["circle"] = dict(theory_matches_all_odd_N=bool((cnt[mask] == pred[mask]).all()), mismatches=int((cnt[mask] != pred[mask]).sum()))
two = mask & (cnt >= 2); one = mask & (cnt == 1)
out["circle"].update(odd_composites_with_2plus_reps=int((two & oddcomp).sum()), odd_primes_with_1_rep=int((one & isp).sum()), composites_with_exactly_1_rep=int((one & oddcomp).sum()),
                     fraction_of_odd_composites_reached_by_euler=float((two & oddcomp).sum() / oddcomp.sum()), fraction_of_hyperbola=1.0 if out["hyperbola"]["lattice_reach_equals_odd_composites"] else None,
                     all_2plus_are_composite=bool((two & isp).sum() == 0))
Ns = np.nonzero(two & (n <= 200000))[0]; fails = 0
for N in Ns:
    g = S.euler_factor(int(N)); fails += (g is None or not (1 < g < N) or N % g != 0)
out["circle"]["euler_factor_checked"] = int(len(Ns)); out["circle"]["euler_factor_failures"] = int(fails)
print("B done", round(time.time() - t0)); sys.stdout.flush()

# ---- C. the tower
ps = [int(p) for p in P[P < 3000]]
red = {"Z": 0, "Zi": sum(1 for p in ps if p == 2 or p % 4 == 1), "Hurwitz": 0, "octonion": 0}
lag = all(any(math.isqrt(p - a * a - b * b - c * c) ** 2 == p - a * a - b * b - c * c for a in range(math.isqrt(p) + 1) for b in range(math.isqrt(p - a * a) + 1) for c in range(math.isqrt(p - a * a - b * b) + 1)) for p in ps[:60])
zi_check = all(S.sum_of_squares_exists(p, 2) == (p == 2 or p % 4 == 1) for p in ps)
red["Hurwitz"] = len(ps) if lag else None; red["octonion"] = len(ps) if lag else None
big = P[P < NB]; hole_zi = float(np.mean(~((big == 2) | (big % 4 == 1))))
out["tower"] = dict(primes_tested=len(ps), reducible=red, lagrange_search_ok_first60=bool(lag), gaussian_rule_ok=bool(zi_check),
                    units={"Z": 2, "Zi": 4, "Hurwitz": 24, "E8": 240},
                    hole_fraction_among_primes_below_1e6={"Z": 1.0, "Zi": hole_zi, "Hurwitz": 0.0, "octonion": 0.0},
                    hurwitz_norm_counts={str(m): S.count_hurwitz_norm(m) for m in (1, 2, 3, 5, 7, 11, 13)},
                    hurwitz_expected_24_p_plus_1={str(p): 24 * (p + 1) if p > 2 else 24 for p in (2, 3, 5, 7, 11, 13)},
                    e8_shells={str(m): S.e8_shell(m) for m in (2, 4, 6, 8)}, e8_expected_240_sigma3={str(m): 240 * sum(d ** 3 for d in range(1, m // 2 + 1) if (m // 2) % d == 0) for m in (2, 4, 6, 8)})
out["tower"]["irreducibles_per_rational_prime_up_to_units"] = {"Z": "1", "Zi": "2 if p ≡ 1 (4); 1 if p ≡ 3 (4) (the element p, norm p²); 1 if p = 2", "Hurwitz": "p + 1 (p odd); 1 (p = 2)", "E8_octonions": "1 + p³ (norm-p elements / 240; the product is not associative)"}
json.dump(out, open(os.path.join(HERE, "nshape_shape.json"), "w"), indent=1)
print("C done", round(time.time() - t0))
for x in xs: r = out["onion"][str(x)]; print(x, "layers 1..9+ share:", [round(f, 3) for f in r["fractions"]], "mode", r["mode"], "mean Ω", round(r["mean_omega"], 2), "Landau ratios", [round(v, 2) for v in r["landau_ratio"]])
print(out["hyperbola"]); print(out["circle"]); print({k: v for k, v in out["tower"].items() if k not in ("irreducibles_per_rational_prime_up_to_units",)})
