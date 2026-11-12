# The Fermat side of the horizon, made exact — 2026-11-12   :CALCULATED

Code `rhp/fermat_curve.py`; run `experiments/fermat_horizon.py` → `.json` (≈ 20 min, n = 2 … 20); test `test_fermat_curve`.

## What the record holds on the inertia half (read 2026-11-12)
- **Text, no code.** `Ainulindale/AgeSecond/Second_Age_Ainulindale_Conjecture.md` Part III: start at the Speed of Causality ceiling, work backwards through the Berry–Keating domain, chasing inertia; the domain floor is A_π = 1/137.035999…, 'the minimum inertia threshold'. Wiki 17 names a 'Fermat prime generator'; no repository code implements it (searched: only `ValaQuenta/modules/constants/maths.py::derive_alpha_fermat`, which holds the causality bound v₁ = αc < c and the cardioid golden angle 360°/φ² = 137.507°, 0.34 % from 1/α, and hard-codes 1/137.035999084).
- Older papers attribute A_π to 'E8/Wyler n-ball geometry'. Wyler (1971): α = (9/8π⁴)(π⁵/(2⁴·5!))^{1/4} = 1/137.0361 (≈ 6·10⁻⁷ from the measured value; numerology by most accounts).
- Cody, 2026-11-12: the original Lambert-W experiment found ONE number (Ω); the four d* and three W values came from later testing; the inertia half is the same chase run on the FLT side, and A_fermat was found there. The λ-family table in docs/omega_dstar is a later-testing layer.
- Scale discrepancy: the 'thermal ceiling' is 1.4·10¹⁷ K (Second Age, 'about 140 quadrillion') in the text, 140e15 °F in `conjecture_proof.py`, and T_Planck = 1.4168·10³² K in the Ω paper.
So the inertia chase cannot be reproduced as run: it was never written as an experiment. Nothing here derives α.

## Exact Fermat-side objects (ESTABLISHED theory, checked numerically)
For p ≡ 1 mod n the affine count N = #{xⁿ + yⁿ = 1 over F_p} = p + Σ_{a,b=1}^{n−1} J(χ^a, χ^b). 139 cases (n = 3 … 7, p ≤ 400): the sum reproduces the brute-force count exactly, and the (n−1)(n−2) = 2g generic terms have |J| = √p to 1e−6 (the other n − 1 terms have modulus 1). So the Fermat curve's prime wave is a FINITE sum of 2g tones of fixed modulus √p; RH for it is a theorem (Weil). The ceiling is 2g√p; the main term p overtakes it at p > (g + √(g² + 3n − 1))² (the Weil threshold, ≈ (n−1)²(n−2)²), beyond which a non-trivial point (xyz ≠ 0) is guaranteed.

## The primes the Fermat equation cannot reach (certified complete up to the Weil threshold)
| n | 2g | Weil threshold | exceptional primes (no solution with xyz ≠ 0 mod p) | largest ÷ threshold |
|---|---|---|---|---|
| 2 | 0 | 5 | 3, 5 | 1.00 |
| 3 | 2 | 16 | 2, 7, 13 | 0.81 |
| 4 | 6 | 56 | 3, 5, 13, 17, 41 | 0.73 |
| 5 | 12 | 171 | 2, 11, 41, 71, 101 | 0.59 |
| 6 | 20 | 433 | 5, 7, 13, 19, 43, 61, 97, 157, 277 | 0.64 |
| 7 | 30 | 940 | 2, 29, 71, 113, 491 | 0.52 |
| 8 | 42 | 1,810 | 3, 5, 13, 17, 41, 113 | 0.06 |
| 9 | 56 | 3,188 | 2, 7, 13, 19, 37, 73, 181, 523, 577 | 0.18 |
| 10 | 72 | 5,242 | 3, 11, 31, 41, 71, 101, 281, 401, 1181 | 0.23 |
| 11 | 90 | 8,164 | 2, 23, 89, 947, 1409 | 0.17 |
| 13 | 132 | 17,500 | 2, 53, 131, 443, 521 | 0.03 |
| 17 | 240 | 57,700 | 2, 137, 239, 443, 1667, 4013 | 0.07 |
| 19 | 306 | 93,748 | 2, 191, 419, 647, 761, 1217, 1901, 2129 | 0.02 |
(n = 12, 14–16, 18, 20 are in the JSON: 21, 17, 25, 12, 22, 24 exceptional primes, largest 1,201 / 2,801 / 2,851 / 1,409 / 3,061 / 4,241.)
- The set is finite, small (2 to 25 primes), and its largest member sits well inside the Weil ceiling: the ratio falls from 1.00 (n = 2) to 0.02–0.1 for n ≥ 13. Largest prime ≈ 1.4 n² (n = 3) … ≈ 11–14 n² (n = 10, 14, 17, 20), not a clean law and below the n⁴ Weil bound.
- For prime n the exceptional primes are 2 and primes ≡ 1 mod n (otherwise every nonzero residue is an n-th power and 1 + 1 = 2 works); for composite n only gcd(n, p−1) matters.
- Even the Pythagorean case n = 2 has holes: 3 and 5 are too small to host a triple.

## Reading (not a result)
This is the Fermat-side counterpart of 'the places factorization cannot reach': a finite set of primes unreached by the bulk of solutions, bounded by a ceiling. The N-shape's extension in n adds tones (2g = (n−1)(n−2)) as it adds structure. Whether any threshold here relates to A_π or to the Ω-side crossing is open; no numerical relation was found or looked for beyond the ratios above.
