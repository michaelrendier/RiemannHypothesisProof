# The Shape — the Extinguishable Bulk around the holes   2026-11-10   :CALCULATED / ESTABLISHED as marked (no claim about RH)

Code `rhp/nshape.py`, run `experiments/nshape_shape.py` → `experiments/nshape_shape.json`, figure `rhp/atlas_riemann_fermat/nshape.svg`, test `test_nshape_shape`. Scope (Cody, 2026-11-10): Fermat's side of the boundary as a description of a continuous Laurelin bulk (everything a factor reaches) with the primes as the holes complete factorisation cannot reach (N-Shape / N-Holes). No J operator, no 0_RB, no α; recursive sieving is not built.

## 1 The Ω-onion (ESTABLISHED asymptotics, counts :CALCULATED, n ≤ 10⁷)
Layer k = {n ≥ 2 : Ω(n) = k}. Share of n ≤ x:
| x | Ω=1 (primes) | 2 | 3 | 4 | 5 | 6 | mode | mean Ω |
|---|---|---|---|---|---|---|---|---|
| 10³ | 0.168 | 0.299 | 0.247 | 0.149 | 0.076 | 0.037 | 2 | 2.88 |
| 10⁵ | 0.096 | 0.234 | 0.256 | 0.187 | 0.112 | 0.059 | 3 | 3.44 |
| 10⁷ | 0.066 | 0.190 | 0.244 | 0.205 | 0.135 | 0.077 | 3 | 3.79 (ln ln x + 1.03 = 3.81) |
The core (the holes) thins like 1/ln x (0.066 against 1/ln 10⁷ = 0.062); the bulk's weight moves outward, mode layer 2 → 3 at ≈ 10⁵. Landau's leading term x(ln ln x)^(k−1)/((k−1)! ln x) converges slowly: measured/leading = 1.07 / 1.10 / 1.02 / 0.92 for k = 1…4 at 10⁷.

## 2 The two factoring conics (exhaustive, N ≤ 10⁶)
- **Fermat's hyperbola** a² − b² = N, lattice points with a − b ≥ 3: enumerated from the lattice, they are exactly the 421,502 odd composites ≤ 10⁶ and no prime; a direct search over a (no factoring) agrees for every odd N ≤ 30,000. A prime has only the trivial point ((N+1)/2, (N−1)/2): in this shape the primes are the holes, and every composite is reached.
- **Euler's circle** x² + y² = N: the number of primitive representations (x > y ≥ 1) equals 2^(ω₁−1) when every prime factor of the odd N is ≡ 1 mod 4 and 0 otherwise, for all odd N ≤ 10⁶ (0 mismatches). Two or more representations occur only for composites (48,606 of the 421,502 odd composites, 11.5 %; share flat at 0.10–0.12 from 10³ to 10⁶); Euler's gcd(N, x₁y₂ − x₂y₁) gave a non-trivial factor for all 9,660 such N ≤ 2·10⁵ (0 failures). Primes ≡ 1 mod 4 (39,175) have exactly one representation, and so do 100 composites (prime powers p^k, k ≥ 2, p ≡ 1 mod 4): the circle cannot tell them from primes.
- The circle reaches a small, structured part of the bulk (composites built only from primes ≡ 1 mod 4 with two or more distinct ones); the hyperbola reaches all of it.

## 3 The tower of division-algebra lattices (ESTABLISHED theta series, checked by enumeration)
| level | lattice | units | norm-n element count (checked) | rational primes that stay holes | irreducible classes per prime p |
|---|---|---|---|---|---|
| ℤ | Z | 2 | – | all (100 %) | 1 |
| ℤ[i] | Z² | 4 | 4 Σ χ₋₄(d); rule p = a² + b² verified for all 430 primes < 3000 | those p ≡ 3 mod 4: 50.09 % of the primes < 10⁶ | 2 (p ≡ 1 mod 4), 1 (p ≡ 3 mod 4: the element p, norm p²), 1 (p = 2) |
| Hurwitz ℍ | D4 | 24 | 24(p+1) for odd p, verified p = 3…13 (96, 144, 192, 288, 336) and 24 for p = 2 by counting elements | 0 (four squares for every p; Lagrange search checked for the first 60 primes) | p + 1 (p odd), 1 (p = 2) |
| integral octonions | E8 | 240 | 240σ₃(n): shells 240 / 2160 / 6720 / 17520 for squared lengths 2, 4, 6, 8 verified by brute force | 0 (eight squares) | 1 + p³ (norm-p elements / 240; not associative, so 'class' is a count, not a factorisation) |
Reading: at each level up, rational primes are absorbed into the bulk (1 → ½ → 0 → 0) while the level's own irreducibles (elements of prime norm in the quaternions and octonions) multiply (1, ≤ 2, p + 1, p³ + 1 per prime): the shape grows holes of its own kind as the rational ones close. At dimension 16 the zero divisors end factorisation (tables to 128 dimensions in `rhp/sedenion_slots.py`; zero-divisor pair counts 42 / 294 / 1,518 / 6,942). The Dedekind zeta ζ_ℚ(i) = ζ·L(χ₋₄) is the level-2 statement: χ₋₄(p) = +1 is a prime that factors (bulk), −1 a hole, 0 for 2.

## 4 What the three views share
In each, **the bulk is whatever a factor reaches, and it is certified by exhibiting one factor** (one multiplication); the holes (primes in ℤ; inert primes in ℤ[i]; prime-norm elements in ℍ) are the elements for which no such witness exists and must be certified by exhausting the factors or by an algebraic witness. The bulk's reach is complete for the hyperbola, partial for the circle, and total for the quaternions; the holes' number depends on the level (1, 1/2, 0 of the rational primes; the level's own irreducibles grow).
Open (stated, not tested): whether the hole at level k is best taken as the rational primes (absorbed one level at a time) or the level's irreducibles (a changing set); both are measured above.
