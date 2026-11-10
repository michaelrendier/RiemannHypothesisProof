# Sedenion slot rule — which prime sits on which slot of T_N?   2026-11-10   :CALCULATED (no claim about RH)

Question carried since 2026-11-08: wiki 130 maps slot k ↔ the k-th prime for k ≤ 15 only. What puts the next primes on slots, and does anything in the algebra select a rule? "numerology until a rule exists." Code: `rhp/sedenion_slots.py`; run: `experiments/sedenion_slot_rule.py` → `experiments/sedenion_slot_rule.json`; tests `test_sedenion_*`.

## 1 Machinery (exact; N = 16, 32, 64, 128)
- Cayley–Dickson sign table e_i e_j = S[i,j] e_{i⊕j} by doubling, convention of `emerger_spectrum.cd_mul`. At N = 16 it matches `cd_mul` on all 256 basis products (0 mismatches). Norm is multiplicative at N = 8 (ratio 1.0000000000000) and not at 16 (1.0216 on a random pair), as it must be.
- L_k = left multiplication = a signed permutation, L_k² = −I and L_k skew for every k ≥ 1 at all four N: every slot is a circle.
- Pairs of slots whose circles do not anticommute ("deforming", the Assessor-type pairs): **42 of 105 (N=16), 294 of 465 (32), 1518 of 1953 (64), 6942 of 8001 (128)** = 40 %, 63 %, 78 %, 87 %. Anticommuting: 63, 171, 435, 1059. (42 and 63 reproduce wiki 130 R2.)

## 2 Capacity: where the slots run out
| T_N | slots | seq0 (slot 0 ↔ 2, wiki 130 table): last slotted prime / first unslotted | seq1 (k-th prime on slot k, k ≥ 1): last / first unslotted |
|---|---|---|---|
| T_16 | 16 | 53 / 59 | 47 / 53 |
| T_32 | 32 | 131 / 137 | 127 / 131 |
| T_64 | 64 | **311 / 313** | 307 / 311 |
| T_128 | 128 | 719 / 727 | 709 / 719 |
311 is the 64th prime, 313 the 65th, 719 the 128th, 727 the 129th. Note seq0 puts prime 2 on slot 0 = the identity e₀, which has no circle (L₀ = I): under that table prime 2 cannot rotate. Under seq1 slot 0 is empty and 15 / 31 / 63 / 127 primes are slotted.

## 3 Stragglers (the composites left over after 313)
Composites with no prime factor among the slotted primes: smallest = (first unslotted prime)², i.e. **313² = 97,969 for T_64** (seq0), 727² = 528,529 for T_128, 59² = 3,481 for T_16, 137² = 18,769 for T_32. Counted to 2·10⁶: 41,719 stragglers for the 64 primes ≤ 311 (2.25 % of all composites up to 2·10⁶), 11,758 for the 128 primes ≤ 719 (0.64 %), 128,463 for the 15 primes ≤ 47 (6.9 %). Literature name: **y-rough** (y-sifted) numbers, counted by the Buchstab function Φ(x, y); the Mertens estimate ∏_{p≤y}(1 − 1/p) is 0.0968 for y = 311 (this is the density of rough numbers, primes included). For fixed y the stragglers never run out; they only thin. A straggler cannot be reached by firing the slotted circles alone.

## 4 The operator Euler product and an exact identity
Each Euler factor on σ = ½ is a circle (1 − p^(−1/2) e^(−it ln p))^(−1); putting the circle's i on slot k (i → L_k) gives F_p = (αI − βL_k)/(α²+β²), α = 1 − p^(−1/2) cos θ, β = p^(−1/2) sin θ, θ = t ln p (O(N) per prime). v_K = F_{p_K} … F_{p_1} e₀. Control: all primes on one slot = the complex partial Euler product to 2·10⁻¹⁶ (checked; remaining components exactly 0).

**Identity (exact, any assignment):** L skew and orthogonal ⇒ ⟨v, L v⟩ = 0 and |Lv| = |v|, so |F_p v| = |v|/√(α²+β²) = |v|/|1 − p^(−s)|. Hence **|v_K| = |ζ_K(½+it)| (the complex partial product) for every assignment of primes to imaginary slots**. The size of the path does not depend on the slot rule; only its direction does. A slot-0 prime (identity) breaks it (its factor is a real scalar with a different modulus) — that is the only source of the "ratio" and "drift" deviations in the JSON (seq0, char), and the huge z values there are artefacts of a null with sd = 0.

## 5 Tests of rules (E1 structure ↔ arithmetic; E2 direction), with a shuffle / random null
Rules: seq0, seq1, mod (slot = p mod N), char (slot = Σ 2^j[(d_j/p) = −1], d = −4, 8, −3, 5, −7, −11, 13: slot(pq) = slot(p) ⊕ slot(q), a homomorphism; primes in the trivial class land on slot 0), random.

**E1.** Phi correlation between the deform graph D[slot_i, slot_j] and three pre-declared pair features over 200 primes (17 … 1,223): both ≡ 3 mod 4 (the reciprocity sign), same residue mod 4, Legendre (p_i/p_j) = −1; null = 1,500 shuffles of the slots among the primes.
- **seq0 and seq1: no signal at any N** (|z| ≤ 1.6, p ≥ 0.11 for all 24 tests). The sequential rule is indistinguishable from a random permutation on every feature tried.
- mod and char: strong associations with "same residue mod 4" (z = −36 … −13) and, for mod, "reciprocity sign" (z = −36 … −22). These are BY CONSTRUCTION (bit 0 of the char slot is χ₋₄(p); slot mod 4 of the mod rule is p mod 4): they show the rules carry arithmetic in, they do not show the arithmetic selects them. The Legendre feature does not separate either (|z| ≤ 2.5, p ≥ 0.011, not significant over 12 tests).

**E2.** Operator Euler path, 400 primes (2 … 2,741), t ∈ {14.1347, 21.0220, 30, 50, 100}, statistics averaged over t, null = 300 random assignments to slots 1..N−1. By §4, size statistics are identical for all assignments. Direction statistics: the scalar-shadow error |v₀ − Re ζ_K|/|ζ_K| is 0.52–0.55 for random assignments (0 for the control) and every rule is within ±1.1 sd of random; the slot-occupancy entropy of v_K is within |z| ≤ 2.1 for every rule at N ≤ 64, and at N = 128 seq0 gives z = −3.5, seq1 −1.9: the sequential rules keep the heavy small primes on the low slots {1, 2, 3, …} (the octonion sub-block, where everything is alternative) so the path occupies fewer slots. That is a property of putting small primes on small indices, not evidence for the rule (32 entropy/shadow comparisons; one at |z| = 3.5, seq0 which also carries the slot-0 scalar).

## 6 Verdict
1. **Nothing here selects a slot rule.** The sequential rule (slot k ↔ k-th prime) is not distinguished from random permutations by structure (E1) or by the operator Euler path (E2); the magnitude of the path is slot-blind by an exact identity.
2. If the slot rule is to be arithmetic, **char** is the one candidate with an algebraic reason: its slot map is a homomorphism, so a product of primes sits on the XOR of its factors' slots (the multiplicative structure of Z maps onto the XOR group of T_N), it has no capacity limit (every prime gets a slot at any N), and its bits are the quadratic characters of the Fermat-side objects ζ_ℚ(√d)(s) = ζ(s)L(s, χ_d). It is a candidate, not a result: the first test that could favour it over random (E1, Legendre) did not.
3. Prime 2 on slot 0 (wiki 130 table) leaves 2 without a circle; seq1 or any rule that leaves slot 0 empty avoids that.
4. Not tested: any quantity that depends on the zeros' sedenion analogue (no zeros of a sedenion Euler product are defined here); the sedenion ring's non-associativity (we used the associative operator ring generated by the L_k, the "ring-circle" reading); a rule picked by the Fermat-side objects directly (L(Δ, s)).
