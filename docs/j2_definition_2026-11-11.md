# The J₂ definition, tested — 2026-11-11   :CALCULATED

Definition under test (wiki 51, wiki 60, Cody 2026-11-11): J₂ swaps Red and Blue; two strokes; the electron-spin 'two rotations to return'.
Code `rhp/j2.py`, run `experiments/j2_probe.py` → `.json`, test `test_j2`. Colour = bit 4 of the index (ptol.c: red {0–3, 8–11}, blue {4–7, 12–15}).

## A. As an operator on T_N (N = 8 … 64) — exact
1. **Every swap made by multiplication is a four-cycle.** J₂ := L₄ (left multiplication by the blue unit, the grading's own shift) swaps red ↔ blue, is an isometry, and has L₄² = −I, L₄⁴ = I: order 4 on vectors, order 2 on rays. Same for every slot (L_k² = −I for all k, all N tested).
2. **The orbit of the identity is the stroke cycle:** e₀ → e₄ → −e₀ → −e₄ → e₀. Two strokes give the sign flip (the half-turn, −1); four return. The involution is the statement about rays and colour; the 4-cycle is the operator. This is the spin structure (−1 after one 2π, +1 after two) and is forced, not chosen.
3. **A strict involution exists only as a bare index swap** P₄ (i → i⊕4, all signs +): order 2, swaps colour, but it is not a product by any unit, and it moves e₀ (the identity, '1') to e₄. Neither L₄ nor P₄ is an algebra (anti)automorphism (defects ≈ 18–55): no swap of red and blue can be, since blue·blue → red and the identity cannot move.
4. **Acting on the Euler circles** (slot rule): L₄ conjugates the circle on slot k to its time-reverse (L_k → −L_k) when L_k and L₄ anticommute ('anti'), leaves slot 4 alone ('axis'), and is neither for the zero-divisor pairs ('deform'). Slot 0 (identity) commutes with everything: a prime there is not reversed.

## B. The Klein group of the functional equation — the established 4-cycle that degenerates to a 2-cycle
s → 1−s and s → s̄ both hold: |ξ(s)−ξ(1−s)| ≤ 6e−32, ξ(s̄) = conj ξ(s) exactly; same for Λ(s, χ₋₄) (1e−31), the Fermat-side n = 2 character. Their product s → 1−s̄ fixes the line. Orbit {ρ, 1−ρ, ρ̄, 1−ρ̄} has 4 points off the line, 2 on it: all 25,000 zeros of the table have orbit size 2, none 4 (this restates the verified line count, no new evidence for RH).

## C. The arithmetic J₂ — the Galois involution
For a quadratic field ℚ(√d) the conjugation acts on the primes above p: **fixed** (p inert = irreducible = a hole of that level, (d/p) = −1) or **swapped** with its conjugate (split, (d/p) = +1). Measured holes share per level, first 20,000 odd primes: 0.501, 0.501, 0.501, 0.502, 0.501, 0.500, 0.501 (d = −4, 8, −3, 5, −7, −11, 13), i.e. ½ at every level (Chebotarev). All seven fields have class number 1, so split primes factor into elements.
**The `char` slot is the Frobenius element** of p in (ℤ/2)^m = Gal(ℚ(√d₁…√d_m)/ℚ): slot bit j = 1 ⇔ the j-th conjugation fixes p. Slot 0 = trivial Frobenius = p splits completely. Share of primes on slot 0: 0.499, 0.248, 0.123, 0.060, 0.029, 0.013, 0.005 for m = 1…7 (= 2⁻ᵐ).
ζ_K check, m = 1…4, ~20,000 primes each: the number of characters with χ(p) = +1 is 2^m when the Frobenius is trivial and 2^(m−1) otherwise, 0 failures (Euler factor (1−x)^(−2^m) against (1−x²)^(−2^(m−1))). **A prime that splits completely contributes the plain circle with multiplicity 2^m, not '1'.**

## D. The slot rules against the J₂ axis (e₄), shares of 2000 primes
| T_16 | zero | axis | anti (reversed) | deform |
|---|---|---|---|---|
| seq1 | 0 | .067 | .534 | .399 |
| mod | 0 | 0 | .503 | .496 |
| char | .053 | .063 | .503 | .381 |
| charx | 0 | .128 | .872 | 0 |
| random | 0 | .062 | .548 | .390 |
At T_32/64/128 the 'anti' share is 0.56/0.34/0.20 for `charx` against 0.32/0.19/0.10 for `char`, `seq1` and random. `charx` is the only rule that keeps slot 0 empty AND puts more primes on the reversible slots; `char` is indistinguishable from random there, and puts 5.3 % → 0.2 % of primes on slot 0.
Caution: `charx` uses the lower slots 1…N/2, and the deforming partners of e₄ lie mostly higher, so part of this is geometry of the slot list and not arithmetic. It is a property of the rule, not a selection of it.

## Reading (not a result)
- J₂ as defined is an involution on colour and rays and a 4-cycle as an operator; the heart's two strokes and the spin sign are the same statement.
- The arithmetic J₂ (conjugation) makes 'hole at level K' = 'fixed point of that level's involution', which is the irreducibles definition, and shows the trivial-Frobenius prime is a hole at no level yet still a circle.
- The crossing (B) depends on the Klein group alone, not on holes or slots.
