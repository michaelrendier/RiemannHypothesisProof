# The norm check at the tower's interfaces — 2026-11-11   :CALCULATED

Machinery: `rhp/normcheck.py` (`mul`, `defect`, `interface_samples`, `level_report`, `basis_pair_sums`, `zero_product_halves`); run `experiments/normcheck_probe.py` → `.json`; test `test_normcheck`.
Defect D(x,y) = |xy|²/(|x|²|y|²) − 1. Zero where the composition law holds; −1 at a zero divisor.

## Results
| level | max |D| over random x,y | LL | LU | UL | UU (random vectors confined to a half) |
|---|---|---|---|---|---|
| T_1, T_2, T_4, T_8 | ≈ 6e−16 (exact) | exact | exact | exact | exact |
| T_16 | 0.65 | **exact** | **exact** | **exact** | **exact** |
| T_32 | 0.74 | 0.61 | 0.58 | 0.65 | 0.70 |
| T_64 | 0.51 | 0.82 | 0.63 | 0.71 | 0.73 |
| T_128 | 0.39 | 0.49 | 0.51 | 0.55 | 0.58 |
| T_256 | 0.32 | 0.37 | 0.36 | 0.44 | 0.43 |
(L = indices < N/2, the old algebra; U = the new directions. Mean |D| ≈ 0.14, 0.15, 0.12, 0.09, 0.07 at T_16…T_256: the typical failure shrinks with N, the worst case does not vanish.)

1. **Law exact to T_8, broken at T_16** (Hurwitz, reproduced).
2. **At T_16 the break is a cross term.** Every product of two elements that each live in ONE half keeps the norm exactly (all four interface classes ≈ 1e−16). The law fails only when a factor has support in BOTH halves.
3. **Above T_16 the failure is inherited.** At T_32 the LL class (the lower half is T_16 itself) already fails, so each doubling carries the old failure up and adds its own.
4. **Exact, on sums of two units** (the objects of `zdtree`, |x|²|y|² = 4): |(e_a±e_b)(e_c±e_d)|² takes only {4} at T_8; only {0, 4, 8} at T_16 and T_32. Zeros and doublings come in EQUAL number: 336 / 336 at T_16, 5,040 / 5,040 at T_32. A zero divisor is the norm collapsing, its partner the norm doubling.
5. **All 336 zero products at T_16 are (L,U | L,U)**: both factors mix old and new directions. A zero divisor is born exactly where a unit from the old algebra is paired with one from the new.

## Reading (not a result)
The interface of a doubling is where the composition law either survives (T_8→ only inside a half) or not. Up to T_8 a mixed factor is harmless; at T_16 it is where every zero divisor lives; above T_16 the damage is already in the lower half. This matches the ZD tree: immortal zero divisors (old pairs never change status) and a lower half that is already broken.

Not done: the norm check at the Hurwitz/E8 integer lattices (norm-n counts are in `rhp/nshape.py`), and a per-prime use of the defect as an engine signal (deferred with the ValaQuenta engine).
