# Filling the gaps before the standing-wave module — step 1 (tone bank, ledger, Guinand–Weil, Weil counts) — 2026-11-12   :CALCULATED

Order agreed with Cody: the gaps one at a time, each building on the last, the five analogous tests folded in. Plan and state:
| step | gap / test | state |
|---|---|---|
| 1 | gap 1 (one object) + gap 4 (executable ledger) + analogous test 1 (Guinand–Weil two-sided identity) | **done here** |
| 2 | analogous test 2: zero-phase scan on the Fermat cubic's a_p (zeros of its L-function from Fermat-side data) | next |
| 3 | analogous test 3: Chebyshev bias, L(s, χ₋₄) zeros against Fermat's two-squares primes | |
| 4 | analogous test 4: linear-recurrence rank (Berlekamp–Massey / Prony), finite (Fermat) against growing (Riemann) | |
| 5 | gap 2: the second coordinate of the two-hemisphere field, defined for both sides | |
| 6 | gap 3 + analogous test 5: zero dynamics (de Bruijn–Newman flow), stability of the critical-line row | |
| 7 | gap 5: violation detectors for RH (what a zero off the line would do to each check) | |
| 8 | the standing-wave module | |

## What was built
- `rhp/tonebank.py`: `ToneBank.zeta(K)` (γ_k, ζ'(ρ_k), ψ and M weights; `psi_wave`, `mobius_wave`, `zero_sum`, `resolution`) and `FermatBank(n, p)` (Frobenius angles with multiplicities, `wave`, `affine_count(k)`). Same attribute names on both sides.
- `rhp/gf.py`: exact finite fields GF(p^k) by log tables and the exact count of xⁿ + yⁿ = 1 over them.
- `rhp/explicit.py`: `guinand_weil(γ, Λ, r₀, σ)` returns both sides and the residual.
- `rhp/ledger.py` (`python -m rhp.ledger` → `experiments/ledger.json`): 14 identities, each with a residual and a tolerance.

## Results
1. **Guinand–Weil, zeros against primes** (25,000 zeros, Λ(n) to 3·10⁶; Gaussian pair, shifted to r₀, width σ): residual **9·10⁻¹⁴** (r₀ = 14.13, σ = 1), 1·10⁻¹³ (σ = 0.7), 1.2·10⁻¹³ (r₀ = 30), 5·10⁻¹³ (r₀ = 100, σ = 0.8), 4.5·10⁻¹² (r₀ = 1000), 9.8·10⁻¹¹ (r₀ = 10,000), 1.25·10⁻¹⁰ (r₀ = 15,000). At r₀ = 14.13 the zero sum is 2.0000000001 (the first zero, counted at ±γ) and the right side is the sum of the Γ term (+1.5581), the −g(0)ln π term (−0.9134, i.e. lnπ term +0.9134 subtracted) and the prime sum (+1.3553): the two descriptions give the same number to 13 digits. The zero table reaches γ = 21,942.59.
2. **Weil explicit formula on the Fermat side, exact**: brute-force counts of xⁿ + yⁿ = 1 over GF(p^k) agree with p^k − Σ(−J_p)^k to 4·10⁻¹³ for (n, p) = (3, 7), (4, 5), (3, 13), (5, 11) and k up to 4 (3, 7: 6, 60, 321, 2328; 4, 5: 8, 40, 104, 664). The finite tone sum reproduces the point count at every power.
3. **The ledger**: 14 of 14 identities pass in 5.7 s: Ω four routes 2.5·10⁻³⁹; closure 0; prime-spacing crossing 7·10⁻⁴¹; digits per step 2.5·10⁻¹⁵; Klein group 1.3·10⁻³¹; Guinand–Weil ×3 (9·10⁻¹⁴, 4.5·10⁻¹², 1.25·10⁻¹⁰); convolution Σ M_K(x/n) = 1 at 1.1·10⁻⁴; Möbius explicit-formula RMS 1.1·10⁻³; ξ/1/ξ nodal sign agreement 0; Weil counts 4.6·10⁻¹³; |J| = √p 8·10⁻¹⁴; F real on the line 1.5·10⁻¹⁵. Table: `docs/ledger_table.md`.

## Reading
- Guinand–Weil is the two-sided identity in its sharpest form: the SAME number computed once from the zeros and once from the primes; the Weil count is its finite-field twin, exact at every power. These are the two 'data sets' the standing wave must reproduce, and any engine that claims to carry both sides has to pass them.
- The ledger is the executable form of the 'says the same thing in every maths' rule. It does not yet contain a 0_RB residual: no concrete boundary check has been written down for 0_RB.

## Step 2 — zero-phase scan of the Fermat cubic's a_p (x³+y³=z³, conductor 27)   :CALCULATED
Exact zeros of L(E,1+it) from the AFE (mpmath, t<30): 21, first 4.043044, 6.048935, 8.21765, 9.429199, 10.908728.
a_p-only scan A_E(t) at X=2·10⁶ (resolution ≈2π/ln X ≈ 0.43): peaks at 12 of the first 12 zeros within 0.082 (4 within 0.03); heights 0.71–0.97 (mean 0.86).
Also 15 extra peaks (9.51, 10.2, 12.04, …) with no zero beside them: side-lobes of a ~0.4-wide window, not matched. Random-sign control: median 0.60, max 1.67, 34 peaks — it is NOT cleaner than the real scan by peak count, so the evidence is the zero-matched heights and positions, not peak count. Verdict: the Fermat side's zeros are readable from its a_p alone, at the coarse resolution; extra-peak rate needs a window/apodisation before it can be called a detector.
