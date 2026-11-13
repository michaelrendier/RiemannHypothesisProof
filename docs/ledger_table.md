| # | check | side | residual | tol | ok |
|---|---|---|---|---|---|
| 1 | Ω_ζΣ: four independent routes agree (closed form, contraction, Newton, bisection on the digit curves) | horizon | 2.50e-39 | 1e-35 | ✔ |
| 2 | closure d*·ln 10 = Ω_ζΣ (the Translator) | horizon | 0.00e+00 | 1e-28 | ✔ |
| 3 | prime-spacing crossing: mean gap ln X = 1/X at X = 1/Ω; gap = ln 10 × digits | gap | 6.89e-41 | 1e-30 | ✔ |
| 4 | digits of the fixed point gained per iteration from the ceiling = d* (to 1e-8) | horizon | 2.49e-15 | 1e-08 | ✔ |
| 5 | functional equation and Schwarz reflection (Klein group): ξ(s) = ξ(1−s), ξ(s̄) = conj ξ(s); Λ(s,χ₋₄) | riemann+fermat | 1.25e-31 | 1e-28 | ✔ |
| 6 | Guinand–Weil: Σ_ρ h(γ) (25,000 zeros) = smooth terms − 2ΣΛ(n)n^{−½}g(ln n) (primes), window r₀ = 14.13, σ = 1 | boundary | 8.88e-14 | 1e-12 | ✔ |
| 7 | Guinand–Weil, window r₀ = 1000, σ = 1 | boundary | 4.49e-12 | 1e-10 | ✔ |
| 8 | Guinand–Weil, window r₀ = 15,000, σ = 0.8 (near the top of the zero table, γ_max = 21,943) | boundary | 1.25e-10 | 1e-09 | ✔ |
| 9 | convolution identity Σ_{n ≤ x} M_K(x/n) = 1 (zero-tones are annihilated by convolution with 1), x = 200, K = 25,000 | riemann (1/ζ side) | 1.09e-04 | 1e-03 | ✔ |
| 10 | Möbius explicit formula: RMS of (M_K − M)/√x for x ∈ [10³, 3·10³] with K = 25,000 zero-tones | riemann (1/ζ side) | 1.07e-03 | 2e-03 | ✔ |
| 11 | ξ and 1/ξ share their nodal lines: sign agreement of Re/Im on the strip grid (1 − fraction) | riemann (both hemispheres) | 0.00e+00 | 1e-12 | ✔ |
| 12 | Weil explicit formula on exact counts: #{xⁿ + yⁿ = 1 over GF(p^k)} = p^k − Σ(−J)^k, (n,p) ∈ {(3,7),(4,5),(3,13),(5,11)}, k ≤ 3 (brute force over the field) | fermat | 4.55e-13 | 1e-09 | ✔ |
| 13 | Fermat tones have modulus √p: max | |J| − √p | over the 2g generic Jacobi sums, (n, p) = (5, 31), (7, 29) | fermat | 8.17e-14 | 1e-06 | ✔ |
| 14 | F(s) is real on the critical line (finite-tone completed function): max |Im F| / max |Re F|, (n, p) = (5, 31) | fermat | 1.53e-15 | 1e-06 | ✔ |
