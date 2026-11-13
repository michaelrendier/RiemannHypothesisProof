# The flow of the factors past the prime holes; the finite-tone wave — 2026-11-12   :CALCULATED

Code `rhp/factor_flow.py` (was written over `rhp/flow.py` by mistake, see log 58); runs `experiments/flow_tests.py`, `experiments/fermat_mobius_growth.py` (→ `.json`); test `test_flow_tests`. Symbols follow the new glossary section (σ_ζ, Ω_ζΣ, A_fermat; α_fs is not used).
Question (Cody): do the factors flow past the primes; is it turbulent; are the primes rocks piercing the surface; what does the Laplacian over the Cartesian ζ say?

## F1 — the Laplacian over the Cartesian ζ (ESTABLISHED, checked)
Δ log|ζ(σ_ζ + it)| = 2π Σ_ρ δ(s − ρ) − 2π δ(s − 1). Flux of ∇log|ζ| through small circles, divided by 2π: **−1.00000 around s = 1; +1.00000 around each of zeros 1–5; 0.00000 around an ordinary point (½ + 20i); 2.0000 around a circle holding two zeros.** The zeros are point SOURCES of the potential log|ζ| (and vortices of arg ζ), the pole a sink: 'islands' in the strip are the zeros. Between them log|ζ| and arg ζ are harmonic: smooth, no extrema. (Euler product: for σ_ζ > 1 log ζ = Σ p^{−s}/… is regular, no sources: the prime currents alone give a laminar field; the sources appear only on continuation to σ_ζ ≤ 1.)

## F2 — islands against rigid points (spacing statistics; unit mean spacing)
| point set | spacing variance | P(s < 0.1) | number variance Σ²(5) | Σ²(10) |
|---|---|---|---|---|
| zeros 1–25,000 (unfolded by θ(t)/π) | 0.157 | 0.001 | 0.39 | 0.33 |
| zeros 20,001–25,000 | 0.160 | 0.000 | 0.40 | 0.33 |
| primes 10⁶–10⁷ (unfolded by li(p)) | 0.693 | 0 | 3.15 | 5.75 |
| primes 10⁷–10⁸ | 0.723 | 0 | 3.37 | 6.21 |
| primes 10⁸–3·10⁸ | 0.741 | 0 | 3.54 | 6.37 |
| Poisson control | 0.992 | 0.096 | 4.90 | 9.78 |
Zeros: repulsion (almost no small spacings), variance near the GUE 0.18 and number variance saturating at ≈ 0.35 (Berry's saturation) well below Poisson: a rigid, crystal-like row. Primes: spacing variance 0.69 → 0.74 and rising slowly toward 1 with x, number variance below Poisson but 10–20× the zeros' at L = 10: scattered nearly independently (the minimum normalised spacing 2/ln p ≈ 0.1–0.15 is arithmetic discreteness, not repulsion). **Primes are scattered islands; the zeros are the rigid structure.**
Depth field: h(n) = d(n) − 2 (non-trivial divisors) is zero exactly at the primes: the points no factor stream-line crosses; the composites lie under water to depth h(n).

## F3 — the Fermat finite-tone wave
p ≡ 1 (mod n); the 2g = (n−1)(n−2) Frobenius eigenvalues √p e^{iθ_j}; F(s) = p^{gs}∏(1 − ω_j p^{−s}) is real on σ_ζ = ½ (checked), with zeros at t = θ_j/ln p (mod 2π/ln p).
| curve | tones 2g | distinct angles | multiplicities | RMS of W_F (independent tones would give √2g) |
|---|---|---|---|---|
| n = 3, p = 7 | 2 | 2 | 1, 1 | 1.414 (√2) |
| n = 4, p = 17 | 6 | 2 | 3, 3 | 4.24 (√6 = 2.45) |
| n = 5, p = 31 | 12 | 4 | 3, 3, 3, 3 | 6.00 (√12 = 3.46) |
| n = 7, p = 29 | 30 | 8 | 6, 6, 3×6 | 11.23 (√30 = 5.48) |
- The standing wave has FAR fewer distinct tones than 2g and its zeros are **multiple** (order 3, 3, 6): the Jacobi sums repeat. The tones add coherently. The nodal lines of Re F, Im F at a zero of order m are 2m rays, not the single rung ⟂ critical line of the simple Riemann zero; for n = 3 (simple) the rung is horizontal to 0.04°; for the degenerate zeros the rungs cross at a different angle (≈ 70–78° from horizontal, extraction at multiple points is approximate).
- **Consequence for the inverse (Möbius) side:** the coefficients of 1/P(T) normalised by p^{k/2}: n = 3 (simple roots) bounded, envelope 1.02 at k ≤ 160; n = 4 (order 3) grows like k²: local exponents 1.94, 2.07, 2.01 (expected 2); n = 5 and 7 grow polynomially (exponents 1.7–2.8 against 2, and 3.6–4.6 against 5; the recurrence is ill-conditioned at large k). The Riemann zeros are simple to 25,000 (smallest |ζ′(ρ)| = 0.243), so its Möbius wave stays bounded (RMS 0.17); the Fermat curve's inverse wave is not.

## Reading (not a result)
- Turbulence: nothing here is a turbulent cascade. Each prime is an independent periodic current (CRT) and each zero-tone an independent rotation (Kronecker flow on a torus): superpositions, no scale-to-scale coupling, a LINE spectrum (zero-phase scan, docs/zerophase). The sources (zeros) are a rigid row; the islands (primes) are scattered.
- The rock-pile picture holds in two exact senses: Λ(n) is supported on prime powers (spikes over the mean level 1), and h(n) = d(n) − 2 vanishes exactly at primes; in the (σ_ζ, t) plane the 'islands' are the zeros (point sources).
- The Fermat finite-tone wave is the same kind of object (a sum of unit rotations, zeros on σ_ζ = ½) with two differences: finitely many tones, and degenerate zeros that make its inverse side grow polynomially.
Not done: the prime-spacing curve (deferred by Cody).
