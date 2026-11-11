# The Möbius / Mertens wave from the 25,000 zeros — 2026-11-11   :CALCULATED

Code `rhp/mertens_wave.py`; table `data/zetaprime_25000.tsv` (ζ'(ρ_k) = −i e^{−iθ(γ_k)} Z'(γ_k), 20 digits, `experiments/zetaprime_at_zeros.py`, checked against mpmath's direct ζ' to 4e−17); run `experiments/mertens_probe.py` → `.json`, `.log`; test `test_mertens_wave`.
The M-wave is the ψ-wave's zero-tones rescaled: ψ(x) = x − Σ x^ρ/ρ …, M(x) = Σ x^ρ/(ρ ζ'(ρ)) − 2 + (trivial-zero terms). Same carriers x^{½+iγ}; the amplitude is multiplied by 1/|ζ'(ρ)| and the phase turned by −arg ζ'(ρ).

## Results
1. **The slope table.** mean |ζ'(ρ)| 6.42; smallest 0.243 at zero 18,859; so the loudest M-tone is 4.1× its ψ-tone. log|ζ'| against log of the gap to the nearest zero: correlation +0.67 — the M-tones are loudest where two zeros crowd.
2. **Size of the wave.** RMS of M(x)/√x per decade 10²…10⁷: 0.189, 0.178, 0.177, 0.171, 0.168; the zero sum predicts √(2 Σ 1/|ρζ'(ρ)|²) = √0.0290 = 0.170 (partial sums 0.0283 at K = 100, 0.0289 at 1000, 0.0290 at 10⁴ and 25,000: converged). max |M|/√x to 10⁷ is 0.894 (the old conjecture |M| < √x holds in the range). Compare ψ: RMS 0.21.
3. **The two hemispheres are in antiphase.** corr(M(x)/√x, (ψ(x) − x)/√x) per decade: −0.756, −0.816, −0.826, −0.819, −0.815. Sign changes of M per decade: 30, 111, 348, 1161, 2214; ψ(x) = x crossings: 138, 539, 1650, 4963, 13118.
4. **Reconstruction.** RMS error of M_K/√x for x in 10³–3·10³: 0.019 (K = 300), 0.011, 0.006, 0.002, 0.001 (K = 25,000); in 10⁴–2·10⁴: 0.020, 0.011, 0.007, 0.004, 0.003.
5. **The μ ping** (unit-window jump of the explicit formula, classified to −1 / 0 / +1). 3-class accuracy with 3,000 tones: 1.00 to n ≈ 1000, 0.65 in 1000–3000; with 10,000 tones 1.00 to n = 3000, 0.62 in 3000–10,000; with 25,000 tones 1.00 to n = 3000, 0.956 in 3000–10,000, 0.56 in 10,000–20,000. Same resolution law as the ψ ping (K tones resolve n ≈ K/3 to 0.4 K).
6. **The convolution identity Σ_{n≤x} M_K(x/n) = 1** (the counting function and M are inverses). Residual with 25,000 tones: −0.001, −0.0001, 0.010, 0.0016 at x = 50, 200, 1000, 3000 (3000 tones: 0.005, 0.0005, 0.058, 0.051). Each zero-tone is annihilated by convolution with 1 because ζ(ρ) = 0: the tone lies in the null space of the counting side. This is a complementarity check the engine can run.
7. **Two questions from the same tones** (400 sampled n per block). 'Is n a hole' (ψ ping, prime power) balanced accuracy and μ 3-class accuracy with 25,000 tones: 2–200: 0.925 / 1.00; 200–1000: 0.970 / 1.00; 1000–3000: 0.980 / 1.00; 3000–10,000: 1.00 / 0.93; 10,000–20,000: 0.62 / 0.56. The hole ping resolves slightly further than the factor-content ping.

## Reading (not a result)
The ψ wave and the M wave are the two descriptions of one set of tones, the T side and the 1/ζ side, in antiphase (−0.82). The μ ping returns factor content (squarefreeness, parity of the number of prime factors), not primality; the ψ ping returns the holes. Cody's question 'is it a factor, not is it prime' is the μ ping's question.
