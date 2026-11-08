# The Heartbeat waveform, derived from the two sets (2026-11-08)

Code: `rhp/wave.py`, figure `rhp/atlas/wave.svg` (`python -m rhp.atlas.wave`), tests `test_wave_*`. Reference engine: `GenerationalLineage/engine/lineage.py::two_trees`
(same partition as `AbrikosovTree/engine/two_trees.py`). Labels: ESTABLISHED / OURS-CALCULATED / READING (not a theorem).

## 1. The two sets (the engine's own partition)
Every integer is Telperion (prime), Laurelin (composite) or Mingling ({0,1}); B(n)+R(n)+M(n)=1 exactly. Checked: engine and an independent sieve agree on
[0, 10⁵]: 9592 + 90407 + 2 = 100001. (`test_wave_partition_matches_engine_and_is_exact`)

## 2. A partition forces exactly opposite deviations (identity, ESTABLISHED-trivial)
Give each integer one unit of weight and split it: Telperion takes Λ(n) (ln p on a prime power, else 0), Laurelin takes 1 − Λ(n).
T(x) = Σ_{n≤x} Λ(n) = ψ(x);  L(x) = Σ_{2≤n≤x} (1−Λ(n)) = ⌊x⌋ − 1 − ψ(x);  T + L = ⌊x⌋ − 1.
Why Λ and not the prime count: Λ(n) ≠ 0 exactly when n is a prime power, i.e. exactly when the smallest and greatest prime factor coincide — the only composites on which
the two trees' orders of approach agree (earlier finding, 78,734 of 10⁶). With Λ the wave has no Chebyshev bias (π(x) − Li(x) carries −½Li(√x) from the prime squares).
Native variable u = ln x; divide by √x = e^{u/2} (this flare is what the ½ of the line does to the amplitude):
W_T(u) = (ψ(x) − x)/√x,  W_L(u) = (⌊x⌋ − 1 − ψ(x))/√x,  W_T + W_L = −({x}+1)/√x, |W_T + W_L| ≤ 2/√x.
Measured on u ∈ [6, 15.9]: max |W_T + W_L| = 0.098 ≤ bound 0.0996. The two trees are exact mirror images (figure a).

## 3. Expansion in zeros (von Mangoldt, ESTABLISHED)
ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²). With ρ = ½ + iγ:
W_T(u) = −Σ_ρ e^{iγu}/ρ − ln 2π e^{−u/2} − … = −2 Σ_{γ>0} Re[e^{iγu}/ρ] = −2 Σ_{γ>0} [½ cos γu + γ sin γu]/|ρ|².
One zero pair = two circles of equal radius 1/|ρ| turning opposite ways (e^{+iγu}/ρ and its conjugate). Their sum is the real tone −2cos(γu − arg ρ)/|ρ|,
a sine of amplitude 2/|ρ| shifted in phase by arg ρ = π/2 − arctan(1/2γ): the ½ is the small cosine admixture ½/γ of the sine.
**H(u) = −2 Σ_{γ>0} cos(γu − arg ρ_γ)/|ρ_γ|** is the foundational function; Telperion carries +H, Laurelin −H.

## 4. Measured
| quantity | value |
|---|---|
| correlation of W_T (from primes) with H, first 10 / 100 / 1000 / 3066 zeros | 0.771 / 0.931 / 0.984 / 0.993 |
| rms residual, same | 0.138 / 0.079 / 0.038 / 0.025 |
| rms of W_T on u ∈ [6,15.9] vs √(Σ_ρ 1/|ρ|²) = √(2+γ_E−ln4π) | 0.2173 vs 0.2149 |
| skew / excess kurtosis | −0.02 / −0.25 (flatter than Gaussian: sine-like tones, not a bell) |
| peak |W_T|, gust factor peak/rms | 0.795, 3.66 |
| envelope (500 tones, Hilbert): mean, sd/mean, max | 0.267, 0.48, 0.74; 9.5% of u below 0.1 (calm), 16% above 0.4 |

## 5. The 'wind' (READING)
- In the native variable the wave is stationary-looking (rms constant, 0.215). In the raw variable the excursion √x·W grows as e^{u/2}: a wind that increases without ever being periodic;
  its gusts are the beating of incommensurate tones (frequencies γ_k; independence assumed, not proved), its calms are envelope zeros.
- One factor p adds a step ln p to ψ at x = p: in W units ln p/√p. This is maximal at p = e² (ln p = 2), value 2/e = 0.736 — the engine's Mingling sits near e² = 7.389.
  The largest envelope seen (0.74) is of the same size as the largest single-factor step. OURS observation; no mechanism claimed.
- Fraction of u-time with |W_T| ≥ one factor's step: p=2: 2%, p=5–11: 0.02%, p=101: 3%, p=1009: 33%, p=10⁴: 68%, p=10⁶: 95%. Small factors need a rare gust (3.3σ at p=5, 11);
  by p ≈ 10³ the ordinary wind already equals a factor's step, so a single factor cannot be seen against the wind without the tones adding coherently at x = p.
- Resolving the prime p's jump needs bandwidth: fraction of the true jump ln p recovered with K tones — p=11 immediately; p=101: 0.34/0.59/0.94 at K=5/100/1000; p=307: 0.12/0.19/1.11 at K=5/100/1000 (Gibbs overshoot, 0.91 at 3066).
  Roughly T_K ≳ a few × p.
- 'Blow a factor off the tree' has no definition in the files yet; the numbers above are for 'wind ≥ the factor's own step'. The shakeable-leaf conjecture (semiprime N = pq, the point where it drops
  defines the factors via ζ_ℚ(√N)) is the other candidate and is not tested here.

## 6. Not shown
Nothing here is evidence about RH. The Fermat side (cusp form L(Δ,s) / ζ_ℚ(√N)) is not in this wave; the wave is the prime-side zero sum split between the trees. That the Telperion/Laurelin
interaction IS this wave is the working reading; what is proven is only §2–§3.
