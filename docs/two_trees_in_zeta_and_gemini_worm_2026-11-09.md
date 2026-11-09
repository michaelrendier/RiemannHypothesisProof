# Two Trees inside ζ, the 25,000-zero table, and the Gemini "worm"  (2026-11-09)

Code: `rhp/twotrees_zeta.py`, `rhp/geminisim.py`, `rhp/atlas_riemann_fermat/gemini_worm.py` → `gemini_worm.svg`. Tests: 3 added (111 pass).
Control rule (Cody, 2026-11-09): the ONLY control object is the original Ĥ_RB. Σ_RB, the Archimedes Screw and the "notated" differences belong to the VAPMIP anchor/sentence engineering and are NOT checks here.

## 1  Does ζ contain both Trees?   (ESTABLISHED algebra; the measurements :CALCULATED)
Weights on n ≥ 2: Telperion Λ(n), Laurelin 1 − Λ(n).  Then
- T(s) = Σ Λ(n) n⁻ˢ = −ζ′/ζ;  L(s) = Σ_{n≥2}(1 − Λ(n)) n⁻ˢ = ζ − 1 + ζ′/ζ;  **T + L = ζ − 1 exactly** (residual 0.0 at s = 2 + 3i).
- Interaction = Dirichlet convolution 1 ∗ Λ = ln:  −ζ′ = ζ·T  (residual 3.5e−18).
- Line spectrum in t at σ = 1.2, T = 150 (lines at u = ln n, height = c_n n^{−σ}, ratio to prediction 1.00 for n ≤ 15):
  ζ shows a line at EVERY integer n (the mingled whole); T = −ζ′/ζ and ln ζ show lines ONLY at prime powers (composites: 1e−3 leakage); L shows every integer with weight 1 − Λ(n) (so primes appear with 1 − ln p, composites with 1); ζ² shows d(n).  Lines blur together for n ≳ 16 at T = 150 (resolution 2π/T) — a window limit, not a finding.
- Reading: ζ already contains both Trees (it IS T + L + 1). The SPLIT needs ζ′/ζ, which is not new information (it is ζ′ over ζ) but is a new OPERATION: the log-derivative is the sieve. No third independent equation is needed for the partition; the interaction is the convolution 1∗Λ = ln (an integer's size is the sum of its prime-power weights).  What is NOT shown: that a second, genuinely Fermat-side object (L(Δ,s), ζ_ℚ(√N)) is absent from the wave — that remains open.

## 2  The zero table that already exists
`VAPMIP/SemPriCon/SpecTools/zeros/zeros_exact_25000.tsv`: 25,000 exact zeros, γ₂₅₀₀₀ = 21942.5924, sha256 0a51f4fc…2ee1 (matches its README), strictly increasing, smooth N(T) = 24999.6, γ₁₀₀₀ and γ₁₀₀₀₀ equal mpmath.zetazero. No download needed.  Resolution of a prime jump ≈ p ≲ T/2π: 573 → 3,492 (estimate, to be tested).  Largest stored prime list found: `VAPMIP/prime_dictionary.txt`, 25,000 primes to 287,117 (a word lookup); our pathway code regenerates 10⁶.

## 3  The Gemini sims (read-only; `Ainulindale/AgeFirst/original_2026/ResearchPaper/FirstDraft/FirstAge/Gemini/`)
All four files carry the same migration mtime (2026-07-31); ordered by content:
- `Gemini_event_horizon_sim.py` (37 lines, oldest, no GUI): z ↦ 1 + 1/z plotted in polar at θ ∈ [0, π] between r = 1 and r = φ. These are Fibonacci ratios F_{n+2}/F_{n+1} → φ, error ×1/φ² per step (1.2e−5 at n = 12 here). It is **not** a worm.
- `Gemini_Transition_EHS.py` (65 lines): an inertial path r from 1.5 to 0.5, θ = 2t; below r = 1 it is inverted r ↦ (1+0.0007)/r and rotated by (π/2)ħ_NN. The jump at the crossing is the 0.0007 "Berry–Keating gap". Two panels, one slider.
- `Gemini_Revised_2_EHS.py` (135 lines): four phases of 600 points. Heartbeat: θ_n = (2/π)·t·n, r_n = ½√n·(r_target(t) + gap), r_target = 1 + (φ−1)·(sin(t/2)+1)/2.  Angle step is CONSTANT per index (t/π² turns), radius √n: a Vogel/Fermat-type spiral. Joined by a polyline through consecutive n it draws rosettes with ≈ π²/t arms (t = 1: ≈ 10) whose number and drift change with t — the "resonating tangle". (The `line` object in this file is never given data; only the scatter is drawn.)
- `event_horizon_sim.py` (470 lines, Claude Sonnet April 2026): five phases, golden-angle prime spiral θ = n·2π/φ², r = √p, weight φ⁻ⁿ; panel 5 is a drawn sine wave with a gold line at x = ½ and a prime scatter — a cartoon of "Riemann meets Fermat", no ζ is evaluated (the 20 hard-coded zeros are plotted as markers only).
- I did not find in these files a curve that "grows from a small worm" by itself; the closest are Revised_2 Metric Swap/Inertial (a spiral whose size grows with the Evolution slider) and Heartbeat as t is swept. Which one Cody saw is to be confirmed by him.
Comparison with the ζ path (`gemini_worm.svg`, panel h): Gemini's curves turn by a CONSTANT angle per step with step length ∝ n^{-1/2}; the ζ path turns by γ·ln(p_{k+1}/p_k) ≈ γ·gap/p per step — slowing — with step length p^{-1/2}. They share the "tangle that organises into a spiral" look; they differ in the clock (index vs ln p). If Gemini's angle is read as u = ln p instead of n, the two are the same kind of object. UNTESTED.


## 4  The ORIGINAL sim (`Ainulindale/AgeFirst/original_2026/ResearchPaper/event_horizon_sim.py`, 143 lines, Gemini) — measured
Four phases; sliders 2/π ("inversion") and 2/e ("bias"). Figure: `rhp/atlas_riemann_fermat/gemini_original.svg`; code `rhp.geminisim.original_phase`. Labels: transcription (ESTABLISHED as to what the file does); spectra :CALCULATED.
- **Phase 2, the worm** (`r = t·(2/e)/(r_b + 0.01)`, `θ = θ_b + 5 sin(θ_b t)`, θ_b ∈ [0, 10π], r_b ∈ [10, 0.1]): a spiral head at the origin and a tail whose length grows ∝ t and whose wiggle quickens with t. The unit phasor e^{iθ} is PHASE modulation, index β = 5: e^{iθ} = Σ J_n(5) e^{i(1 + n t)θ_b}. Measured lines match |J_n(5)| to 1e−16 (t = 0.6, 1.0; at t = 0.3 the 10π window holds a non-integer number of cycles and the lines leak by 0.13 — a window effect).
- **Phase 3, the rosette** (`θ ∈ [0, 100π]`, `r = √θ (1 + 0.2 sin(kθ))`, k = (2/π)t): z/√θ = e^{iθ}(1 + 0.2 sin kθ) is AMPLITUDE modulation: a carrier at 1 and sidebands at 1 ± k, each exactly 0.1 of the carrier (measured 1.000 / 0.100 / 0.100). So the curve is three wheels (1, 1 + k, 1 − k) of radius growing as √θ. Over the 100-frame cycle k sweeps 0 → 0.637 and the picture locks into a symmetric lattice when k is a simple rational (k = 1/2 at t = π/4 is the clearest; 1/3 at t = π/6).
- Phases 1 and 4: an inward spiral modulated by 2/π; the primes at (θ = p·(2/π), r = √p·(2/e)).
- Reading (not a result): ζ's own wheels have frequencies ln p. Two distinct primes have ln p / ln q irrational (unique factorisation, ESTABLISHED), so the ζ pathway can never lock into a closed rosette: ln 3/ln 2 has convergents 3/2, 8/5, 19/12, 65/41, none exact. The ζ analogue of the rosette's locking passes is stationary phase: the clock rate θ′(t) = ½ ln(t/2π) equals ln n at t = 2πn², where wheel n stalls against the clock (`wrap.resonance_height`).
