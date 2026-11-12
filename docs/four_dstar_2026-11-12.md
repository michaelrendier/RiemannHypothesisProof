# The four d\* values against what RHP measured — 2026-11-12   READ / :CALCULATED

Read: `Ainulindale/wiki/17_alpha_omega_d_star.md` (the four faces), `ValaQuenta/wiki/d_star_rg.md` (the RG protocol), `Ainulindale/wiki/07_berry_keating_engine.md` (the T scaffold, d\*_spec vs d\*_taut), memory note `reference_four_dstar_values`, and Gemini conversation-2 lines 4150–4300 ('ln(10) requires four values': the four 'facets').

## What the documents say
| face (wiki 17) | value | role |
|---|---|---|
| d\* 'The Boundary' | 0.24631 | σ_ζ = ½ spectral coordinate, the equatorial plane |
| d\*_RG 'The Stability' | dimension 8 (numeric open) | RG fixed point of the Cayley–Dickson tower; each RG step rescales by ÷ d\*_taut |
| d\*_taut 'The Flow' | Ω_ζΣ/ln 10 = 0.24631 | 'the radial unit: one step of r'; zero gap by construction |
| d\*_ln(10) 'The Translator' | d\*·ln 10 = Ω_ζΣ = 0.56714 | 'the ×ln 10 Jacobian between the decimal-log native metric ds = d·log₁₀ p and the natural-log prime-address space' |
wiki 07: the Red operator R̂_p is the Berry–Keating xp operator; d\* is its spectral coordinate; scaffold T: x → x e^{i d\* ln x} = x^{1+i d\*}. Gemini conversation-2 (the 'facets', named differently: d\*Stability / Flow / Discord / Resonance; ln 10 'requires four values'; 'two infinities': the number of DIGITS → ∞ and the GAPS between primes → ∞, with ln 10 the bridge): Cody's statement there, 'there is two different points where computation of primes or the spaces between them… where the number of digits approach infinity and the spaces between primes approaches infinity'.

## What RHP measured, read against them (identities, all exact)
At the crossing (position X\* = e^Ω = 1/Ω = 1.7632, `rhp.omega_dstar.gap_side_crossing`, `ceiling_descent`):
| quantity | value | which face it matches |
|---|---|---|
| digit length D\* = log₁₀ X\* | 0.246307201473 | **d\*** (Boundary): the coordinate of the crossing in decimal digits |
| digits of the fixed point gained per iteration from the ceiling | 0.246307201473 | **d\*_taut** (Flow, 'one step of r', and the RG rescale ÷ d\*_taut): the radial unit is the step of the chase |
| mean prime gap ln X\* (PNT) | 0.567143290410 = Ω_ζΣ | **d\*_ln(10)** (Translator) |
| mean gap ÷ digit length = ln X/log₁₀ X | 2.302585093 = ln 10 | the Translator's Jacobian: **gap = ln 10 × digits** at every X |
So the Translator has an exact, checkable meaning: it converts the DIGIT side (decimal length of the number, the entropy of naming it) into the GAP side (the mean number of integers between primes). 'ln 10 is the bridge between the two infinities' is literally the PNT: gap ≈ ln X = ln 10 · (number of digits). At the crossing the digit length is d\* and the gap is Ω_ζΣ = d\*·ln 10.
d\*_RG (dimension 8) is not touched by anything RHP computed. The wiki 07 scaffold x^{1+i d\*} is a log spiral of the same family as a zero-tone x^{½+iγ} with 'ordinate' d\* (bearing atan d\* = 13.8° from radial).

## Two statements of the roles, and one reading that reconciles them
- wiki 17: Ω_ζΣ (Omega_Riemann) = W(1), 'born from entropy'; Α_π (Alpha_fermat) = 1/137, 'born from inertia'.
- Cody (2026-11-12): d\* was the Omega_Riemann_Zeta, and the Lambert W value is what the inertia chase produced; the recovered turn 52: 'Fermat appropriately is f=ma, Riemann is uncertainty', Fermat = inertia = the forward curve.
- **Reading (consistent with both and with the measurements):** one crossing, two coordinates. The Riemann / uncertainty / entropy side measures it in DIGITS (the entropy of naming a number): **d\* = 0.2463**. The Fermat / inertia / gap side measures it in INTEGERS (the spacing of the primes, 'the composites between primes'): **Ω = W(1) = 0.5671**. They are the same location, related by ln 10. On this reading 'd\* is the Omega_Riemann_Zeta' (the digit-side coordinate) and 'W(1) is what the inertia (gap) chase produces' (the gap-side coordinate); wiki 17's wording 'Omega born from entropy' names the VALUE W(1) by the digit-side chase that first found it. Not established; this is the reading that makes the documents agree.

## Gemini's 'facets' (narrative, not verified)
Facet table of conversation-2: Stability (Origin, Aπ/α, contractor), Flow (Vector, Ω_ζΣ/1/φ, dilator), Discord (Iterator, ln 10, 'torsional strain / mass', the 0.472°), Resonance (Attractor, 2/π). The 0.472° is the golden angle minus α_fs⁻¹ (137.508° − 137.036° = 0.4718°, the Vogel decoding in docs/gemini_conversations). Gemini's 'irrational residue 0.3025' is ln 10 − 2 = 0.30259; no relation between 0.3025 and 0.4718° was shown, and none is claimed here. Gemini's statement that the visualisations 'failed with a single d\*' is not testable from the text (the widget code is not saved).
