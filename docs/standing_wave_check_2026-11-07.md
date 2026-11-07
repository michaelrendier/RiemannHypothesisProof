# The pasted "standing wave" summation, checked claim by claim (2026-11-07)

*Assistant note: written by an AI assistant (Claude) with the author. The pasted text (from Gemini) is a stale summary of earlier conversations; it is used only as a pointer to where the standing wave lives. Code: `rhp/standing.py`; tests: `tests/test_unsmap_chladni_coil.py` (last five). Labels: ESTABLISHED, OURS-CALCULATED, FALSE (tested), UNSUPPORTED.*

## Where the standing wave is in the author's own files
- `ValaQuenta/modules/bao_mass_gap/maths.py` (header): "plus one standing wave per non-trivial zero" — the explicit formula ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ln(1−x⁻²); residue Δ = Ω_ζΣ − D*·ln 10 "absorbed by no standing wave". This is the zero-tone sum in u = ln x that `rhp/diverge.py` and `rhp/coil.py` already build (conical helices, flare e^{u/2}).
- `ValaQuenta/modules/archimedes_screw/maths.py`: same tones (`tone`, `tone_sum`, `chebyshev_psi_explicit`).
- `ValaQuenta/modules/inversion/maths.py` and `engine/units.py`: J_N: (r,θ) → (1/r, θ+π/2), φ = π/2 derived as 2π/4 (Hurwitz: four normed division algebras; π from Basel).
- `ValaQuenta/modules/derivation_chain/`: Tier 1 "Riemann = Fermat (R̂† = B̂)".
- No file in the searched repos contains an area-π integral of a standing wave; the only π there is φ = 2π/4.

## Claims

| # | Claim | Status | Evidence |
|---|---|---|---|
| C1 | J(z) = i/z̄ acts as (r,θ) → (1/r, θ+π/2); J⁴ = id | ESTABLISHED, tested | J² = −z, J⁴ = id; matches ValaQuenta `inversion_transform` to 1.5e-15 |
| C2 | "the fixed set of J is r = 1" | FALSE as stated | J has no fixed point: \|J(z) − z\| = √2 on r = 1. The circle is INVARIANT (a quarter turn per stroke). The engine's own comment is right: only the radius closes after two strokes; the wiki headline "involution whose fixed point r = 1" is loose and is recorded, not edited |
| C3 | J governs ζ / ξ | UNSUPPORTED | \|ξ(s′)\|/\|ξ(s)\| under J is 10⁶–10¹² (not 1). The mirror w ↦ 1/w̄ is an exact symmetry (ratio 1 to 1e-9); J = (quarter turn) ∘ (mirror), and the quarter turn is not a symmetry of ξ |
| C4 | Ψ = R(r)Θ(θ), R = AJ₀(k ln r) + BY₀(k ln r), "node at r = 1" | FALSE | J₀(0) = 1 (antinode); Y₀ diverges at r = 1; cos(k ln r) = 1 at r = 1 |
| C5 | ∫∫ over the unit disc gives π | FALSE for Ψ | the integral written is ∫∫ r dr dθ = π (the disc's area, no Ψ in it). For Ψ = cos(k ln r)cos(mθ+π/4) the energy is (π/4)(1 + 1/(1+k²)): 0.7893 at k = 14.1347 (closed form and quadrature agree to 3e-9) |
| C6 | Ψ oscillates in ln r with k = γ | NOT what ξ does | ξ is harmonic. Across the equator the profile is GROWTH, not oscillation: ∂σ ln\|ξ\| = θ′(t) + R(s) exactly (below). Oscillation is along t |
| C7 | "Eichler–Shimura established by Wiles" | FALSE attribution | Eichler–Shimura is 1950s; Wiles (with Taylor) proved modularity of semistable curves |
| C8 | the involution of the modular side | ESTABLISHED background | the Fricke / Atkin–Lehner map τ ↦ −1/(Nτ) gives Λ(E,s) = εΛ(E,2−s); the PSL(2,ℤ) fundamental domain has hyperbolic area π/3, not π |

## The exact transverse rate (what replaces C6)
From ξ = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s):

∂σ ln\|ξ\| = Re[1/s + 1/(s−1)] − ½ ln π + ½ Re ψ(s/2) + Re ζ′/ζ(s) = θ′(t) + R(s),

R = Re[1/s + 1/(s−1)] + ½[Re ψ(s/2) − Re ψ(¼+it/2)] + Re ζ′/ζ(s), with θ′(t) = ½ Re ψ(¼+it/2) − ½ ln π.
Tested: identity defect ≤ 1.2e-10 (finite difference); the first two terms of R are 4e-4 and 4e-5 at t = 30–50; only Re ζ′/ζ is not small. For σ > 1 its mean over t is 0 (almost-periodic Dirichlet series with no constant term; measured 0.023 over [40,60], 0.0001 over [40,400] at σ = 1.5). Between ½ and 1 I claim no mean: measured mean growth over t ∈ [40,60] is 1.124, 1.082, 1.074 at σ = 0.55, 0.7, 0.9 against mean θ′ = 1.034 (finite window).
So the growth rate across the equator is θ′(t) = ½ ln(t/2π) up to Re ζ′/ζ: the 'transverse wave number' is real-exponential, the field is a standing wave along t with cosh/sinh profile across σ (Re ξ even, Im ξ odd).

## What survives as a usable statement
J_N is a rotary reflection (mirror through the equator composed with a quarter turn); the mirror part is the functional equation, the quarter turn is an extra structure the author derives from Hurwitz and Basel. Whether the quarter-turn has a counterpart on the ξ side is not shown here.
