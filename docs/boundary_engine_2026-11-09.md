# Boundary Condition Engine — patterns, transitions, and how the ζ spiral turns  (2026-11-09)

Code: `rhp/boundary.py` (patterns, transition mappings, diagnostics), `rhp/boundary_engine.py` (viewer), `rhp/viewport.py` (from the nes-viewport skill, unchanged).
Outputs: `rhp/atlas_riemann_fermat/boundary_engine/` (poster, strips, 15 transition strips, 4 GIFs). Tests: 4 added (117 pass).
Run: `python -m rhp.boundary_engine` (window) · `--export DIR` · `--selftest`.
Sources: the five Gemini files (original `ResearchPaper/event_horizon_sim.py`; `Gemini/event_horizon_sim.py` (Claude); `Gemini_event_horizon_sim.py`; `Gemini_Revised_2_EHS.py`; `Gemini_Transition_EHS.py`) — read, not edited.

## Design (the menu is a layer over the viewport)
Map layer = a 5 × 5 grid of tiles (rows = source files, columns = phases; 19 patterns: O1–O4, C1–C5, R1–R4, G1, X1, W1, Z1–Z3), content `cw × ch`, scrolled by a `Viewport` (3 × 2 tiles, screen-snap archetype). Overlay layer = scrollbars (derived: `vbar`/`hbar`), the menu bar (generated from `MENU`), the HUD; `App.overlay()` is a pure function of state. Modes MAP → STAGE → TRANSITION. `selftest()` drives the menu by key presses and checks the map layer's recorded data against the formulas: 52 checks (cursor always inside the viewport after every key, `assert_viewport_sound` after every move with the offset restored, thumb follows offset, tile = pattern, play moves the map each tick and τ rises, pause freezes it, parameter keys change the map, the worm's turn-around law, and all three transition mappings run α 0 → 1 through the menu). It found two real defects on the way (my harness left the viewport displaced by the destructive soundness sweep; Home/End moved the viewport but not the cursor) — both fixed.

## Patterns isolated (static and in motion)
Worm (original ph.2) = phase modulation β = 5; rosette (original ph.3) = amplitude modulation 0.2 (see the 2026-11-09 doc §4). New here: the **wheel family** `wheel_draw(Th, β, ω, d, k, m, a)` contains both exactly (checked point for point at 5000 points) and everything between.

## Transition mappings (`boundary.transition`)
- `polar` — Gemini's TransitionManager: index-matched LERP of θ and r.  `cartesian` — LERP of the points.  `wheel` — LERP of the wheel-family STATE (only for O2, O3, W1; otherwise falls back to polar and says so).
- O2 → O3 by `wheel`: structured, concentric intermediate states; by `polar`: a tangle at α = 0.25 (see `transition_O2_O3_*.png`).
- Turn-arounds along O2 → O3 (wheel), τ = 0.6: 6 angular reversals at α = 0, rising to ≈ 160 near α = 0.3–0.4 (Th grows 10π → 100π so many more wobble cycles fit), then ZERO from α ≈ 0.45. The cut-off is exact: reversals need β·ω > 1, and β(α)·ω(α) = 5τ(1−α)², which is 1 at α = 1 − 1/√(5τ) = 0.423 for τ = 0.6.  The transition is where the turn-arounds are created and destroyed.

## What makes the worm turn around (EXACT)
θ = θ_b + β sin(ωθ_b) ⇒ dθ/dθ_b = 1 + βω cos(ωθ_b). It runs BACKWARDS in angle exactly when βω > 1, for the fraction arccos(1/(βω))/π of every cycle; two reversals per cycle. For the original (β = 5, ω = τ) the worm never turns back for τ ≤ 0.2 and has 2·(5τ) reversals above it (0 at τ = 0.15, 6 at τ = 0.6, 10 at τ = 1.0 — measured).

## How the ζ Euler path turns (log P_K at s = ½ + it; measured, 20,000 primes)
1. **Tangent rotation is ONE-sense.** The k-th step is −ln(1 − p^{-s}) = p^{-s}(1 + p^{-s}/2 + …): direction −t ln p, so consecutive steps turn by −t·ln(p_{k+1}/p_k). That is one sense for every k once t·ln(p_{k+1}/p_k) < π. Last prime where the sense flips: 2 at t = 3, 31 at t = γ₁, 131 at t = 30, 547 at t = 100. After that the tangent never reverses (curvature reversals: 6 at γ₁, all among the first primes).
2. **The ζ phase-modulation index decays.** Step angle = −t ln p + ½ p^{-½} sin(…): measured rms of the residual 0.00577 / 0.00184 / 0.00041 on p ∈ [10³,10⁴) / [10⁴,10⁵) / [10⁵,3·10⁶) against the analytic ½p^{-½}/√2 = 0.00549 / 0.00173 / 0.00037, at t = γ₁ and 30 alike (within 5 %). The worm's index is the constant 5. So the worm's mechanism (βω > 1) cannot operate in ζ beyond the first few primes (β_ζ·ω < ½).
3. **Loops.** Winding −20.2 / −38.5 / −104.6 at t = γ₁ / 30 / 100 (one tangent turn per 2π of t ln p: loops are spaced geometrically, ratio e^{2π/t} in p). Angular reversals about the origin (retrograde): 47 / 91 / 244, i.e. 2.33 per turn, constant across t; the worm has 2.0 per cycle. That is close, but any loop displaced from the origin has 2 per turn — generic geometry, NOT evidence of a shared mechanism. In ζ they come from the fluctuation of the partial sums about the smooth spiral (|z|/R(p) spreads 0.4 – 2.3, R = p^{½}/(|½−it| ln p)).
4. A "tone share" metric for the tangent rate (fraction of variance in one Fourier line) did NOT discriminate (worm 0.58 at τ = 0.15 but 0.02 at τ ≥ 0.6; ζ 0.10–0.27) — recorded as a failed metric, not used.

## Verdict for the Riemann Hypothesis Engine (a reading, with the numbers above)
- The worm's turning law is NOT a law of the ζ pathway: ζ's tangent turns one way, its wobble decays as ½p^{-½}.
- The ζ pathway and the Gemini curves share a LANGUAGE: sums of wheels with phase/amplitude modulation, locked into rosettes when frequency ratios are rational (never for ln p/ln q). The engine should hold a general wheel-sum object (frequencies ln p, radii p^{-σ}) and this viewer as its idealization sandbox; the Gemini worm/rosette are the 2–3-wheel rational special cases.
- (The u = ln p test below settled the open item that stood here.)

## The u = ln p clock test (`boundary.zeta_loops`, `boundary_uclock.svg`, 1 test added; 118 pass)
Gemini's worm turns by a CONSTANT angle per step of its parameter θ_b. Reading θ_b ≡ t·u (u = ln p), which is exactly the ζ carrier −t ln p, three comparisons:
1. **Clock — ALIGNED.** ζ loops (one 2π of the unwrapped tangent angle) are uniform in u: spacing 0.4444 against 2π/t = 0.4445 at t = γ₁ (rel. sd 2.4 %); 0.2095 against 0.2094 at t = 30 (4.6 %). In prime index n the same loops go from 7 to 56,535 steps apart (geometric, ratio e^{2π/t}): Gemini's clock fits as u and fails as index.
2. **Radius law — NOT aligned.** Radius growth per loop. ζ away from zeros is a log-spiral with the constant factor e^{π/t}·(u_j/u_{j+1}): measured 2.19 / 1.60 / 1.26 at t = 3 / 5 / 10 against the smooth law 2.29 / 1.64 / 1.28 (1–4 %, from |z|/R(p) = 1.04–1.05, R = p^{½}/(|½−it| ln p)). On or near a zero the growth is suppressed (1.03 at γ₁ against 1.19; |z|/R = 3.9) — the zeros cancel the smooth growth. The worm is hyperbolic: 1.25, 1.33, 1.49, 1.95, 19 (the inversion r → 1/r of a linear inward spiral), independent of τ. (The worm's first-loop factor 1.247 is within 0.2 % of e^{π/γ₁} = 1.249: one number, noted as a coincidence, not used.)
3. **Wobble — NOT aligned.** ζ's tangent wobble touches ½ p^{-½} as an envelope (panel c); the worm's is 5.
Result: of the worm's three ingredients — clock, radius law, wobble — only the clock (the carrier t·ln p) is shared with ζ. Gemini's constant-angle spiral is the ζ carrier in u; everything that makes it a worm is its own. For the Riemann Hypothesis Engine: keep the carrier clock u = ln p and the wheel-sum language; do not carry the worm's radius law or PM index as ζ mechanics.

