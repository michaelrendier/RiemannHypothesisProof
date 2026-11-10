# The nodal-line extractor   2026-11-10   :CALCULATED (no claim about RH)

Code `rhp/nodal_extract.py`, run `experiments/nodal_extract_run.py` (→ `experiments/nodal_extract/{r1,r2}.json`), figure `rhp/atlas_riemann_fermat/nodal_extract.svg`, test `test_nodal_extractor_core_and_results`. Cymatics: a Chladni figure is the zero set of a field made of standing waves; sand collects on it. The extractor takes a sampled field and returns its zero set, the crossings of two zero sets, and the sign domains.

## Extractor (no plotting library)
`zero_set_2d` marching squares with linear edge interpolation and the cell-mean rule at saddles, segments stitched through shared edge crossings into polylines (open / closed, length); `seg_intersections` (exact, deduplicated); `sign_domains` (Courant's nodal domains, 4-connected); `sign_nodes` / `refine_nodes` (1-D nodes with direction); `exact_prime_nodes` (below); `zero_tone_wave` (the wave rebuilt from K zero-tones, see GLOSSARY); `rice_rate`; `sand`.
Checks: unit circle → one closed polyline of length 6.28311 (2π to 1.2·10⁻⁴); Chladni mode cos3x cos2y − cos2x cos3y total nodal length 12.5476 = matplotlib's; z² − 1 → zeros at (±1, 0).

## R1 — the prime wave W(u) = (ψ(x) − x)/√x: nodes are the prime crossings ψ(x) = x
Exact nodes: ψ is constant between integers and jumps up by Λ(n) at n. On (n, n+1) the function ψ(n) − x decreases and has a node at x = ψ(n) iff n < ψ(n) < n+1 (a **slide**); at x = n the jump is a node iff ψ(n−1) − n and ψ(n) − n differ in sign (a **jump**). Nodes therefore come in jump + slide pairs less than one unit apart (20 jumps, 19 slides in [1000, 1500]).
Exact nodes: 39 in x ∈ [1000, 1500], 19 in [3000, 3500], 10 in [10000, 10500], 0 in [30000, 30500] and [100000, 100200] (the wave stays on one side over a window that short).
Nodes of the zero-tone wave against the exact ones (a match = within 0.25 in x; chance level 0.02–0.10):
| window | K = 100 | 1,000 | 3,066 | 25,000 (exact→found / found→exact) |
|---|---|---|---|---|
| [1000, 1500] | 0.05 | 0.23 | 0.41 | **1.00 / 0.91** (47 found vs 39 exact) |
| [3000, 3500] | 0 | 0.05 | 0.16 | **1.00 / 0.78** (27 vs 19) |
| [10000, 10500] | 0 | 0 | 0.10 | 0.50 / 0.62 (8 vs 10) |
The extractor finds every exact prime crossing out to x ≈ 3,500 with 25,000 zero-tones; the match falls to 0.5 at 10⁴ — the same edge as the ping (perfect to ~10⁴, log item 28). The node count saturates near the exact count as K grows (47 against 39) while **Rice's law for a Gaussian band-limited wave does not**: it fits for K ≤ 100 (3.5 / 9.1 against 4 / 9 found) and overpredicts by 3× at K = 25,000 (134.5 against 47), because every added zero-tone adds derivative amplitude (a_kγ_k = 2) without adding amplitude (a_k = 2/|ρ_k|), so Rice's rate grows as √K while the true node set converges. Rice is a bound on chatter, not the node count of the prime wave.

## R2 — ξ on the strip σ ∈ [−0.5, 1.5], t ∈ [0.5, 80] (grid 3181 × 201, mirror symmetry ξ(σ+it) = conj ξ(1−σ+it))
- Re ξ = 0: **21 open lines, one per zero (21 zeros below 80), every one running across the whole strip from σ = −0.5 to 1.5, none closed**; each crosses the equator once, at a zero (heights agree with the table to 1.8·10⁻⁴ max, the grid step 0.025).
- Sign domains of Re ξ: 11 positive + 11 negative = 22 = 21 rungs + 1 bands: the rungs partition the strip into alternating bands.
- Im ξ = 0 off the equator (σ > ½; the equator itself is Im ξ ≡ 0): **20 lines, none closed, each leaving the equator at an extremum of the heartbeat** (F′(t) = 0, 20 of 20; median offset 0.0013 in t at σ = 0.51) and running to σ = 1.5, tilting by −0.15 … −0.34 in t (mean −0.23) over the half-strip. Zeros (rungs) and extrema (Im lines) alternate along the equator: 21 + 20.
- **Off-line crossings of Re ξ = 0 and Im ξ = 0: none** (in this region). That is the extractor's check against the established fact; it is not evidence about RH beyond what is already computed.
Structure in words: the standing wave's nodal set on the strip is a ladder: red rungs through the zeros, blue lines through the extrema, alternating, all running to the edge of the strip; nothing closes into a loop.
