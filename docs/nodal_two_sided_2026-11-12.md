# Nodal lines with both sides mapped — 2026-11-12   :CALCULATED

Code `rhp/two_sided.py`; run `experiments/nodal_two_sided.py` → `.json` (reuses the earlier `xi_grid.npz`); test `test_nodal_two_sided`.
Status before this run: item 31 extracted the nodes of the prime wave (ψ(x) = x) and the ξ strip ladder, BEFORE the Möbius side existed. This run puts the two sides together. (Cody, 2026-11-12: the capital Alpha_fermat is the structural constant; the lowercase fine-structure α is only the source of the 'says the same thing in every maths' error-check rule and is not an ingredient here.)

## R3 — 1-D, exact data to 10⁷: the ψ wave and the M wave
Nodes = sign changes of ψ(n) − n and of M(n). Locking = the normalised position of a node of one wave inside the bracketing gap of the other (0.25 = uniformly random).
| x | ψ nodes | M nodes | ψ nodes per M node | M node inside ψ gap | ψ node inside M gap | opposite-sign fraction |
|---|---|---|---|---|---|---|
| 10³–10⁴ | 471 | 111 | 4.2 | 0.211 | 0.214 | 0.780 |
| 10⁴–10⁵ | 1,500 | 348 | 4.3 | 0.225 | 0.208 | 0.800 |
| 10⁵–10⁶ | 4,567 | 1,161 | 3.9 | 0.234 | 0.205 | 0.812 |
| 10⁶–10⁷ | 12,224 | 2,214 | 5.5 | 0.224 | 0.211 | 0.782 |
(null for the two locking columns 0.249–0.252.) The two waves are in antiphase about 80 % of the time (corr −0.82, docs/mertens_wave), the ψ wave has 4–5 nodes for each M node, and the nodes are **weakly but significantly locked**: each M node sits 10–15 % closer to a ψ node than chance (≈ 8σ for the M-in-ψ column, ≈ 30σ for ψ-in-M). It is a resonance, not a lock.

## R4 — 2-D: ξ and 1/ξ on the strip (σ ∈ [0, 1], t ∈ [0.5, 80], 3,181 × 201 grid)
- **Same nodal lines.** Re(1/ξ) = Re ξ/|ξ|² and Im(1/ξ) = −Im ξ/|ξ|²: the signs agree at 100 % of the 639,381 grid points (both), and the extracted Re = 0 line sets match (21 lines each, total length 50.70 vs 50.86, the 0.3 % being extraction noise near the poles). The T side and the 1/ζ side share the nodal structure exactly; what differs is the weight at the crossings: ξ has a zero there, 1/ξ a pole of residue 1/ξ′(ρ) (the M-tone amplitude 1/|ζ′(ρ)|).
- **Right angle at every zero.** ξ is real on the line, so near a zero ξ ≈ ξ_t (δt − iδσ): the Re ξ = 0 rung is horizontal (δt = 0) and the critical line is Im ξ = 0 (δσ = 0). Measured angle of each rung from the horizontal at the equator, all 21 zeros below t = 80 (42 vertex hits, each counted twice): median |angle| 0.38°, maximum 0.80°. The nodes of the 2-D standing wave are perpendicular crossings of a rung with the critical line.

## Reading (not a result)
The standing wave in the arena has two compatible descriptions with the SAME nodal skeleton: rungs ⟂ the equator at the zeros (2-D), and, along the prime axis, an antiphase pair W_ψ, W_M whose nodes interleave with a weak resonance. Not done: the Fermat-side finite-tone wave (docs/fermat_horizon) has not been put through the same extraction, and the prime-spacing curve (docs/gap_side) is not a wave.
