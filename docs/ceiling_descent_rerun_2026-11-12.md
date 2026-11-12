# The ceiling-descent experiment, rebuilt from claude.ai's primer — 2026-11-12   :CALCULATED

Source: the primer claude.ai wrote on request (private Artifact `https://claude.ai/artifact/7tmDnxd2m2ecb2wUfvn4hA`, saved verbatim as `ContextPlease/claude/hist_prime/RiemannHypothesisProof/PRIMER_2026-11-12_CEILING_DESCENT_FROM_CLAUDE_AI.md`). Code `rhp/polar_field.py`, `rhp/omega_dstar.py`; tests `test_polar_field`, `test_omega_dstar`.

## What the primer establishes (claude.ai's account of its own history)
- The **entropy chase** (April 9) is the digit-curve crossing; the step S/k_B = T/T_P ⇒ T/T_P = e^{−T/T_P} is **assumed**: the Bekenstein form is linear in T and the substitution only fixes the radius r = l_P/2π for every T (checked here: 2.57235·10⁻³⁶ m, exact).
- **Inertia, version A** (April 9, RG form): g e^{g/g0} = g0 ⇒ g = g0·Ω is the SAME equation u e^u = 1 as the entropy chase (u = g/g0); no β-function on record.
- **Inertia, version B** (the Gemini 'α line' from the fine-structure floor to the speed of light, whose intersection with Ω's line gave d* ≈ 0.246) is **not in any Claude record**; it was a Gemini session, never supplied. Unrecovered. Not reconstructed here.
- 'Tangency, K = 0 at d*' was stated in two papers and never computed. 'Fermat prime generator': no record under any name.

## Rerun here (independent port; agrees with the primer's RERUN numbers)
| | primer | here |
|---|---|---|
| min \|S − I\|, widget defaults (α 0.137, depth 1.57, ln) | 0.1471267 at r 0.2563 | 0.1471267 at r 0.25633 |
| log₁₀ | 0.1645164 at r 0.6268 | 0.1645164 at r 0.62685 |
| α = 1/137.036, depth π/2 | 0.1750564 at r 0.2346 | 0.1750564 at r 0.23460 |
| base terms S₀ − I₀: crossings on (0, 50]; min \|S₀ − I₀\| | none; 0.3132 | none; 0.3132 |
| 1-D forms at d = 0.246: S, I | 0.0016037, 0.3622851 | 0.0016037, 0.3622851 |
| K of z = S − I, ring mean at r = 0.1, 0.2, 0.246, 0.3, 0.5, 1.0 | −4.0316, +0.2124, +0.8184, +0.5972, −0.4712, +0.0045 | identical to 4 digits |
S and I never meet; the 'd* ring' is the radius of closest approach, not an intersection or a tangency; log₁₀ moves it to 0.627, not Ω = 0.567.
**New (the record left the surface open):** K of z = S: ring means −7.70, −1.26, −0.79, −0.49, −0.15, −0.018 at the same radii; of z = I: −15.5, −2.47, −1.37, −0.79, −0.21, −0.038. Both are negative at every radius tested; only S − I changes sign (mean > 0 near r = 0.2–0.3), and none is zero at r = 0.246 on any surface. 'K = 0 at d*' is not supported by any of the three.

## The ceiling scale (resolves a discrepancy)
- Cody (2026-11-12): the information ceiling he chose is the Hagedorn ceiling, far below the Planck-wavelength temperature; beyond it signal and noise are indistinguishable.
- The code states it as **140 quadrillion °F** (`conjecture_proof.py` 140e15; `FirstDraft/FirstAge/Gemini/event_horizon_sim.py`: `OMEGA_K = 77.8e15  # Noether limit; from 140Q °F × 5/9`) = **7.78·10¹⁶ K**. The Second Age text's '≈ 1.4·10¹⁷ K' and wiki 17's 'Kelvin' drop the conversion: a unit slip, now explained.
- The April-9 run on record used **T_Planck = 1.4168·10³² K** (primer §1.2a). The dimensionless results do not depend on the scale: Ω = 0.5671432904…, d* = |log₁₀ Ω| = 0.2463072014734. Only T* does: Ω·T_Planck = 8.0352·10³¹ K; Ω·(7.78·10¹⁶ K) = **4.411·10¹⁶ K**. Which scale the original run used for T* is Cody's to confirm.
- 'ω_H = e^π' (the later dimensionless Hagedorn notation) was not an input to the April runs (primer).

## Status of the Fermat-side chase
Version B cannot be built from what exists. To recover it: the Gemini conversation of ~2026-04-09/10 (Google Takeout → 'My Activity' → Gemini Apps), or ask Gemini for the same kind of primer (forward curve, inverse curve, both ceilings with units, the crossing equation, the numbers at every step, the code). Until then the Fermat side stands on the exact finite-field objects (docs/fermat_horizon_2026-11-12.md) and not on a derivation of A_π.
