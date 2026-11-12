# RHP glossary — words that collide

Maths in Unicode. Each entry says what the word means **in RHP** and where the code is. Labels: ESTABLISHED = textbook; OURS = this repo's naming.

| word | meaning in RHP | where | collides with |
|---|---|---|---|
| **zero-tone** (was written bare "tone", "tones", "K tones") | one conjugate pair of zeros of ζ, as ONE real standing-wave component of the prime wave: −2 cos(γu − arg ρ)/|ρ|, u = ln x, frequency γ, amplitude 2/|ρ|. "K zero-tones" = the first K pairs (γ₁ … γ_K). A frequency component, in the sense of a Fourier term. ESTABLISHED (explicit formula); the name is OURS. | `rhp/wave.py` `tone(u, γ)`, `heartbeat(u, zeros)` (the sum of K zero-tones); `docs/heartbeat_waveform_derivation_2026-11-08.md` | **tonne** (1000 kg, mass; nothing here is a mass). Also the **tone** of the invisible-Unicode tone/tense channel on words (SemPriCon, "tone quasi-particles") and **tone** in the cymatics sonification engine: different objects, same word. |
| **heartbeat** | H(u) = Σ over zero-tones; the wave that rebuilds (ψ(x) − x)/√x | `rhp/heartbeat.py`, `rhp/wave.py` | the Gemini sim's "Heartbeat" θ = (2/π)tn (a spiral clock) |
| **wind** | |W_T(u)|, the size of the prime wave; its gusts and calms | `rhp/wave.py` `wind_stats` | meteorology |
| **ping** | per-integer test: the recovered jump ψ_K(n+½) − ψ_K(n−½) from K zero-tones, read against Λ(n); "prime power" if above ½ ln 2. A bash-style ping of n, not the Modulus ping. | `experiments/zeros25k_rerun.py` `jump_recovered` | the Modulus ping (VAPMIP) |
| **straggler** | a composite with NO prime factor among the primes that hold slots (the y-rough composites for y = the last slotted prime). OURS (Cody); literature: y-rough / y-sifted numbers, Buchstab function Φ(x, y) | `docs/sedenion_slot_rule_2026-11-10.md` | |
| **slot** | a basis index k of T_N (N = 2^m), e_k; L_k = left multiplication by e_k, L_k² = −I for k ≥ 1 | `rhp/sedenion_slots.py` | slot in Python/hardware |
| **carrier** | the clock t·ln p shared by the ζ path and a wheel (frequency ln p in u = ln p) | `rhp/boundary.py` | radio carrier (the same idea) |
| **wobble index** | the phase-modulation index of a step-angle wobble: ½p^(−½) for ζ, 5 for the Gemini worm | `rhp/boundary.py` `zeta_wobble` | |

Rule from here: write "zero-tone" (hyphenated) for the object above and never bare "tones" in RHP text. The earlier documents (log items 14, 20, 25, `heartbeat_waveform_derivation`) say "tones" and mean zero-tones.

## Symbols that collide (added 2026-11-12; Cody: 'the same symbol' and the three σ were mixing up the 'find the curve' context)
Source for σ: `Ainulindale/wiki/123_sigma_rb_mass_gap_calibration.md`. From now on RHP text writes the subscripted form; a bare Ω, α, d\* or σ is not used.
| symbol | the objects it has named | RHP writes |
|---|---|---|
| **σ** | (1) the real part of s: Re s, the Dirichlet/Euler coupling exponent, critical line σ = ½ (`h_rb_hat`, ζ); (2) the σ-strata σ₀…σ₄ of the sedenion, ordinal nested-inclusion levels (how many of the 16 components are on; σ₄ = 'string / ground state'); (3) σ_RB[k] = ψ[k]·conj ψ[k⊕4], the oblique-gear tilt (Re) and axis (Im), and the derived scalar σ_self | **σ_ζ**, **σ-stratum**, **σ_RB**. In this repo σ_ζ is the one used (`rhp/nodal*.py`, strip pictures); σ_RB and the strata are not used |
| **Ω** | (1) **Ω_ζΣ = W(1) = 0.567143…**, the Lambert fixed point of x = e^{−x}, the entropy-chase constant (First Age Ω paper, `rhp/omega_dstar.py`); (2) **ω_H = e^π = 23.1406926…** (Gelfond; the recovery note proposes the names OMEGA_W and OMEGA_GELFOND; chat 49156a9e, Ĥ_AΩ document), Gemini's 'Hagedorn ceiling' in the Berry–Keating scaling 2/ln ω_H = 2/π (`conjecture_proof.py`, `OMEGA`); (3) the 'Entropy (Ω)' slider of the Gemini widgets; (4) Ω(n), the number of prime factors with multiplicity (`rhp/nshape.py`) | **Ω_ζΣ**, **ω_H**, 'entropy slider', **Ω(n)** |
| **α** | (1) **α_fs** = 1/137.035999, the fine-structure constant, in this project the SOURCE of the rule that a quantity must say the same thing in every version of the maths (error-check discipline); (2) **Α_π / A_fermat** (capital Alpha, subscript π or fermat), the Fermat-side structural constant (the 'domain floor'); (3) **α_NN**, the neural fine-structure constant, defined to say the same thing at every tower layer; (4) the widget slider 'Inertia (α)' and the 'entropy coupling α' of the polar field | **α_fs**, **A_fermat**, **α_NN**. A_fermat is NOT α_fs and α_fs is not an ingredient of the horizon mechanism (Cody, 2026-11-10 and 2026-11-12). **Original meaning (recovered note, claude.ai chat 2a6155dc, 2026-04-09): 'ALPHA' = the d = 0 endpoint, T = T_P, the start of the forward curve (normalised ceiling x = 1)**; the 1/137.036 identification is later. A third object to confirm: the field tensor **Α^μν** of the Red-Blue Hamiltonian (chat 49156a9e, 04-11/12) |
| **d\*** | five values: d\*_taut = Ω_ζΣ/ln 10 = 0.24631 (exact, = −log₁₀ Ω_ζΣ); d\* 'Boundary' 0.24631; d\*_RG (dimension 8, numeric open); d\*_ln10 = Ω_ζΣ; and the measured d\*_spec = 0.24600 | always with its subscript; wiki 17, `reference_four_dstar_values` |
| **W** | Lambert W with its branches (W₀, W₋₁); W(0), W(1) = Ω_ζΣ, W(−1/e) = −1 (the fold) | `rhp/omega_dstar.py` |

### Which number is which (added 2026-11-12, after the question 'is Ω_ζΣ Omega_Riemann and d\* Alpha_fermat?')
| name | what it is | x = T/T_P | digit coordinate d | source |
|---|---|---|---|---|
| **Ω_ζΣ** (Omega_Riemann) | W(1), the crossing of the forward and inverse curves | 0.5671 | 0.2463 | wiki 17; `derive_omega_zs`; recovered 04-09 code |
| **d\*** (d\*_taut) | the same crossing in decimal digits, −log₁₀ Ω_ζΣ = Ω_ζΣ/ln 10; equivalently the DISTANCE in digits from the d = 0 endpoint to the crossing | not a separate point | 0.2463 | recovered code; docs/omega_dstar |
| **Α_π / A_fermat, original 04-09 code ('ALPHA')** | the d = 0 endpoint, T = T_P, start of the forward curve (the ceiling) | 1 | 0 | recovery note, chat 2a6155dc |
| **Α_π, later usage** | 1/137.035999, the 'domain floor' (E8/Wyler) | not on this scale | n/a | wiki 17; `derive_alpha_fermat`; First/Second Age |
| prime-spacing crossing (RHP) | x ln x = 1: mean gap ln x = 1/x; the inversion x → 1/x image of the Ω crossing | 1/Ω = 1.7632 | +0.2463 | docs/gap_side |
OPEN: wiki 17 puts the domain at [Α_π = 0.0073, Ω_ζΣ = 0.567] and says d\* sits inside it; the original derivation puts Alpha at x = 1. The two descriptions of the Fermat-side constant are not reconciled (Cody to say which is meant).

### The four d\* faces against RHP's measurements (2026-11-12; docs/four_dstar_2026-11-12.md)
d\* (Boundary) = the digit length D\* of the crossing; d\*_taut (Flow) = digits of the fixed point gained per iteration from the ceiling ('one step of r'); d\*_ln(10) (Translator) = the mean prime gap ln X\* at the crossing = d\*·ln 10 = Ω_ζΣ, and mean gap ÷ digit length = ln 10 at every X (PNT); d\*_RG (Stability, dimension 8) not touched. READING: one crossing, two coordinates: digits (Riemann / entropy side) d\* = 0.2463 and gap count (Fermat / inertia side) Ω = W(1) = 0.5671, related by ln 10.

