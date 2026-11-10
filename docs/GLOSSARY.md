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
