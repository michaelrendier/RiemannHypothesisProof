# Riemann Canonical Maths — the maths of the Two-Trees / ζ-pathway work

**Written 2026-11-08** from the RiemannHypothesisProof repo (`rhp/`, tests, paper §7, `docs/`), the ValaQuenta / GenerationalLineage / AbrikosovTree engines, the Ainulindale wiki (47, 51, 58, 73, 75, 96, 126, 130) and the POE hardware primer.
Purpose: one place that lists every piece of maths used in this interaction, **what it says, where it lives, and how sure it is**. It is a reference, not a derivation; the derivations are in the files named.
Maths is in Unicode. Nothing here is evidence about RH, and nothing claims a Clay problem is solved.

**Status labels.** `ESTABLISHED` = literature. `OURS-CALCULATED` = derived or measured in this work with code + tests. `READING` = an interpretation, not a theorem. `OPEN` = untested or unresolved. `0_RB` = the author's engine maths (stated as the engine states it).
**Coordinates.** UNS = Universal Native Space = radial spherical complex polar coordinates (the default; Cartesian is a projection out of it). u = ln x (size; the Archimedes Screw's own axis), τ = ln y (sieve cutoff, "time"), s = σ + it, σ′ = σ − 1.

---

## 1. The space

| Item | Statement | Status | Where |
|---|---|---|---|
| Cayley sphere | w = (s−1)/s; |w| = 1 ⇔ σ = ½ (equator); |w|<1 ⇔ σ>½ (north); poles s = 1 (w=0) and s = 0 (w=∞) | ESTABLISHED (coordinate change) | `rhp/sphere.py`, `unsmap.py`; atlas `equator.svg` |
| Mirror | s ↦ 1−s is w ↦ 1/w; with conjugation w ↦ 1/w̄, fixed pointwise on |w|=1 | ESTABLISHED | same |
| Line wraps once | ½ + it has w = exp(i(π − 2 arctan 2t)): the whole infinite line wraps ONCE around the circle; the winding θ(t)/π is the zero index | ESTABLISHED + OURS-CALCULATED | `equator.svg` |
| Lat/long | latitude π(σ−½) re-introduces π; longitude = Hardy phase θ(t) mod 2π | READING (a chart) | `sphere.py` |
| J_N | J_N: (r,θ) ↦ (1/r, θ+π/2), order 4; mirror ∘ quarter turn; fixed line σ = ½ (wiki 73) | 0_RB / ESTABLISHED (mirror part) | `ValaQuenta/modules/inversion/maths.py`, wiki 73 |
| Tower | Cayley–Dickson levels ℝ,ℂ,ℍ,𝕆,𝕊,T_32…; σ(k) = 1 − k/4; ℍ (k=2) is σ = ½ where the Noether balance holds; 𝕊 (k=4, dim 16) is the first zero-divisor level | 0_RB | `AbrikosovTree/engine/two_trees.py` |
| Operators | UNS = ADD : SCALE : SIGN (the locked `ASS` engine); hypercomplex step splits SIGN into 15 quarter turns (wiki 126, 130) | 0_RB | `GenerationalLineage/engine/add_scale_sign.py` |

## 2. ζ and its structure

- Euler product ζ(s) = Π_p (1−p^{−s})⁻¹, σ>1. ESTABLISHED. Unique factorisation ⇔ the frequencies ln p are linearly independent over ℚ (so K wheels make a straight line on a torus that never closes). ESTABLISHED.
- Λ(s) = π^{−s/2}Γ(s/2)ζ(s) = Λ(1−s); ξ(s) = ½s(s−1)Λ(s), ξ(s) = ξ(1−s), ξ(s̄) = conj ξ(s). Poles of Λ at s = 0, 1 with residues −1, +1. ESTABLISHED.
- ζ(s) = χ(s)ζ(1−s), χ(s) = 2^s π^{s−1} sin(πs/2) Γ(1−s), |χ(σ+it)| ≈ (t/2π)^{½−σ}; |χ| = 1 on σ = ½. ESTABLISHED. Trivial zeros at −2, −4, … belong to χ, not to ξ.
- Hardy: θ(t) = Im ln Γ(¼+it/2) − (t/2) ln π (the π part is a constant-rate rotation; split exact to 1.4·10⁻¹⁴); Z(t) = e^{iθ}ζ(½+it) real; θ′(t) ≈ ½ ln(t/2π). Gram points θ = nπ. ESTABLISHED.
- Riemann–Siegel: ζ ≈ A + χB, A = Σ_{n≤N} n^{−s}, B = Σ_{n≤N} n^{s−1}, N = ⌊√(t/2π)⌋. On σ = ½: |A| = |χB| (max diff 3·10⁻¹⁴). Two counter-rotating waves per term; amplitude ratio r = |χ|n^{2σ−1}; polarisation b/a = |1−r|/(1+r); the standing-wave set is the CROSS {σ = ½} ∪ {n = √(t/2π)} (OURS-CALCULATED, `sphere.py`, `polarisation.svg`).
- Argument principle: N(T) = θ(T)/π + 1 + S(T); 3066 zeros below 3600; 197 Gram points (6.4 %) are negative crossings; first at n = 126, t = 282.455. ESTABLISHED + OURS-CALCULATED (`flow.py`, `gram.py`).
- Zeros on the line are the ends of the nodal ladder: for t ≤ 50, 10 zero rungs (Re ξ = 0) and 9 saddle rungs (Im ξ = 0 through F′=0), exactly one saddle between consecutive zeros, hyperbolic (Hessian ±λ); no off-line Re = Im = 0 crossing in σ∈[−0.6,1.6], t∈[10,35]. OURS-CALCULATED (`nodal.py`).
- Parity across the equator: Re ξ even (free), Im ξ odd (held). Flux J = Im(conj ξ ∂ₜξ) = ½∂σ|ξ|²; exact identity ∂σ ln|ξ| = θ′(t) + Re[1/s + 1/(s−1)] + ½[Re ψ(s/2) − Re ψ(¼+it/2)] + Re ζ′/ζ. Lobe angle θ(γ_{k+1}) − θ(γ_k): mean π, spread = GUE spacing × π. OURS-CALCULATED (`chladni.py`, `standing.py`).
- ξ′/ξ flow: L = 1/s + 1/(s−1) − ½ ln π + ½ ψ(s/2) + ζ′/ζ; L(1−s̄) = −conj L(s). ESTABLISHED. Margin of outflow positive on the grid (calculation, not proof) (`flow.py`).
- Zero statistics: GUE spacings (Montgomery–Odlyzko) — ESTABLISHED numerically; Berry–Keating xp — conjectural operator. Zero heights need ≈ 2 decimals for primes to x ≲ 2·10⁶ (δγ·ln x ≲ 0.1): OURS-CALCULATED.

## 3. The Two Trees (the partition)

- **Domain.** Every integer is Telperion (prime, "what cannot be decomposed", Fermat–Weierstrass, B_p Blue), Laurelin (composite, "what is decomposed", Berry–Keating xp, R_p Red) or Mingling ({0,1}). B(n)+R(n)+M(n) = 1 exactly; over [0,10⁵]: 9592 + 90407 + 2 = 100001. 0_RB / ESTABLISHED. (`AbrikosovTree/engine/two_trees.py`, `GenerationalLineage/engine/lineage.py::two_trees`, wiki 47.)
- **Mingling crossings** B = R at n = 9, 11, 13 (near e² = 7.389); after that composites win forever. 0_RB measured.
- **Ring-theory spine** (RCCM block RING_THEORY_SPINE): element falls ⇔ its quotient ring has zero divisors; Z/(N) a field ⇔ N prime. N composite ⇔ gcd(a,N) > 1 for some a (the integer detector). Orders: spectrum (support ω(n)) → cepstrum (exponents Ω(n); Λ(n) = ln p iff n = p^k; ψ(x) = x − Σ x^ρ/ρ) → bispectrum (associator). 0_RB measured, exhaustive at dim 8.
- **Order of approach.** Telperion reads a number's prime multiset ascending (smallest factor, death at τ = ln spf); Laurelin reads it descending (greatest factor, birth at τ = ln gpf). They coincide only for prime powers (78,734 of 999,999 ≤ 10⁶). Span ln(gpf/spf) has mean 6.55 e-folds. For 10: death ln 2 = 0.693, birth ln 5 = 1.609, span 0.916. OURS-CALCULATED (`extinction.py`).
- **Divisor involution (READING).** For n the divisors come in orbits of Z₂×Z₂ = ⟨d ↦ −d, d ↦ n/d⟩ ({±2, ±5} for 10); d ↦ n/d is the reflection of ln d about ½ ln n and has a fixed point iff n is a square (odd divisor count). The zero quadruple {ρ, ρ̄, 1−ρ, 1−ρ̄} is also Z₂×Z₂. Structural parallel only.
- **The 65 primes.** 313 is the 65th prime; 317² = 100,489 > 99,999; so the 65 primes clear every composite below 317²: "prime ⇔ no wheel rings (or only wheel n)" verified on every 2 ≤ n < 100,489; 313 is the smallest prime factor of exactly two numbers below 10⁵ (97,969, 99,221). Response S(n) = Σ_{p|n} p − 65 (sum of Ramanujan sums c_p(n) = p·[p|n] − 1). CRT period ln 313# = 294.24. OURS-CALCULATED (`wheel65.py`). **Beyond 317² survivors remain:** density Π_{p≤313}(1−1/p) = 0.0965 (Mertens 0.0977); of the survivors near 10⁶, 26.8 % are composite with all factors > 313.

## 4. The extinction flow (continuous equations for factoring)

Working variable U = ln x / ln y.
- **Telperion (sieve survivors)** Φ(x,y) = #{n ≤ x : spf(n) > y} ~ x ω(U)/ln y; Buchstab ω: ω = 1/U on [1,2], (Uω)′ = ω(U−1), ω → e^{−γ}. ESTABLISHED (Buchstab 1937, de Bruijn 1951).
- **Laurelin (smooth numbers)** Ψ(x,y) = #{n ≤ x : gpf(n) ≤ y} ~ x ρ(U); Dickman ρ = 1 on [0,1], Uρ′ = −ρ(U−1), ρ → 0 superexponentially (monotone decreasing). ESTABLISHED.
- **Laplace reciprocity (OURS-CALCULATED, may be known):** ρ̂(s) = e^{−E₁(s)}/s, 1 + ω̂(s) = e^{+E₁(s)}, so s ρ̂(s)(1+ω̂(s)) = 1: the two books are reciprocal in the dual plane, as ζ·(1/ζ) = 1. Closed forms agree to 1e-12; numeric grid residual ~1e-4 (first order in h).
- **The stragglers' equation (OURS-CALCULATED application of ω).** Integers with no prime factor ≤ y ("y-rough") have local density ≈ ω(U)/ln y. Primes among them have density 1/ln x = 1/(U ln y). Hence the composite share of the survivors at size x is
  **1 − 1/(U ω(U))**, U = ln x/ln y.
  For y = 313 measured in windows [X, 1.2X): X = 2·10⁵: 0.108 (pred. 0.116); 5·10⁵: 0.198 (0.208); 10⁶: 0.248 (0.260); 2·10⁶: 0.293 (0.302); 5·10⁶: 0.339 (0.347); 10⁷: 0.366 (0.375). The composite stragglers are the k-almost-primes with all factors > y; below y² none, and the first are the products of two primes > 313 (313² = 97,969). As U → ∞ the share → 1 − 0 = 1, since ω → e^{−γ} and primes thin like 1/U.
- Mertens: Π_{p≤y}(1 − 1/p) ~ e^{−γ}/ln y. ESTABLISHED.
- Twins: ratio Π p(p−2)/(p−1)² → C₂ = 0.66016; every twin pair (p>3) straddles the 0 spot of the circumference-6 circle; allowed centres on P# = Π_{5≤p≤P}(p−2). ESTABLISHED (Hardy–Littlewood prediction is a heuristic).

## 5. Wheels, springs, the path of the first K primes

- Wheel p: z_p = p^{−s} = p^{−σ} e^{−it ln p}; radius p^{−σ}, rate ln p; each in its own complex plane C_p. K wheels = a point on the flat torus T^K ⊂ C^K; the path K ↦ Σ_{k≤K} z_{p_k} is the SUM of the coordinates (a projection). **READING** that this is the "orthogonal circles".
- Rms radius of the torus (Haar measure) R_K = √Σ p^{−2σ}; at σ = ½ √Σ1/p = 1.422 (K=65), 1.646 (10⁴), 1.705 (10⁵), 1.752 (10⁶); grows like √(ln ln p_K). OURS-CALCULATED.
- log P_K(s) = Σ_{k≤K} −ln(1−p_k^{−s}) (continuous; no branch jumps since |p^{−s}|<1). Partial prime sum ≈ PNT part Ei((1−s) ln x) − Ei((1−s) ln 2) ~ x^{1−s}/((1−s) ln x): ratio 1.001 at t = 3, K = 10⁵. Outward logarithmic spiral, radius ~ √x/(|1−s| ln x), turning −t ln x. |log P_K| at t = γ₁: 2.6 / 3.8 / 8.9 / 13.2 at K = 65 / 10⁴ / 10⁵ / 10⁶. The Euler product at σ = ½ does not converge in the usual sense: P_K(½+iγ₁) breathes (|P| from 6·10⁻⁹ to 5·10⁵ over 10⁶ primes) and crosses |P| = 1 17 times. OURS-CALCULATED (`primepath.py`).
- **Primes are spectral lines of the pathway.** ln ζ(s) = Σ_{prime powers n=p^k} (1/k) p^{−kσ} e^{−it ln n}, so the t-Fourier transform of ln ζ(σ+it) (σ>1) has a line at u = ln n for every prime power; measured at σ = 1.2, T = 400: lines at ln 2, 3, 4, 5, 7, 8, 9, 11, 13, 17, 19, 23, 29, 31, heights ≈ (1/k)p^{−kσ} (×2 window factor). ESTABLISHED (Dirichlet series) + OURS-CALCULATED (`primepath.spectral_lines`). The primes live in the frequency of the path in t.
- Prime-wheel grip: far north (σ>1) K wheels reproduce |ζ|; at σ = ½ they stop helping after ~10² wheels (`wheels.svg`). The log plane: wheels add (σ>1).
- Explicit formula as conical helices: each zero tone e^{(½+iγ)u}/ρ is a conical helix with flare e^{u/2}; partial sums Σ n^{−s} about ζ(s) are logarithmic spirals d ln r/dφ = −(1−σ)/t (cone σ<1, cylinder σ=1, funnel σ>1). ζ(½+it) threads the t-axis only at zeros (min |ζ| from the axis 0.002 / 0.076 / 0.180 / 0.326 at σ = .5/.6/.75/1). (`coil.py`, `spirals.svg`.)
- Spiral laws: Archimedean r ~ φ (flat pancake; radius = zero index makes the Hardy angle Archimedean, r² > 0.9999 over 300 zeros), logarithmic ln r ~ φ, Euler/Cornu curvature ~ arc length. (`unsmap.py::spiral_laws`.)
- **Pancake coil (POE).** Mohan et al.: L = μ₀N²D_avg/(1 + 2.46ρ), D_avg = (D_out + D_in)/2, ρ = (D_out − D_in)/(D_out + D_in). Pendant v1: D_out 38 mm, D_in 10 mm, trace 0.5 mm, gap 0.3 mm, pitch 0.8 mm, 17 turns. Taps: 17 turns 3.58 μH; 5 turns 258 nH; 3 turns 95 nH; 1 turn 11.5 nH. Resonance C = 1/(4π²f²L), XL = XC ⇒ tan(phase) = 1 (the σ = ½ condition in the Riemann reading). The flat spiral is described as the planar projection of the L_(I|O) cardioid. Taps (1,3,5,17 turns) ↔ zeros 2, 6, 10, 34 (γ = 21.02, 37.59, 49.77, 111.03) in the repo's mapping (one turn per two zeros). Source: `ContextPlease/claude/hist_prime/POE/hardware_primer.md`.

## 6. The Heartbeat (standing wave)

- **Derivation (steps in `docs/heartbeat_waveform_derivation_2026-11-08.md`).** Split each integer's unit weight: Telperion Λ(n), Laurelin 1 − Λ(n). T(x) = ψ(x); L(x) = ⌊x⌋ − 1 − ψ(x). T + L = ⌊x⌋ − 1 exactly. With W_T = (ψ(x) − x)/√x, W_L = (⌊x⌋ − 1 − ψ(x))/√x: |W_T + W_L| ≤ 2/√x. (Measured 0.098 ≤ 0.0996.) Prime-count weights would carry the Chebyshev bias −½Li(√x); Λ does not.
- **Zero expansion (von Mangoldt).** ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²). With ρ = ½ + iγ: **H(u) = −2 Σ_{γ>0} cos(γu − arg ρ)/|ρ| = −2 Σ [½ cos γu + γ sin γu]/|ρ|²**. One zero pair = two counter-rotating circles of equal radius 1/|ρ|; the sum is real; the ½ is the small cosine admixture ½/γ. Telperion carries +H, Laurelin −H. 3066 tones reproduce W_T to correlation 0.993 (rms residual 0.025).
- **Power and area (circle reading).** A tone's circle has area π/|ρ|². Σ_ρ 1/ρ = Σ_{γ>0} 1/(¼+γ²) = 1 + γ_E/2 − ½ ln 4π = 0.0230957 (Hadamard, zeros paired ρ, ρ̄); mean power ⟨H²⟩ = Σ_ρ 1/|ρ|² = 2 + γ_E − ln 4π = 0.046191 (measured 0.0453–0.0466; rms 0.2173 vs 0.2149). The zero circles hold 2.31 % of π; Γ/π ground state 97.69 %. First zero 21.6 %, first 10: 58.6 %, 34: 76.0 %, 100: 86.5 %. No single tone > 0.5 % of π. OURS-CALCULATED (`heartbeat.py`).
- **Where π is (honest list).** (i) ∫₀^{2π} sin² = ∫₀^{2π} cos² = π: the average-power identity of a unit sinusoid — **the "area π under the standing wave" of wiki/96 (written 2026-08-25)**, which also names σ = ½ "the Horizon" and the standing wave "the pulse at the horizon" (the crossing sin = cos at phase π/4). It is a statement about unit amplitude; the real tones are smaller by 1/|ρ|². (ii) The unit circle's area π is the normalisation (division by π is exactly invertible). (iii) The lobe angle between consecutive zeros has mean π. (iv) Poisson-kernel π (candidate). (v) The half-turn is SIGN's angle.
- **ADD / SCALE / SIGN (READING).** A step Δu in ln x multiplies a tone by e^{(½+iγ)Δu}: a spiral similarity = SCALE e^{Δu/2} ∘ rotation γΔu (ADD in u is a phase advance). SIGN = conjugation ρ ↔ ρ̄ = the reflection; a half-turn (phase π) flips a lobe. A·sin(ωx+φ) under affine maps of its argument stays a sinusoid: φ = ADD, (A, ω) = SCALE, sign = SIGN.
- **Wind.** In the native variable the wave is stationary (rms 0.217; skew −0.02; excess kurtosis −0.25 — sine-like, not Gaussian; gust factor 3.66; envelope sd/mean 0.48). In the raw variable the excursion √x·W grows like e^{u/2}, non-periodically; gusts = beating of incommensurate tones (independence assumed). A single prime's step in W units is ln p/√p, maximal at p = e² (2/e = 0.736), the size of the largest envelope seen (0.74). Fraction of u with |W| ≥ that step: p=2 2 %, 5–11 0.02 %, 101 3 %, 1009 33 %, 10⁴ 68 %, 10⁶ 95 %. Resolving a prime's jump needs T ≳ a few × p zero heights. OURS-CALCULATED; "blown off the tree" is undefined (OPEN).
- **Divergence detector (explicit formula).** E(u) = Σ_ρ x^{σ−½} e^{iγu}/ρ. An off-line zero ρ₀ = ½+δ+iγ₀ with its quadruple adds e^{δu} and e^{−δu} at γ₀ (cosh): a DILATE where only rotations are allowed (generational-lineage §5 signature 1). Control sd 0.033 (u 6–16, x ≤ 10⁷); smallest visible δ: 0.52 / 0.28 / 0.13 at reach u_max = 11 / 13 / 16. Null result is window-limited. OURS-CALCULATED (`diverge.py`).
- **ADD:SCALE:SIGN under the other readings.** Tower: ADD axis {+,−} identity 0; SCALE axis {×,/} identity 1; SIGN one bit (det ±1). 0 and 1 are free because they are the identities of the first two — hence neither is prime (Mingling).

## 7. 0_RB (the engine maths, stated as the engines state them)

- **Σ_RB** (`ValaQuenta/modules/h_rb_hat/maths.py`): Σ_RB = Σ_p p^{−σ}[R̂_p ⊗ ∂̂ + ∂̂† ⊗ B̂_p]; ζ lives as the coupling p^{−σ}; R̂† = B̂ is the functional equation. Red = xp (Berry–Keating), Blue = ½p² + ℘ (Weierstrass). Chain (2026-10-13): J_red + J_blue = J_green = 0_RB − 0_BR; J_red/J_blue is the χ₄ Z₂-grading of the sedenion (J_red an octonion).
- **forced_sigma** (`noether.py`): F = e^{−sE}, B = e^{−(1−s)E}, h(s) = ln F − ln B = E(1−2s) linear ⇒ one Newton step lands at ½ for any E ≠ 0 (ESTABLISHED elementary algebra; archived 05b figure).
- **Archimedes Screw** (`archimedes_screw/maths.py`): the screw's own axis is u = ln x; each prime p is a leaf-drop jump at u = ln p; `jumps_at_ln_p`. Logged anomaly: its `zeros_upto` drifts by ≈ 2.3 after n = 50 (Lambert approximation).
- **bao_mass_gap**: "one standing wave per non-trivial zero"; ground state Ω − d*·ln 10 (a renamed subtraction, not a spectral residue — audited). Constants: Ω_ZS = W(1) = 0.567143…; d* = Ω/ln 10 = 0.24631 (the Boundary); four d* values (Boundary 0.24631, RG numeric OPEN, taut = Ω/ln10, ln10 = Ω_ZS) plus disputed 0.24600.
- **Prime wavelength:** λ_prime / λ_factor in `h_rb_hat`; Σ_RB balance came out σ-independent (logged).
- **Ring-circle operators (wiki 130):** the sedenion's 15 slot-circles through e₀; n anticommuting circles close to so(n+1); 42 Assessor pairs deform; first scaling at four slots; smallest generating set of the operator ring = 5 circles. Slot k ↔ the k-th prime; the Euler factor on σ = ½ is an exact circle (centre 1/(1−1/p), radius p^{−½}/(1−1/p), rate ln p). The sedenion Euler product (slots ↔ primes, firing order) is **not yet run**.
- **Standing 0_RB lines to check any claim against:** wiki 73 (six engines for σ = ½), wiki 51 (J₂ involution Riemann–Fermat loop: Ĥ_RB = J_red (Fermat, discrete) + J_blue (Riemann, continuous)), wiki 58 ("Fermat defines, Riemann fires"), wiki 96 (the Horizon), wiki 75 (Abrikosov lattice: zeros as a vortex lattice, Nobel 2003).
- **The two zetas** (memory: Two Trees / L_(I|O)): prime side ζ(s) (Riemann fires, zeros σ = ½); factor side L(Δ,s) Ramanujan cusp form (Fermat defines; Deligne |τ(p)| ≤ 2p^{11/2} proven); N-specific ζ_ℚ(√N) = ζ·L(s,χ_N). The shakeable-leaf conjecture: a semiprime N = pq is a marginally bound leaf; where it drops defines the factors. OPEN; the Fermat half is not in the heartbeat wave of §6.

## 8. Reading guide: what is proven vs what is reading

| Claim | Status |
|---|---|
| partition B+R+M = 1; T + L = ⌊x⌋−1; W_L = −W_T up to 2/√x | identity |
| H(u) is the zero expansion of W_T | ESTABLISHED (von Mangoldt) |
| 65 primes decide primality below 317² | ESTABLISHED/exhaustively verified |
| stragglers' composite share 1−1/(Uω(U)) | OURS-CALCULATED to ≈ 1 % at y = 313 (U 2.1–2.8) |
| Telperion/Laurelin interaction = H | READING |
| ADD:SCALE:SIGN = phase / amplitude-frequency / half-turn | READING |
| divisor involution ↔ zero quadruple | READING (structural) |
| n-ball volume ↔ factor curve | tested, naive form fails (§9); OPEN for other forms |

## 9. Test log: the n-ball hypothesis (2026-11-08)
Candidate: the curve that "gets larger then trails off" is the volume of the unit n-ball, V_n = π^{n/2}/Γ(n/2+1), peak at n = 5 (V = 5.264). Compared with the Laurelin book (Dickman ρ, monotone decreasing, ρ = 1 on [0,1]): ln ρ(u) against ln V_{k u}(1) for k = 1…4 on u ∈ [2, 8] differ by up to e^{19}, e^{16}, e^{11}, e^{5}; ρ has no maximum. The Telperion book (Buchstab ω → e^{−γ}) does not decay at all, so it cannot match either. Both ρ and V decay like Γ⁻¹-type factors, so the leading u^{−u} behaviour is shared; the details differ. **Result: the naive form (R = 1, n = ku) does not match; not excluded: variable-radius balls, n = n(τ), or a different quantity (the torus radius R_K, or the number of wheels).** The straggler equation of §4 is the equation for the factors still present after 65 primes.

## 10. Anomalies logged, not repaired
Archimedes Screw `zeros_upto` off by ≈ 2.3 past n = 50; Σ_RB balance independent of σ; the "σ = ½ saddle λ ≈ 6" (oblique gear) definition not located; `1e0657.clustergal.dat` (Bullet Cluster) is a saved 404 page — unrelated to this repo; ten-item anomaly log in notebook 03.

## 11. Where each item lives
`rhp/` modules: sphere, unsmap, wrap, extinction, nodal, chladni, coil, standing, heartbeat, wave, primepath, wheel65, diverge, zerorb, flow, bigpath, clock, echo, harness, lehmer, balance, primegame. Tests: `tests/` (107). Notebooks 01–06. Paper §7.4–7.20. Atlas: `rhp/atlas/` (see `docs/IMAGE_CATALOG.md`). Experiments log: `experiments/LOG_2026-11-07.md` (items 1–17). Engines: `GenerationalLineage/engine`, `AbrikosovTree/engine/two_trees.py`, `ValaQuenta/modules/{h_rb_hat,archimedes_screw,bao_mass_gap,inversion,spectral_primes,berry_keating,box_kite,oblique_gear,prime_gauge_field,add_scale_sign}`. Wiki: Ainulindale 47, 51, 58, 73, 75, 96, 126, 130.
