# RiemannHypothesisProof — image catalog (2026-11-08)

Every image in the repository, with a full accounting of what each section (panel) shows, how to read it, the maths behind it, the numbers that matter, and what it does **not** show.
Paths are relative to the repository root. `rhp/atlas_riemann_fermat/*` are redrawn by `python -m rhp.atlas_riemann_fermat.<name>` (needs the repo venv; `chladni`, `wheels` take ≈ 1 min).

**Labels used throughout.** ESTABLISHED = textbook/literature. OURS-CALCULATED = derived or measured here with code and tests. READING = interpretation, not a theorem. OPEN = untested.
**Conventions.** u = ln x (size, the Archimedes Screw's own axis). τ = ln y (sieve time, the cutoff). s = σ + it; σ = ½ is the critical line = the equator of the Cayley sphere w = (s−1)/s (north pole s = 1, south pole s = 0).
Telperion = primes / sieve / smallest prime factor (spf) / death. Laurelin = composites / un-sieve / greatest prime factor (gpf) / birth. Colours: blue/north = Telperion side (σ>½), orange/south = Laurelin side, red = highlight, grey = reference.
Nothing in any figure is evidence about RH; a picture here shows a space or a machinery.

Index: 1 `equator` · 2 `polarisation` · 3 `wrap` · 4 `wheels` + `wheels_animated` · 5 `twins` · 6 `flow` · 7 `flow_player.html` · 8 `net` · 9 `nodal` · 10 `diverge` · 11 `uns_map` · 12 `chladni` · 13 `spirals` · 14 `heartbeat` · 15 `wave` · 16–18 `primepath_*` · 19–24 `figures/*` (older) · 25–28 `archive/…`.

---

## 1. `rhp/atlas_riemann_fermat/equator.svg` — the critical line as the equator
Code: `rhp/sphere.py`. Notebook 04. Paper §7.10 (equatorial construction).
- **a — the Cayley sphere.** The unit sphere with north pole s = 1 (Λ residue +1) and south pole s = 0 (residue −1). Faint grey = parallels/meridians; thick red = the equator |w| = 1 ⇔ σ = ½; dashed = far side. Black dots = the first 80 zeros placed at longitude θ(t) mod 2π (open = far side).
  They scatter about two antipodal carrier points (sd 0.25 half-turn; a uniform scatter would be 0.29): a weak bunching, **not** two spots. The pair marked (t = 20) shows s ↦ 1−s̄: same longitude, opposite latitude. Label: coordinate change (ESTABLISHED) + measured scatter (OURS-CALCULATED).
- **b — the same thing flat.** The strip 0 ≤ σ ≤ 1, t ≤ 80. Blue = north (Euler product converges for σ>1), orange = south (its mirror), red line σ = ½ with the first 40 zeros. Top axis re-introduces π: latitude π(σ−½) in degrees (−90…+90). σ = 0 and σ = 1 are the two poles; the strip's edges collapse to points on the sphere.
- **c — π factored out of the phase.** θ(t) = Im ln Γ(¼+it/2) − (t/2) ln π. Blue = Γ part (carries the t ln t growth); orange = π part (a constant-rate rotation); black = their sum = mpmath siegeltheta. Max |sum − θ| = 1.4·10⁻¹⁴ (so the split is exact).
- **d — what the equator forgets.** Left axis (red): Cayley longitude/π of ½+it: the whole infinite line wraps ONCE around the circle and then crowds at w = 1. Right axis: Hardy winding θ(t)/π — one half-turn per zero (up to 3066 at t = 3600). The circle is the flat projection; the winding is the zero index.
- Not shown: any statement about where zeros are; this is a change of coordinates.

## 2. `rhp/atlas_riemann_fermat/polarisation.svg` — two counter-rotating waves per term
Code: `rhp/sphere.py` (abs_chi, amplitude_ratio, polarisation, ellipse_axes). Notebook 04.
Each term n of the approximate functional equation is a pair: forward n^{−s} and mirror χ(s) n^{−(1−s)}. Amplitude ratio r = |χ|·n^{2σ−1}; the two phasors trace an ellipse of axial ratio b/a = |1−r|/(1+r) (0 = linear = pure standing wave, 1 = circular = pure travelling wave).
- **a — term n = 1, t = 100.** Five small ellipses at σ = 0.5, 0.4, 0.3, 0.15, 0.0 with r = 1.00, 1.32, 1.74, 2.63, 3.99 and b/a = 0.00, 0.14, 0.27, 0.45, 0.60. On σ = ½ it is a line (standing wave).
- **b, c — b/a over the (σ, n) plane at t = 100 and t = 1000.** Pale = standing. Two pale sets form a CROSS: the equator (every term) and the horizontal line n = √(t/2π) (3.99 at t = 100; 12.62 at t = 1000) — the Riemann–Siegel cut-off, where forward and mirror balance at every latitude.
- Label: ESTABLISHED (approximate functional equation, Riemann–Siegel) + OURS-CALCULATED (the plots).

## 3. `rhp/atlas_riemann_fermat/wrap.svg` — the primes on the two springs
Code: `rhp/wrap.py`. Notebook 04. Two different springs, kept apart.
- **a, b — the sieve spring.** The integers wrapped on circles of circumference 6 and 30; primes marked; every twin pair (p, p+2), p > 3, straddles the 0 spot of the circumference-6 circle (checked: all of them).
- **c — allowed-centre fraction on the primorial circle** at each level P#: Π_{5≤p≤P}(p−2)/P#: the 0 spot thins at every level.
- **d — the wheel spring.** Prime wheel p is the rotation e^{−it ln p}, radius p^{−σ}. Twin wheels drift apart in phase by t·ln(1+2/p) ≈ 2t/p; the GAP is carried by the beat, not by the values.
- **e — the zeros' spring tightens:** height of one 2π turn = 2π/θ′(t) against t (log axis): it shrinks like 1/ln t.
- **f — the primes' spikes:** prime spikes of the explicit formula sit at t = ln n and their density per unit t grows like eᵗ/t (log axis). Zeros tighten logarithmically; primes crowd exponentially: that is the tension.

## 4. `rhp/atlas_riemann_fermat/wheels.svg` (static) and `wheels_animated.svg` (SMIL animation)
Code: `rhp/wrap.py` (log_abs_euler, euler_grip). ≈ 1 min.
- **a — error of K prime wheels against |ζ|, by latitude σ.** Far north the primes win (the error falls with K); at the equator they stop helping after ~10² wheels (the Euler product does not converge on σ = ½).
- **b — K\*(σ)**, the number of wheels that does best, against latitude.
- **c — the first zero emerging from the wheels:** ln|P_K| against ln|ζ| near γ₁ = 14.1347.
- **Animated.** The first six prime wheels (radius p^{−½}, rate ln p), tip-to-tail, with the trace of their sum.

## 5. `rhp/atlas_riemann_fermat/twins.svg` — the two books, the circles' twin constant, the clock–wheel resonance
Code: `rhp/wrap.py`, `rhp/uns.py`.
- **a — Telperion × Laurelin.** Every n ≤ 10⁵ at (π(spf), π(gpf)): when the sieve first kills it, when the wheels first complete it. Primes sit on the diagonal.
- **b — the circles' ratio Π p(p−2)/(p−1)²** against the largest circle P: converges to the Hardy–Littlewood twin constant C₂ = 0.66016.
- **c — twin pairs below x** counted against 2C₂ Σ 1/ln²x (a heuristic prediction, not a proof).
- **d — gaps between consecutive twin pairs** in units of the local mean (ln p)²/(2C₂), against an exponential law.
- **e — wheel n against the clock:** θ(t) − t ln n is stationary at t = 2πn² (n = √(t/2π)): the one wheel at rest relative to the clock at each height.

## 6. `rhp/atlas_riemann_fermat/flow.svg` — the extinction flow and its Laplace dual
Code: `rhp/extinction.py`, `rhp/wrap.py`. Notebook 05. Paper §7.16.
- **a** survivor fraction of the sieve against the cutoff time τ = ln y, with Mertens' e^{−γ}/τ. After the first 65 primes (τ = 5.746) the fraction is 0.0965 (Mertens 0.0977).
- **b — Telperion.** Φ(N,y)·ln y/N against U = ln N/ln y for four N, on the Buchstab function ω(U) ((Uω)′ = ω(U−1), ω = 1/U on [1,2], → e^{−γ}). Finite-N curves sit off ω by the known O(1/ln y).
- **c — Laurelin.** Ψ(N,y)/N against U for four N, on the Dickman function ρ(U) (Uρ′ = −ρ(U−1), ρ = 1 on [0,1], → 0 superexponentially).
- **d — the gain of the sieve book in the dual plane:** ln|1+ω̂(s′)| = Re E₁(s′), s′ = σ′+it, with a log branch at s′ = 0 (= ζ's pole, the north).
- **e — reciprocity** s′ρ̂(s′)(1+ω̂(s′)) = 1 checked from the numeric ω and ρ along σ′ = 1; residual ~10⁻⁴ is the grid's first-order error (closed forms agree to 10⁻¹²). The two books are reciprocal under the Laplace transform in U — the same relation as ζ·(1/ζ) = 1. Label: OURS-CALCULATED (may turn out to be known).

## 7. `rhp/atlas_riemann_fermat/flow_player.html` — the extinction flow, interactive
Self-contained page, no network. Drag τ = ln y (or press play). Panel 1: every n ≤ N at (u = ln n, status): blue = alive (spf > y), grey = dead (spf ≤ y), orange ring = smooth (gpf ≤ y); primes (black) die and are born at the same instant u = τ.
Panel 2: survivor fraction against Mertens with a marker at the current τ. This is "how long it takes for composites to fall off the tree": integer n dies at τ = ln spf(n) and is born at τ = ln gpf(n); for 10: dies at ln 2 = 0.693, born at ln 5 = 1.609, span 0.916 e-folds (mean span over n ≤ 10⁶: 6.55).

## 8. `rhp/atlas_riemann_fermat/net.svg` — the six faces of the process (a hypercube net)
Code: `rhp/wrap.py`, `rhp/extinction.py`. Notebook 05. Axes: u = ln n, τ = ln y, σ′ = σ−1, t. Each square is one pair of axes; only the six together show the process (observer-position: change one parameter at a time).
Top row (squares containing u): (u,τ), (u,σ′), (u,t). Bottom row: (τ,σ′), (τ,t), (σ′,t).
- **(u,τ)** each integer n ≤ 3000 at its two clocks: death ln spf(n) (blue), birth ln gpf(n) (orange); primes on the diagonal.
- **(u,σ′)** energy of the first n ≤ eᵘ integers Σ n^{−2σ} (log₁₀): finite for σ>½, logarithmic at ½, divergent below.
- **(u,t)** the dual lattices: primes are VERTICAL lines u = ln p (the leaf-drops of the Screw), the first zeros are HORIZONTAL lines t = γ.
- **(τ,σ′)** sieve survivors Π_{p≤e^τ}(1 − p^{−1−σ′}) (log₁₀): settles for σ′>0, falls like 1/τ at σ′ = 0 (Mertens), faster below.
- **(τ,t)** phase of the Euler partial product P(½+it) with the wheels up to e^τ switched on — the spring seen from the side.
- **(σ′,t)** gain of the sieve book Re E₁(s′) (the dual plane of flow.svg).

## 9. `rhp/atlas_riemann_fermat/nodal.svg` — the nodal map of ξ on both hemispheres
Code: `rhp/nodal.py`, `rhp/sphere.py`. Notebook 05. Red = Re ξ = 0, blue = Im ξ = 0.
- **a — the strip σ ∈ [−1.5, 2.5], t ∈ [0, 50].** The blue set contains the whole equator σ = ½ (Im ξ = 0 there); the red RUNGS cross the equator only at the zeros, level with t = γ. Mirror symmetric about σ = ½ (ξ(s) = ξ(1−s)). For t ≤ 50: 10 zero rungs, 9 saddle rungs; hyperbolic saddles (Hessian ±λ); no off-line Re = Im = 0 crossing in σ∈[−0.6,1.6], t∈[10,35] (OURS-CALCULATED).
- **b — zoom on the first five zeros** with the heartbeat's extrema ◆ (F′ = 0): saddles of the analytic landscape on the equator; exactly one saddle between consecutive zeros.
- **c — the heartbeat F(t) = ξ(½+it)** with its zeros (rungs) and extrema (saddles).
- **d — log₁₀|χ(σ+it)|** (ζ(s) = χ(s)ζ(1−s)): zero on the equator, the weighting that turns the north's picture into the south's.
- **e — ζ itself on the real side:** Re ζ = 0 and Im ζ = 0 for σ ∈ [−9, 2], t ∈ [0, 14] — the trivial zeros at −2, −4, −6, −8 are on the real axis in the south; ξ has none there (Γ(s/2) cancels them), so they belong to χ, not to the equatorial structure.

## 10. `rhp/atlas_riemann_fermat/diverge.svg` — the divergence detector
Code: `rhp/diverge.py`. Notebook 05 §8. Tests: control, injection, limits. Explicit formula E(u) = Σ_ρ e^{iγu}/ρ (u = ln x). An off-line zero adds a DILATE e^{δu} and its mirror e^{−δu} (a cosh): a length change where only rotations are allowed (generational-lineage emergence signature 1).
- **a — CONTROL on the real ψ (x ≤ 10⁷).** Slope of ln A against u at the first 20 true zeros; scatter about 0 is the false-alarm floor (sd 0.033, band ±3 sd).
- **b — injected off-line quadruple at γ₁** for δ = 0, 0.1, 0.2, 0.3: lock-in amplitude along u; the mirror partner makes it a cosh, not a line.
- **c — recovered slope against δ** (dots) with the 3-sd floor, and the smallest δ seen at reach u_max = 11, 13, 16: 0.52 / 0.28 / 0.13.
- **d — frequency scan:** amplitude against γ for the real ψ (the true zeros are the peaks) and with δ = 0.3 injected at γ = 17.5 (not a zero): it is found there.
- Honest limit: a null result on the real ψ says only that no growth ≥ the floor exists within u ≤ u_max.

## 11. `rhp/atlas_riemann_fermat/uns_map.svg` — flat Cartesian → UNS, with no numbers in it
Code: `rhp/unsmap.py`. Notebook 06. Artefact detectors.
- **a** lines σ = c and a circle in the plane (the chart s = σ+it). **b** their images on the sphere (Cayley, then stereographic): all circles; σ = ½ (thick) is the great circle = the equator; the others are small circles.
- **c** a cusp in the plane. **d** its image: still a cusp (the map is conformal, so a genuine cusp survives; flat and sphere speeds vanish together).
- **e** a back-and-forth trace (reversals) that is a UNIFORM ROTATION seen edge-on: the Hilbert companion supplies the missing coordinate (constant amplitude, linear phase).
- **f** three spiral laws on test spirals (R² of each law's straight line): Archimedean r ~ φ (pancake), logarithmic ln r ~ φ (loxodrome), Euler/Cornu κ ~ s.

## 12. `rhp/atlas_riemann_fermat/chladni.svg` — the two hemispheres as one oscillation
Code: `rhp/chladni.py`, `rhp/nodal.py`. ≈ 1 min.
- **a — nodal lines on the sphere** in the chart (latitude π(σ−½), longitude θ(t) mod 2π): Re ξ = 0 (red) and Im ξ = 0 (blue) are meridian arcs — a ladder wrapped around the axis, stepping by the lobe angle (NOT a double helix — corrected).
- **b — parity at t = 30:** Re ξ even (free, ∂σ = 0 at the equator), Im ξ odd (held, = 0 there): two polarisations with different boundary conditions.
- **c — the flux,** scale-free ∂σ ln|ξ| = J/|ξ|², J = ½∂σ|ξ|² = Im(conj ξ ∂tξ): zero on the equator, opposite in the two hemispheres, positive everywhere in the north in this window. Exact identity: ∂σ ln|ξ| = θ′(t) + Re[1/s + 1/(s−1)] + ½[Re ψ(s/2) − Re ψ(¼+it/2)] + Re ζ′/ζ.
- **d — inject an off-line zero pair** (δ = 0.08 at γ = 25): the flux reverses in a pocket — the footprint of an off-line zero (a test of the detector, not a search).
- **e — the lobe:** angle swept by the Hardy phase between consecutive zeros: mean π, spread = the GUE zero-spacing distribution × π (sd ≈ 1.18). One of the candidate homes of "π" in the pathway.

## 13. `rhp/atlas_riemann_fermat/spirals.svg` — not flat, not static: conical helices
Code: `rhp/coil.py`, `rhp/unsmap.py`. Notebook 06.
- **a — ζ(½+it) in the cylindrical UNS chart** (angle arg ζ, radius (2/π)arctan|ζ|, height t): the pathway THREADS the t-axis (ζ = 0) at each zero.
- **b — the same at σ = ½, ¾, 1:** off the equator the axis is a hole that widens; min distance from the axis 0.002 / 0.076 / 0.180 / 0.326 (σ = .5 / .6 / .75 / 1); pitch tightens with t in all three.
- **c — partial sums Σ n^{−s} about ζ(s), height u = ln N:** a CONE at σ = ½ (flare e^{u/2}), a CYLINDER at σ = 1, a FUNNEL at σ = 3/2 — logarithmic spirals with d ln r/dφ = −(1−σ)/t (predicted, verified at four latitudes).
- **d — the prime side:** the conical helices of the explicit formula (first zeros, flare e^{u/2}) and their sum; prime heights u = ln p marked on the axis.
- **e — a flat pancake coil** (POE: 17 turns, 0.8 mm pitch, taps at 1/3/5/17 turns ↔ zeros 2/6/10/34 at γ = 21.02/37.59/49.77/111.03) with the zeros as beads, one turn per two zeros: with radius = zero index the Hardy angle is Archimedean (r² > 0.9999 over 300 zeros). It keeps angle and index; it loses height and radius.
- **f — reading.**

## 14. `rhp/atlas_riemann_fermat/heartbeat.svg` — what the standing wave contains under its curve
Code: `rhp/heartbeat.py`, `rhp/diverge.py`. Paper §7.18. **Not** the "area π" of the unit circle: that is the normalisation (wiki/96: ∫₀^{2π} sin² = ∫ cos² = π).
- **a — the unit circle's area π split:** the zero circles (π/|ρ|² each, all zeros) hold 2.31 % (ESTABLISHED Hadamard sum Σ_ρ 1/ρ = 1 + γ_E/2 − ½ ln 4π = 0.0230957); the Γ/π ground state the rest (97.69 %).
- **b — who holds the zero-circle area:** first zero alone 21.6 %, first 10: 58.6 %, the 34 on the pancake coil: 76.0 %, first 100: 86.5 %; 3066 zeros + tail estimate = 0.02309574 vs exact 0.02309571.
- **c — mean power** of the real wave on the exact ψ, in windows of u, against the Parseval value 2 + γ_E − ln 4π = 0.04619 (measured 0.0453–0.0466).
- **d — one standing wave as two counter-rotating circles:** equal radii 1/|ρ|, the sum is real; off the line (δ = 0.15, dashed) both radii grow as e^{δu}.
- **e — ADD / SCALE / SIGN on a tone:** a step Δu in ln x multiplies the tone by e^{(½+iγ)Δu} = a spiral similarity (SCALE e^{Δu/2}, rotation γΔu); SIGN = conjugation (the reflection ρ ↔ ρ̄). At γ₁, Δu = 1: scale 1.6487, rotation 1.5684 rad.
- **f — reading.** Honest limit: no single tone holds > 0.5 % of π.

## 15. `rhp/atlas_riemann_fermat/wave.svg` — the Heartbeat waveform from the two sets, and its wind
Code: `rhp/wave.py`. Derivation: `docs/heartbeat_waveform_derivation_2026-11-08.md`. Paper §7.19. Foundation function H(u) = −2 Σ_{γ>0} cos(γu − arg ρ)/|ρ|, Telperion carries +H, Laurelin −H.
- **a — W_T (primes) and −W_L (composites) from the SETS** (no zeros used): W_T = (ψ(x)−x)/√x, W_L = (⌊x⌋−1−ψ(x))/√x. They coincide: |W_T + W_L| ≤ 2/√x (measured 0.098 ≤ 0.0996). An exact identity, because the two trees partition the integers.
- **b — the sieve wave against the sum of the first 10 / 100 / 3066 zero tones.** Correlation 0.771 / 0.931 / 0.984 (1000) / 0.993 (3066).
- **c — the wind:** |W| (grey) and its Hilbert envelope (500 tones, red) against the rms 0.217 (predicted √(2+γ_E−ln4π) = 0.215). Gust factor 3.66; envelope sd/mean 0.48; 9.5 % of the range below 0.1 (calm), 16 % above 0.4. Skew −0.02, excess kurtosis −0.25.
- **d — one factor's step ln p/√p** (a prime adds ln p to ψ at x = p, in W units) against the rms wind. It peaks at p = e² (7.389) with value 2/e = 0.736 — where the engine's Mingling sits. The largest envelope seen (0.74) is the same size (observation; no mechanism claimed).
- **e — fraction of u-time that |W| ≥ one factor's step,** against p: p=2: 2 %, p=5–11: 0.02 %, p=101: 3 %, p=1009: 33 %, p=10⁴: 68 %, p=10⁶: 95 %.
- **f — zero tones needed to resolve a prime's jump** (fraction of the true ln p recovered): p=11 at once; p=101: 0.34/0.59/0.94 at K=5/100/1000; p=307: 0.12/0.19/1.11 (Gibbs overshoot) — bandwidth T ≳ a few × p.

## 16. `rhp/atlas_riemann_fermat/primepath_flat.svg` — the pathway of the first K primes, flat
Code: `rhp/primepath.py`. Paper §7.20; `docs/prime_pathway_2026-11-08.md`. The path is log P_K(s) = Σ_{k≤K} −ln(1−p_k^{−s}) at s = ½+it, the logarithm of the partial Euler product; red = first 65 primes (up to 313), blue = the rest, orange dot = the end.
- **Row 1 (t = γ₁ = 14.1347, the first zero), columns K = 65 / 10⁴ / 10⁵ / 10⁶:** |log P| = 2.6 / 3.8 / 8.9 / 13.2. After a jagged tangle from the first 65 primes the path becomes a smooth OUTWARD spiral centred near (−2.7, 0) that never closes.
- **Row 2 (t = 3):** |log P| = 1.6 / 9.5 / 27.3 / 78.5 — the same law, clearer. The prime sum tracks the prime-number-theorem part Ei((1−s) ln x): ratio 1.001 at t = 3, K = 10⁵. Radius grows like √x/(|1−s| ln x); angle turns at −t ln x: a logarithmic spiral.
- Read-across: the 65 primes are the TIP of the spiral; the spiral does not stop at prime 65.

## 17. `rhp/atlas_riemann_fermat/primepath_3axis.svg` — the same path in three axes
Axes (Re log P, Im log P, ln p) for 10⁴ primes at t = γ₁; ln p is the sieve clock τ.
- **Isometric / top (Re, Im) / front (Re, ln p) / side (Im, ln p):** a funnel (conical helix) opening with ln p; red first 65 primes = the jagged tip. Orthographic projections.
- **e — |Σ p^{−s}| against K** (log–log), t = γ₁ (red) and t = 3 (blue); dashed = the PNT smooth part. The path keeps growing.
- **f — the sphere radius R_K = √Σ1/p holding the torus of K wheels** (each wheel is a circle in its own orthogonal complex plane; K wheels = a point of a flat torus T^K ⊂ C^K): 1.42 (K=65), 1.65 (10⁴), 1.70 (10⁵), 1.75 (10⁶). It grows like √ln ln p — without bound, absurdly slowly.

## 18. `rhp/atlas_riemann_fermat/primepath_uns.svg` — the partial product on the Riemann sphere
P_K = exp(log P_K) stereographically projected: |P| = 1 is the equator, P = 0 the north pole, ∞ the south pole. Panels 1–3: 10⁴ primes at t = γ₁ — isometric, from the north pole, side with the equator edge-on (red = first 65 primes). Over 10⁴ primes the product sits near the north pole (|P| ≈ 0.04–0.08).
Panel 4: ln|P_K| against K to 10⁶: oscillates with growing amplitude (±19 by 10⁶), crossing |P| = 1 17 times; the product does NOT converge to ζ(½+iγ₁) = 0 within 10⁶ primes (the Euler product on the critical line is not convergent in the usual sense — ESTABLISHED).
The s-point itself stays on the equator for every K; the primes act through the VALUE P_K.

---
## Older figures (`figures/`) — made before the Riemann–Fermat Atlas; notebooks 01–03

## 19. `figures/flattening_map.svg` — notebook 01, paper §7.4
Orientation: log-polar in the ζ-plane (radius ρ = log₁₀|ζ| + 3, angle arg ζ) with t explicit; aperture: σ = ½, continuous in t, t ≤ 50.
- **A flat:** Re ζ against Im ζ (t not drawn); hub: the 10 zeros are one point; ζ(½) = −1.4604; dots on the real axis = Gram points gₙ (θ = nπ). **B UNS end-on:** the axis is the zero set; ρ = 0 at the floor (a rendering choice). **C** log₁₀|ζ| against t (zeros = spikes to the floor).
- **D the lost clock:** −θ(t)/π (Gram points at integers); turning point θ′ = 0 at t = 2π; rate −θ′ = −½ ln(t/2π).
- **E, F, G the spring:** front [90,0,0] height = t; iso [55,0,25]; top [0,0,0] (= panel B). Colour = t (violet 0 → yellow 50). Orthographic.

## 20. `figures/machinery.svg` — notebook 02, paper §7.6
- **a** pathway ζ(σ+it), t ≤ 50, flat, at σ = 1.5, 1.1, 0.8, 0.65, 0.5. **b — the log plane: wheels add:** log ζ(σ=1.5) against the sum of the first 3 and first 10 prime wheels (max error over t ≤ 60 after K = 1, 3, 10, 100 wheels: 0.52, 0.22, 0.070, 0.011 at σ = 1.5; 1.7, 1.2, 0.79, 0.43 at σ = 1.1).
- **c — σ = ½: |A| = |χB|** exactly (the forward wheels A and the mirror χB have equal modulus, max difference 3·10⁻¹⁴ on 25 points in t∈[20,120]). **f — |χB|/|A|:** balance on the line; off the line one half dominates (median 1.17 at σ = 0.4, 0.85 at σ = 0.6).
- **d — primes from zeros:** exact ψ(x) − x against the first 10 / 100 / 3066 zeros (rms error on [2,60]: 0.74, 0.37, 0.095). **e — a continuous sieve:** a Gaussian e^{−(γw)²/2} on each zero gives a smoothed sieve that sharpens to the discrete jumps as w → 0 (w = 0.02, 0.005 shown).

## 21. `figures/flow.svg` — notebook 02 (the ξ′/ξ flow; module `rhp/flow.py`)
Complex potential Φ = log ξ(s); velocity (u, v) = (Re L, −Im L), L = ξ′/ξ. Zeros are sources + vortices.
- **Left — margin of outflow:** log₁₀ u/(σ−½) (finite = same sign everywhere); black dots = zeros. **Middle — direction of (u,v):** u mirrors across σ = ½ (odd), v does not (even). **Right — histogram** of log₁₀(u/(σ−½)): > −∞ everywhere ⇒ outflow (a calculation on the grid, not a proof). Gauss's-law accounting: argument principle N = θ/π + 1 + S.

## 22. `figures/long_pathway.svg` — notebook 02, 3065 zeros
Riemann–Siegel (first correction), zeros polished with mpmath, t = 20…3600; orthographic; panels A, B overlay 1532 full turns.
- **A flat; B UNS end-on** (ρ = log₁₀|ζ| + 2.5, angle arg ζ). **C the hub fan:** S(γₙ) mod 1, one dot per zero; direction at the origin = π/2 + π·S(γₙ) (mod π). **D histogram of that fan:** mean +0.004, sd 0.254. **E S(t)** = N(t) − θ(t)/π − 1 between zeros: sd 0.20, max |S| 0.54. **F spring, front:** height = t. **G spring, iso:** [55,0,25], t = 20…300. Colour = t (violet 20 → yellow 3600).

## 23. `figures/clock_echo_harness.svg` — notebook 03, paper §7.9
- **a — the θ-clock:** Z(t) against θ/π with the n = 1 carrier 2cos θ: zeros at half-integers of θ/π, jittered by S. **b —** zero − carrier in spacings (3066 zeros within 0.81 spacing of the carrier, mean 0.22). **c — Gram points:** 1 = negative real-axis crossing (197 of 3066 below 3600 = 6.4 %; first at n = 126, t = 282.455).
- **d — echoes at t = ln(prime power):** Σ cos(γt) over 3066 zeros (Γ = 300) against the explicit formula: spikes at the primes. **e — unfolded spacings:** zeros against GUE (unfolded) and Poisson. **f — PRIMEGAME registers,** first 600 steps (rows = prime).

## 24. `figures/balance_margins_primegame.svg` — notebook 03
- **Left — off-line balance crossings of the truncated halves** (curves, σ against t) with the zeros (black): where the truncated forward and mirror sums balance off the line. **Middle — Lehmer margins:** relative extremum height of Z in each gap. **Right — PRIMEGAME (Conway):** steps against p, slope 2.89. A ten-item anomaly log accompanies notebook 03.

---
## Archive (`archive/2026-11-03_pre-rewrite/`) — pre-rewrite, kept for lineage; NOT authoritative
## 25. `figures/05a_redblue_balance.png`
RedBlue balance E_Red − E_Blue along the isoenergy surface xp = 1 (x from 0.5 to 5): the curve is above zero only near x ≈ 1.3–2.8 and strongly negative outside; shading marks where Red or Blue dominates. 0_RB's Red = xp (Berry–Keating), Blue = ½p² + ℘ (Weierstrass).
## 26. `figures/05b_forced_sigma.png`
forced_sigma: iteration σ ← (Fσ + B(1−σ))/(F+B), F = e^{−σE}, B = e^{−(1−σ)E}: every start (0.01…0.99) converges to ½. The cleaner statement (canonical maths FORCED_SIGMA): h(s) = ln F − ln B = E(1−2s) is linear, so one Newton step reaches ½ exactly.
## 27. `figures/05c_fermat_extinction.png`
Gradient descent minimising |E_Red − E_Blue|² with the real ℘′: 25 random starts (red dots) flow to a locus (orange survivors) near x ≈ 1.2–1.4 / 3.0. "Extinction" here is a dynamical flow on the 0_RB balance surface, not the sieve extinction of flow.svg.
## 28. `images/Gemini_Generated_Image_Riemann_Proof.png`
An AI-generated illustration ("The Riemann Sphere — Geometry of the Zeta Resonance"), 768×1376, used as the old README header. It is an illustration, not data; the summarised standing-wave text from the same source is stale and was used only as a pointer.


## gemini_worm.svg  (rhp/atlas_riemann_fermat/gemini_worm.py)
Eight panels, 2 × 4. a: the oldest Gemini sim's flow z ↦ 1 + 1/z in polar (13 points converging on r = φ). b: Revised_2 Inertial (θ = n, r = 0.02·n·t, t = 1). c: Revised_2 Metric Swap (r = 1/(0.02·n·t + 0.1), θ = (2/π)n + π/2). d: the Transition sim: the inertial path (blue) crosses r = 1 and is inverted (red) with the 0.0007 gap. e–g: Revised_2 Heartbeat at t = 1, 3, 7 (600 points, polyline through consecutive n; rosettes of ≈ π²/t arms). h: the ζ path log P_K at s = ½ + iγ₁, K = 2000, first 65 primes in red. Reading: constant-angle √n spiral (Gemini) against slowing-angle p^{-1/2} spiral (ζ).


## gemini_original.svg  (rhp/atlas_riemann_fermat/gemini_original.py)
2 × 5. Top row, columns 1-4: phase 2 of the original event_horizon_sim.py (the worm) at cycle time t = 0.1, 0.3, 0.6, 1.0 — a spiral head with a tail that lengthens and wiggles. Bottom row, columns 1-4: phase 3 (the rosette) at t = 0.2, π/6 (k = 1/3), π/4 (k = 1/2), 1.0 (k = 2/π) — a √θ spiral with radial modulation; symmetric lattice at k = 1/2. Top right: line spectrum of the worm's phasor against Bessel |J_n(5)| (phase modulation). Bottom right: rosette sidebands, carrier 1 and 1 ± k at 0.1 for k = 1/3, 1/2, 2/π (amplitude modulation).


## boundary_engine/  (directory of PNG and GIF; `rhp/boundary_engine.py --export`)
`poster_all_patterns.png` — the whole 5 × 5 pattern grid (rows: original, Claude's five phases, Revised_2, φ-flow/crossing/wheel family, ζ whole/Telperion/Laurelin). `strip_{O2,O3,W1,Z1,R3}.png` — each pattern at τ = 0.1, 0.3, 0.6, 1.0. `transition_<A>_<B>_{polar,cartesian,wheel}.png` — A → B at α = 0, 0.25, 0.5, 0.75, 1 for O2→O3, O1→O2, R2→R3, G1→Z1, O2→Z1. `motion_worm.gif`, `motion_rosette.gif`, `motion_transition_wheel_O2_O3.gif`, `motion_zeta_height_sweep.gif` (log P_K at t = 3…30). Reading: wheel-state LERP keeps structure between worm and rosette; Gemini's polar LERP tangles at α ≈ 0.25.
