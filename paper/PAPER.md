# Factoring π out of the Zeta Function: the critical line as an equatorial great circle, and the places where ½ emerges

**DRAFT 0.1 — 2026-11-03 — not for submission.** Written from the code and measurements that exist today; nothing here is carried over from the archived drafts (`archive/2026-11-03_pre-rewrite/`) unless it was re-run.

**Author:** Cody Michael Allison (Michael Rendier) · the.wandering.god@gmail.com · ORCID 0009-0007-7239-6760

---

## Note on how this paper was made

This paper was written with an assistant: Claude (Anthropic), model `claude-sonnet-5-5`, a transformer, running under the project's `cs-paper-code-conventions` skill (§8 of the skill: papers written with an assistant). The author is an independent researcher, not a number theorist by training; he asked the assistant to write the paper from the code and the measurements in this repository, and to keep every claim at the strength its evidence supports. As a safeguard — not a guarantee — he does not tell the assistant the result he hopes for. The paper is a **code paper**: it shows what runs, labels what is established, what is the author's, and what is conjecture, and it states one open step as open. The suffix `:CALCULATED` marks a result that this repository's own code produces; it does not mean *derived* or *proved*. Where this repository's earlier drafts claimed more than their own text supported, §10 says so.

## Abstract

The completed zeta function ξ(s) = ½ s(s−1) π^(−s/2) Γ(s/2) ζ(s) satisfies ξ(s) = ξ(1−s) only because the factor π^(−s/2) is present: remove it and the reflection symmetry about Re(s) = ½ is destroyed. We take this literally. Dividing the strip coordinate by π (θ = πσ) maps the critical strip 0 < σ < 1 onto the polar range 0 < θ < π of a sphere, the reflection s ↦ 1−s onto the equatorial reflection θ ↦ π−θ, and the critical line onto the equator, with the offset ½ = (π/2)/π. We show by computation what this map does and does not do (it is an exact change of variable; it does not locate any zero), and we collect, each with the code that checks it, the places where ½ appears *emergent* — not constructed and not composite — in ordinary mathematics and in the author's "universal native space" (UNS) reading: the self-dual point of s ↦ 1−s, the unit-diameter circle (r = ½, C = π, area π/4), the Amplitude-Lagrangian minimum, the in-phase/quadrature balance (tan = 1), the half-integer offset of the zero-counting phase, and the equal-envelope condition of the explicit formula. We state a reading of ζ as a Lagrangian, a map of what the flat plot of the zeta spiral overlays and what it keeps when the same curve is read in the UNS log-polar chart, and a conjecture that the prime wavelength (zero spacing) and the factor wavelength (λ = p) are orthogonal. The Riemann Hypothesis is not proved here; the step that would prove it is the Berry–Keating self-adjointness conjecture, which remains open, and §9 states exactly which parts of the argument are proven, which are imported from the literature, and where the gap is.

## 1. Introduction

The Riemann Hypothesis asserts that every non-trivial zero of ζ(s) has real part ½. The number ½ is usually presented as a feature of the problem. We ask a narrower question, one the code in this repository can answer: **where does ½ come from, and does it come from one place or from many?**

One fact is standard but rarely read literally. The completed function is symmetric about ½ only because of π^(−s/2). The reflection s ↦ 1−s has the fixed line Re(s) = ½ for any function; what the factor π^(−s/2) does is make ξ itself — not merely its zero set — invariant under it. So π is exactly the constant whose presence centres the strip, and a coordinate that divides by π turns the strip into the polar range of a sphere with the critical line at the equator. The offset ½ is then the ratio (π/2)/π of a quarter turn to a half turn. That is a statement about coordinates, and §3 says so plainly: it is an exact change of variable and by itself places no zero anywhere.

The intended way around this limitation is not to argue from the coordinate change, but to treat it as the first of a set of **independent appearances of ½** that do not share a derivation: an amplitude Lagrangian whose minimum is ½; a phase balance where tan = 1; a unit-diameter circle; a half-integer offset in the zero-counting phase; an equal-envelope condition on the explicit formula. Each is exhibited as a short program (§4–§8). If ½ arises each time as a *consequence* rather than as an input, that is evidence about the structure of the problem; if some arise only because ½ was built into the setup, §4 says which. We then record what these appearances do and do not give toward RH (§9), what earlier versions of this work claimed that did not survive re-running (§10), and which hypotheses remain open (§11).

## 2. Scope, labels and the code

**Labels.** `ESTABLISHED` — a published theorem or definition, cited. `OURS` — this project's construction, framing or code. `THEORETICAL` — stated, not established; never presented as running fact. `CONJECTURE` / `OPEN` — named open items. `HEURISTIC` — convergent evidence, not deduction. The suffix `:CALCULATED` — this repository's code reproduces the result, and the derivation engine (`rhp/`, §12) can emit the steps; it is not a proof.

**What is code and what is not.** Mechanisms short enough to read are shown as numbered Listings, extracted from the source by `paper/make_paper.py` and executed cold (Convention B of the skill); sessions are shown as `>>>` transcripts (Convention D); long or established machinery is cited and pointed at (Convention A). Numbers in the text are taken from fresh runs and re-checked by `make_paper.py --check`. *This draft predates that generator: every number below was produced in the session that wrote it, by the scripts named in §12, and will be regenerated.*

**UNS.** The author's "Universal Native Space" is radial complex spherical-polar coordinates in which quantities are read on the decimal-log scale; ln 10 is the Jacobian between the decimal-log Native metric ds = d·log₁₀p and natural-log prime-address space, and the four d\* values (the Boundary 0.24631, d\*_RG, d\*_taut = Ω_ZS/ln 10, d\*_ln10 = d\*·ln 10 = Ω_ZS) are the roles of the chart — centre, equatorial plane, radial unit, metric bridge — and close on d\*·ln 10 = d\*_ln10 = d\*_taut·ln 10 = Ω_ZS (`Ainulindale/wiki/17`; with Ω_ZS = W(1) = 0.56714…; the d\*_spec ≈ 0.24600 reading is disputed and not used). For the pathway ζ(½+it) the chart used in §7 is: radius log₁₀|ζ|, angle arg ζ, with t kept as an explicit third coordinate (the spring), and, on the sphere, polar angle πσ (equator σ = ½) with azimuth t mod 2π and the winding kept separately. This is our reading of the author's definition and is `OURS`; the exact UNS form of ζ(s) itself is still to be stated by the author (§11).

## 3. The spine: factoring π out of ζ

### 3.1 π is the factor that centres the reflection `:CALCULATED`

Define Ξ(s) = π^(−s/2) Γ(s/2) ζ(s) (ξ without the polynomial prefactor) and N(s) = Γ(s/2) ζ(s) (π removed). `half_is_the_equator.py` (T1) evaluates |Ξ(s) − Ξ(1−s)| and |N(s) − N(1−s)| at 30 digits:

```text
                 s     |Xi(s)-Xi(1-s)|     |NoPi(s)-NoPi(1-s)|
      (0.3 + 7.1j)            3.27e-33                  0.0124
      (0.8 + 2.5j)                 0.0                    0.35
    (0.25 + 14.0j)             3.1e-36                 1.15e-5
     (0.62 + 0.4j)            2.96e-31                    1.54
```

With π^(−s/2) the symmetry holds to working precision; without it, it does not. This is the functional equation (`ESTABLISHED`, Riemann 1859); the observation that π is the *specific* factor responsible is the spine of this paper (`OURS` as framing).

### 3.2 The coordinate θ = πσ is an exact change of variable `:CALCULATED`

Let θ = πσ. The reflection σ ↦ 1−σ becomes θ ↦ π−θ, the equatorial reflection of a sphere with polar angle θ ∈ [0, π]; its fixed set is θ = π/2, i.e. σ = ½; σ = 0 and σ = 1 are the poles.

```text
    theta=   0deg   s = 0.0000   reflect-> 1.0000
    theta=  45deg   s = 0.2500   reflect-> 0.7500
    theta=  90deg   s = 0.5000   reflect-> 0.5000   <- FIXED
    theta= 135deg   s = 0.7500   reflect-> 0.2500
    theta= 180deg   s = 1.0000   reflect-> 0.0000
```

**Stated limit.** This map holds for *any* function with a reflection symmetry about ½. It is bookkeeping, not evidence about zeros. What it contributes is a geometry — strip = sphere's polar range, critical line = equator — in which the questions of §5–§9 can be asked. The phase offset ½ = (π/2)/π is the ratio of a quarter-turn to the period π; it is exact and, taken alone, carries no arithmetic content.

### 3.3 The J_N picture, and where it stops `THEORETICAL`

The map J_N(z) = i/z̄ is an anti-Möbius involution with J_N⁴ = id and the unit circle as its fixed boundary (`ESTABLISHED`, direct computation). Under w = 2s−1 followed by stereographic projection, Re(s) = ½ is the equatorial circle. If ζ transformed under J_N as the l = 1, m = 0 spherical harmonic Y₁⁰ = cos θ, Courant's nodal-domain theorem (1923) would confine its zeros to the nodal circle — the equator. **That identification (earlier called C1) is not established here and is the geometric form of the open step** (§9).

## 4. Where ½ emerges — a catalogue

For each appearance: *emergent* means ½ is the output of a computation whose input did not contain ½; *constructed* means ½ was imposed (for example by symmetrising); *composite* means it is a combination of other constants. The table is the paper's main empirical content; each row has a program (`rhp/`, §12) and a status.

| # | where ½ appears | emergent / constructed | check | status |
|---|---|---|---|---|
| 4.1 | fixed line of s ↦ 1−s̄ | constructed (the symmetry is s ↦ 1−s) | `half_is_the_equator` T2 | `ESTABLISHED` |
| 4.2 | π is the factor making ξ symmetric; (π/2)/π | emergent (from Γ, π and ζ's own functional equation) | T1 | `:CALCULATED` |
| 4.3 | unit-diameter circle: d = 1 ⇒ r = ½, C = π, A = π/4 (π-units: C = 1, A = ¼) | emergent | `GenerationalLineage` `circle`: `lines.descend("circle", 1.0, role="diameter")` | `:CALCULATED` |
| 4.4 | Amplitude-Lagrangian minimum: L(σ,E) = e^(−σE) + e^(−(1−σ)E), ∂L/∂σ = 0 ⇒ σ = ½, ∂²L/∂σ² = E²·L > 0 | constructed by the symmetric choice of L (see §5) | `noether_derivation.verify_derivation()` | `:CALCULATED` |
| 4.5 | in-phase/quadrature balance: L_(I\|O) = Q/I = e^((2σ−1)E) = 1 at σ = ½, phase arctan 1 = π/4 | follows from 4.4 | same | `:CALCULATED` |
| 4.6 | explicit formula: every zero's term has magnitude x^σ; one shared envelope iff all σ = ½ | conditional on RH, not evidence for it | T3 | `:CALCULATED` (arithmetic) |
| 4.7 | origin reached iff R = e: min\|R e^(iφ) + e e^(−ikφ)\| = \|R−e\| | emergent (trigonometric identity) | T4 | `:CALCULATED` |
| 4.8 | zero-counting phase: θ(γₙ)/π − (n−1) fluctuates about −½ | emergent | §6 | `:CALCULATED` |
| 4.9 | Γ(½)² = π; the Gaussian is its own Fourier transform | emergent, classical | — | `ESTABLISHED` |
| 4.10 | LC resonance XL = XC ⇔ tan = 1 ⇔ σ_self = ½ (the pancake-coil resonance) | emergent in the circuit; the identification with σ is `THEORETICAL` | §8 | `THEORETICAL` |

**Reading the table.** Rows 4.1, 4.4 and 4.6 are constructed or conditional: they cannot be offered as evidence that zeros lie on ½. Rows 4.2, 4.3, 4.7, 4.8 and 4.9 are emergent in the defined sense. The paper's claim is limited accordingly: ½ is over-determined by independent routes in the surrounding mathematics; **that is a statement about the problem's structure, not a proof that the zeros obey it.**

## 5. Zeta as a Lagrangian

### 5.1 The Amplitude Lagrangian `:CALCULATED`

The Euler product gives the amplitude of prime p at s = σ+it as |p^(−s)| = e^(−σE), E = ln p, and the functional equation pairs it with e^(−(1−σ)E). `noether_derivation.py` takes their sum as an *amplitude Lagrangian*

    L(σ, E) = e^(−σE) + e^(−(1−σ)E),

derives the Noether current J = −∂L/∂σ = E(e^(−σE) − e^(−(1−σ)E)), shows J = 0 ⇔ σ = ½, ∂²L/∂σ² = E²·L > 0 (a unique global minimum), and records the product J_forward × |J_backward| = e^(−E), independent of σ. Running it today: `ALL CHECKS PASS`.

**Stated limit.** L is symmetric under σ ↔ 1−σ *by construction*. Its minimum at ½ is the AM–GM equality case, a consequence of the symmetric choice of L, not of any property of ζ's zeros. The derivation shows the *form* of the currents is forced once L is chosen; it does not show that zeros are critical points of L. (§9, step 4.)

### 5.2 "ζ is a Lagrangian" `THEORETICAL`

The stronger reading: ζ(s) is a partition function Z(β) = Σ n^(−β) (each integer an energy level E = ln n); the Euler product is its factorisation over prime modes; the Hagedorn ceiling is the pole at s = 1, and the critical line lies past it. The partition-function identity is `ESTABLISHED` (Julia 1990; Bost–Connes 1995); the claim that the action whose on-shell value is ζ has the zeros as its stationary configurations is the author's reading and is **not shown** here. The project's operator Σ_RB = Σ_p p^(−σ)[R̂_p ⊗ ∂̂ + ∂̂† ⊗ B̂_p] is a candidate form (`rhp.h_rb_hat`, `THEORETICAL`); its σ is not a scalar in the full algebra (§10).

## 6. The spiral and the half-integer offset

On the critical line, ζ(½+it) traces a curve whose phase is governed by the Riemann–Siegel theta function θ(t) = arg Γ(¼ + it/2) − (t/2) ln π. Gram's relation reads N(T) = θ(T)/π + 1 + S(T). The earlier drafts of this repository asserted "θ(tₙ) ≈ nπ". Re-computing at 25 digits for the first 30 zeros (`mpmath.zetazero`, `siegeltheta`) gives instead:

```text
n   γₙ            θ(γₙ)/π      θ(γₙ)/π − (n−1)
1   14.13472514   -0.550253    -0.550253
2   21.02203964    0.570211    -0.429789
3   25.01085758    1.393598    -0.606402
4   30.42487613    2.671065    -0.328936
5   32.93506159    3.317271    -0.682729
…
28  95.87063423   26.197483    -0.802517
29  98.83119422   27.488727    -0.511273
30  101.317851    28.584204    -0.415796
```

So θ(γₙ) ≈ (n − 3/2)π with fluctuations of order S(T) — the "½" here is the midpoint convention of S(T) at a jump, exactly as Gram's relation requires. The mean spacing γₙ₊₁ − γₙ over the first 30 zeros is 0.981 of the Riemann–von Mangoldt value 2π/ln(γₙ/2π). **The earlier "≈ nπ" is corrected in §10.**

## 7. The flat pathway and its UNS map: what flattening overlays `:CALCULATED`

**Setting.** The familiar picture of ζ(½+it) (the "polar graph", t ≈ 0–50) is the image of the line σ = ½ in the ζ-plane. It has no cusp; it passes through the origin ten times below t = 50, once at each zero. Orientation: log-polar in the ζ-plane (radius log₁₀|ζ|, angle arg ζ) with t explicit; aperture: σ = ½, continuous in t. Code: `rhp/uns.py`, `rhp/flatten.py`, `rhp/render.py`; notebook `notebooks/01_flattening_map.ipynb`; figure `figures/flattening_map.svg`.

**7.1 The flat curve is a signed-radius polar curve (`ESTABLISHED` form).** Hardy's ζ(½+it) = e^{−iθ(t)} Z(t), Z real. The flat point has angle −θ(t) and signed radius Z(t), so a passage through the origin continues on the opposite side: the petals have no corners. Checks: max |Im(e^{iθ}ζ)| = 6·10⁻¹⁵ on 4001 points; Z agrees with `mpmath.siegelz` to 4·10⁻¹⁶.

**7.2 Direction at the hub (`:CALCULATED`).** The curve crosses the origin at the nth zero in the direction arg(iζ′(½+iγₙ)) mod π, and this equals −θ(γₙ) mod π (n = 1…10, agreement to 10⁻⁶; Z changes sign at every one). The fan of ten directions at the hub is the phase θ sampled at the zeros; the flat plot superimposes the ten passes and keeps none of it.

**7.3 Angular velocity (`:CALCULATED`).** From ζ = e^{−iθ}Z, Im(ζ̄ζ′) = −θ′(t)Z(t)², so dφ/dt = −θ′(t) independently of the radius, with θ′(t) = ½ Re ψ(¼+it/2) − ½ ln π ≈ ½ ln(t/2π). The winding is one-handed for t above about 2π (θ′ changes sign at t = 6.30 on the grid), and each zero is a pass with zero angular momentum. Verified at six values of t to 5 digits and against the numerical angular rate (max deviation 4·10⁻⁴ at spacing 0.0125).

**7.4 Axis crossings (`:CALCULATED`).** The curve meets the real axis where θ = nπ, the Gram points (Im ζ ≤ 5·10⁻¹⁵ at g₀ … g₅ = 17.846, 23.170, 27.670, 31.718, 35.467, 38.999).

**7.5 Artefact classification.** Writing the pathway as the spring (ρ = log₁₀|ζ| + 3, φ = arg ζ, z = t) and the flat plot as its end-on image with the radius exponentiated:

| artefact | mechanism | measured, 0 ≤ t ≤ 50 |
|---|---|---|
| hub | all zeros project to one point | 10 passes, C(10,2) = 45 pairwise crossings, all at the origin |
| false crossings | two different t land on one flat point | 21 away from the origin (66 self-intersections in all), with |Δt| from 3.4 to 47; a near-triple point at (0.443, −0.306) from t = 2.01, 42.99, 49.29/49.30 |
| petal sign | sign of Z folded into the angle as a π jump | sign flips at all 10 zeros |
| petal size | linear radius hides log structure | petal maxima span 0.34–3.66 linear, −0.47–+0.56 in log₁₀ |
| lost clock | t is not drawn | exact rate −θ′(t); −θ(t)/π is shown in Fig. 7, panel D |
| missed petal | Gram's law fails | none below t = 50; first at the Gram interval n = 125 ([280.80, 282.45) holds no zero), then n = 126, 133, 134 up to t = 300 (`ESTABLISHED`: the first Gram-law violation, reproduced) |
| mirror half | t < 0 is the conjugate | ζ(½−it) = conj ζ(½+it), difference 0 at four tested t |

The 21 and 45 are counts at 4001 samples (`:CALCULATED`, regression-tested). A viewport sweep of 500-sample windows covers every sample (clamp/flush/sweep assertions of the nes-viewport model), and the spring passes the canonical-pose framing check of the scad-spatial camera.

**7.6 The machinery (`:CALCULATED`; notebook 02, `figures/machinery.svg`).** For σ > 1, log ζ(s) = Σ_p Σ_k p^{−ks}/k is a sum of prime wheels (radius p^{−kσ}/k, frequency k·ln p), so in the UNS log plane the pathway is the sum of the wheels, where in the flat plane they multiply. The maximum error over t ≤ 60 after K wheels is 5.2·10⁻¹, 2.2·10⁻¹, 7.0·10⁻², 1.1·10⁻² for K = 1, 3, 10, 100 at σ = 1.5, and 1.7, 1.2, 0.79, 0.43 at σ = 1.1: the number of wheels needed diverges as σ → 1, and on σ = ½ the series does not converge. There the curve is drawn by the approximate functional equation (`ESTABLISHED`) ζ ≈ A + χB, with A = Σ_{n≤N} n^{−s} (forward wheels), B = Σ_{n≤N} n^{s−1}, N = ⌊√(t/2π)⌋. On σ = ½, |χ| = 1 and B = Ā, so ζ = A + e^{−2iθ}Ā and the two halves have equal modulus (max ||A| − |χB|| = 3·10⁻¹⁴ on 25 points, t ∈ [20, 120]); off the line one half dominates, with median |χB|/|A| = 1.17 at σ = 0.4 and 0.85 at σ = 0.6, against (t/2π)^{½−σ} = 1.27 and 0.79 at the median t. This is the pathway's own two-hemisphere structure: the equator is where the forward wheel sum and its functional-equation mirror balance. The primes appear in the other direction through the explicit formula: the zeros are the spectrum of ψ(x) − x (rms error on [2, 60] with the first 10, 100, 3066 zeros: 0.74, 0.37, 0.095), and a Gaussian factor e^{−(γw)²/2} on each zero gives a continuous sieve that sharpens to the discrete jumps as w → 0.

**7.7 The flow of ξ′/ξ and the accounting (`:CALCULATED`; `figures/flow.svg`).** Let Φ = log ξ and L = ξ′/ξ; the velocity (u, v) = (Re L, −Im L) has the zeros of ξ as sources with vortex. From ξ(s) = ξ(1−s) and ξ(s̄) = conj ξ(s): L(1−s̄) = −conj L(s), so u is odd and v is even under σ ↦ 1−σ (residuals 3·10⁻¹⁶), and |L| is even: the hemispheres carry equal energy density, move oppositely in σ and identically in t. On a 60 × 80 grid (σ ∈ [0.05, 0.95] off the line, t ∈ [0.5, 60]) u·(σ−½) > 0 at every point (0 violations in 4800; minimum u/(σ−½) = 0.046): outflow from the line in both hemispheres. The sign condition is equivalent to RH (the statement recalled as Lagarias's positivity equivalence; to be verified by title page), so this is a picture of the equivalence, not evidence beyond the known verification. The argument principle is Gauss's law for this flow: the number of zeros in 0 < σ < 1, 0 < t < T equals θ(T)/π + 1 + S(T), with S followed continuously along 2 → 2+iT → ½+iT, and the number on the line is the number of sign changes of Z. At 30 heights T ≤ 3256 the two agree with |gap| ≤ 5·10⁻¹³ (floating point; Turing's method with interval arithmetic makes this rigorous, as in Platt–Trudgian).

**7.8 The long pathway (`figures/long_pathway.svg`).** 3065 zeros in t ∈ [20, 3600] (3066 with γ₁ = 14.13) are 1533 full turns of the angle −θ. Riemann–Siegel with the first correction (max |Z − Z_mpmath| ≈ 6·10⁻⁴ at t ≤ 3600) draws the curve; zeros are polished with mpmath (cached in `data/zeros_1_to_3600.csv`). The hub fan of §7.2 is S(γₙ): the direction at the origin of the nth pass is π/2 + π·S(γₙ) (mod π); over the 3065 zeros S(γₙ) mod 1 has mean +0.004 and sd 0.254, and S(t) between zeros has mean 0.000, sd 0.195, max |S| = 0.54. The envelope of the flat overlay is cardioid-like; it is the envelope of a statistical overlay of 1533 loops, not the curve, and no fit was made.

**7.9 The θ-clock and the real axis (`:CALCULATED`; notebook 03, `figures/clock_echo_harness.svg`).** In ζ ≈ A + χB the n = 1 pair is 1 + e^{−2iθ}, a circle of radius 1 centred at 1 through the origin, reached when θ ≡ π/2 (mod π); its hits E_k, θ(E_k) = (k − 3/2)π, carry the Weyl density. The 3066 zeros below 3600 lie within 0.81 mean spacings of the carrier (mean 0.22), the displacement being S(γ_k) (θ(γ_k) = π(k − 3/2 − S)). The pathway meets the real axis only at the Gram points θ = nπ, where ζ = (−1)ⁿZ(g_n): the crossing is negative exactly when Gram's law fails. Below 3600: 3066 Gram points, 197 negative crossings (6.4%; by height blocks 2.7%, 5.7%, 7.4%, 7.6%, 7.6%, 7.4%), the first at n = 126 (t = 282.455); the 393 Gram intervals with ≠ 1 zero are all adjacent to a failing point.

**7.10 The echo spectrum (`:CALCULATED`).** With a Gaussian window (Γ = 300) the zeros' Fourier sum F(t) = Σ_{γ>0} e^{−(γ/Γ)²/2} cos γt equals the explicit-formula prediction to 8·10⁻¹¹ (max |F| = 42): the pole term cosh(t/2) (the path of least primes), a smooth ψ-term, and negative Gaussian bumps of area ½Λ(n)/√n at t = ln n (measured against predicted heights: −29.33/−29.33 at ln 2, −43.05/−43.07 at ln 5, −16.62/−16.64 at ln 49). The trace of cos(Ht) for a hypothetical H with spectrum {γ} therefore has the primes as fixed echoes.

**7.11 An operator harness (`:CALCULATED`; `rhp/harness.py`).** A candidate spectrum is scored against the zeros on four columns: deviation from the smooth count, variance of unfolded spacings, number variance Σ²(5), and mean distance to the actual zeros in spacings. Reference zeros: 0.81, 0.149, 0.34, 0. Controls: the carrier 0, 0, 0, 0.22; a ladder (n+½)ħ (the spectrum of ValaQuenta's H_NN as coded, and of xp on a periodic u-box) 202, 0.03, 1.39, 137; Poisson 31, 0.97, 4.7, 13; GUE unfolded onto the zeros' density 2.5, 0.18, 0.49, 0.64. No candidate operator is claimed: a real one (a boundary condition on xp, or a Σ_RB construction) must beat the carrier and GUE on location while keeping the statistics. The Berry–Keating Hamiltonian is not reproduced here; xp on the whole u-line has continuous spectrum, and the boundary condition that selects the zeros is the open content.

**7.12 Margins, balance and PRIMEGAME (`:CALCULATED`).** (i) Near-double zeros: the smallest relative extremum height of Z in a gap below 3600 is 0.016, at n = 1496, 1497 (γ = 1977.1739, 1977.2714; gap 0.089 spacings; |Z| extremum 0.0221, confirmed by mpmath inside the gap). (ii) On the truncated halves of the approximate functional equation, M = [ln|χ| + ln|B| − ln|A|]/(½ − σ) = ln(t/2π) + c with |c| ≤ 2 ln N would make a sign change impossible, but ln(t/2π) − 2 ln N = ln(t/2πN²) vanishes at t = 2πN²: the crude bound has zero slack where a new wheel switches on (56.5, 100.5, 157.1, …); 664 off-line sign changes of M, on 110 smooth branches, occur for t ∈ [20, 400] and σ ∈ (0.02, 0.49), with |ζ| at the branch midpoints 0.4–1.5 (they are not zeros). Their distance to the nearest zero (0.198 local spacings) is smaller than for random t (0.253, null sd 0.023): a weak (2.4σ) association whose mechanism is unidentified. For the exact functions the balance is an identity (ζ(s) = χ(s)ζ(1−s)), so only the truncated balance is informative. (iii) Conway's PRIMEGAME (14 fractions) is a register machine on the exponents of 2, 3, 5, 7 (counters) and 11, 13, 17, 19, 23, 29 (six 0/1 flags); it tests every smaller divisor, i.e. whether a rotation by d on a circle of n positions closes early (period n/gcd(n,d)); the p-th prime costs about 2.25·p^{2.89} steps (steps for 2, 3, 5, 7, 11: 19, 69, 280, 707, 2363). A published figure set (281, 710, 2375) equals these steps plus the number of firings of 15/14 (1, 3, 12) exactly, a counting convention.

**7.13 Anomalies chased.** Ten anomalies were chased to a verdict (notebook 03 §9): two METHOD errors found and fixed (an under-resolved window in a first echo check; a Lehmer search window that overlapped neighbouring lobes), one counting convention, one video imprecision, two expected effects, one finite-height effect, one correction of an earlier suggestion of the assistant (a "nodal line of the modulus balance" is vacuous for the exact functions), one expected control failure (H_NN), and one UNRESOLVED association (7.12 ii).

**7.14 What is not claimed.** Nothing in §7 places a zero or tests RH beyond the verification already known; it maps what the flat picture overlays. That the zeros are sources whose outflow is the factor/prime exchange of the author's Two Trees is `THEORETICAL`: the conserved quantity that carries it (J_fwd + J_bwd = 0) is a different one and is not computed here. The sphere chart as a picture is `OPEN`; the double-zero margin map is §7.12 (i); no operator candidate is claimed (§7.11). The limaçon sifting fit of earlier drafts (loops classed by |b/a|; the cardioid at |b/a| = 1 a classification threshold, not the curve) has no code in this repository and is not used; the earlier hypothesis that a cusp is a flattening artefact is withdrawn, the pathway having none.

## 8. The prime wavelength, and the pancake coil

Two wavelengths, one per ordering of the same set of primes:

* **Prime wavelength** λ_prime(n) = γₙ₊₁ − γₙ, the spacing of consecutive zeros. `ESTABLISHED` (Riemann–von Mangoldt; Montgomery 1973; Odlyzko 1987). The name is a relabelling of a heavily studied object.
* **Factor wavelength** λ_factor(p) = p, the exact period with which a prime, acting as a factor, recurs in the integers. `ESTABLISHED` mathematics (the sieve; Mertens; Legendre's φ(x)/x); nearest modern prior art Torquato–Zhang–de Courcy-Ireland (2018) and Dyson (2009). Original only as a framing: the sieve read as a wave equation, the Fermat-ordinal conjugate of λ_prime.

**Conjecture (author, 2026-07-23; `OPEN`).** λ_prime and λ_factor are orthogonal in the E ⟂ B sense, as the number-theoretic instance of the two-octonion orthogonality of the sedenions (L(eₖ) skew-symmetric with eigenvalues exactly {+i, −i}, k = 1…15; `ESTABLISHED`, checked). **No orthogonality relation between the two wavelength quantities has been derived or computed**; the algebraic parallel alone is not evidence. The test is §11.

**The coil.** The prime wavelength was not an expected outcome: it appeared while analysing the resonance of a pancake coil (the POE design), where XL = XC ⇔ tan(phase) = 1 ⇔ σ_self = ½. The circuit fact is elementary (`ESTABLISHED`); that the coil's tuning wavelength corresponds to the zero spacing is `THEORETICAL` and is recorded here as the provenance of the idea, not as a result.

## 9. What is proved, what is imported, and where the gap is

| step | content | status |
|---|---|---|
| 1 | ξ(s) = ξ(1−s) | `ESTABLISHED` (Riemann 1859); π's role `:CALCULATED` (§3.1) |
| 2 | the reflection has fixed line Re(s) = ½; the strip is the sphere's polar range under θ = πσ | `ESTABLISHED` / exact change of variable (§3.2) |
| 3 | an amplitude Lagrangian built from the Euler product has Noether currents that balance only at σ = ½ | `:CALCULATED` for the stated L (§5.1); depends on the choice of L |
| 4 | the zeros are stationary points of that Lagrangian / stable equilibria of H = xp | **`OPEN`** — the Berry–Keating self-adjointness conjecture (Berry–Keating 1999; Hilbert–Pólya) |
| 4′ | ζ transforms as Y₁⁰ under J_N (geometric form of step 4) | **`OPEN`** (earlier "C1") |
| 5 | Re(ρ) = ½ for all non-trivial zeros | follows from 1–4; not established because 4 is not |

**An equivalence that is a restatement.** Earlier drafts defined a zero to be *stable* when J(σ,E) = 0 and then proved "RH ⇔ every zero is stable". Since J = 0 ⇔ σ = ½ for E ≠ 0, "stable" is a re-statement of "σ = ½"; the equivalence is true and carries no content beyond RH. The real content of step 4 is an argument that zeros *must* satisfy the balance, which this paper does not supply.

## 10. Corrections kept in the record

1. **"θ(tₙ) ≈ nπ"** (earlier drafts): the computation of §6 gives (n − 3/2)π with S(T) fluctuation.
2. **"Forced σ = ½ exactly"** (the semantic engine's output `σ = 0.500000`): `forced_sigma` solves the linear balance E(1−2σ) = 0, whose solution is ½; the output demonstrates the algebra of §5, not a measurement on zeros. (An earlier softmax iteration was replaced by a one-step Newton solve; it had failed for large E.)
3. **BAO "spectral residue = mass gap exactly"**: disproven as stated — the residue is the same Ω − d*·ln 10 subtraction renamed; Δ = 1/(1000√2) is a 0.035 % coincidence; Δ > 0 and Δ ≈ Ω − d*·ln 10 stand.
4. **σ in ∅_RB is not a scalar** (the Oblique Gear): σ is a tilt (Re) plus an axis (Im), two 8-vectors; σ = ½ is a hyperbolic saddle (λ ≈ 6). A scalar σ is its projection.
5. **d\* has four recorded values** (the Boundary 0.24631; d\*_RG, dimension 8, numerically open; d\*_taut = Ω/ln 10; d\*_ln10 = Ω_ZS/ln 10) and a disputed fifth, d\*_spec = 0.24600. They are not interchangeable.
6. **"We prove the Riemann Hypothesis"** (an earlier README abstract) is not supported by that document's own body, which names the mode identification as open. This paper does not make that claim.
7. **A tension recorded, not smoothed:** the Hagedorn-ceiling reading (the critical line lies past the pole at s = 1) and the Lambert-W reading of the constants are not reconciled here (`OPEN`).

## 11. Open items and the tests they need

| item | what would settle it |
|---|---|
| step 4 / Berry–Keating self-adjointness | a self-adjoint operator whose spectrum is the zeros' imaginary parts, with a boundary condition that makes xp self-adjoint (the obstruction of Berry–Keating 1999 and the literature that followed) |
| a candidate operator passing all four harness columns (§7.11); σ ≠ ½ sweep as a picture; the sphere chart; the association of balance branches with zeros (§7.12 ii) | ζ(σ+it) at σ = 0.3 … 0.7 in the same three charts: finite log-radius dips replacing the spikes, and the change of winding across σ = ½ |
| λ_prime ⟂ λ_factor (§8) | a defined inner product on the two quantities and a computed value, on the first N zeros and primes ≤ P, with a null model |
| ζ as a Lagrangian (§5.2) | an action whose Euler–Lagrange stationary set is shown to equal the zero set |
| exact UNS form of ζ(s) (§2) | the author's explicit definition; §7 uses the log-polar/spring reading stated in §2 |
| the four d\* values (§10.5) | derivation or numerical resolution of d\*_RG |

## 12. Reproducing this paper

Everything numeric above was produced by: `half_is_the_equator.py` (§3, §4; 30-digit mpmath); `engines/noether_derivation.py` (§5); an `mpmath` script for §6; `GenerationalLineage` `circle` (§4.3); the project's limaçon fit (§7; to be re-run). The rest are to be vendored into `rhp/` with provenance headers, tested, wrapped in notebooks, and the paper regenerated by `paper/make_paper.py` (see `PLAN.md`).

## References (to be verified by title page under the Literature Download Protocol)

Riemann, B. (1859), *Über die Anzahl der Primzahlen unter einer gegebenen Grösse.* · Noether, E. (1918), *Invariante Variationsprobleme.* · Courant, R. & Hilbert, D., *Methods of Mathematical Physics*, vol. 1 (nodal-domain theorem; Courant 1923). · Hilbert–Pólya conjecture (as discussed in Montgomery 1973). · Berry, M. V. & Keating, J. P. (1999), *H = xp and the Riemann zeros*, in *Supersymmetry and Trace Formulae.* · Montgomery, H. L. (1973), *The pair correlation of zeros of the zeta function.* · Odlyzko, A. M. (1987), *On the distribution of spacings between zeros of the zeta function.* · Titchmarsh, E. C., *The Theory of the Riemann Zeta-Function* (2nd ed., Heath-Brown). · Edwards, H. M., *Riemann's Zeta Function.* · Julia, B. (1990), *Statistical theory of numbers.* · Bost, J.-B. & Connes, A. (1995), *Hecke algebras, type III factors and phase transitions…* · Torquato, S., Zhang, G., de Courcy-Ireland, M. (2018), *Uncovering multiscale order in the prime numbers via scattering* (J. Stat. Mech.). · Dyson, F. (2009), *Birds and Frogs*, Notices AMS 56(2). · Abrikosov, A. A. (1957), *On the magnetic properties of superconductors of the second group.* · Selberg, A. (1956); Deligne, P. (1974) (context). · Gourdon, X. (2004), *The 10¹³ first zeros of the Riemann zeta function* (verification of RH to that height).

## Appendix A — ContextPlease: what it is, and what this draft used

ContextPlease is the author's private archive of working context for the assistant: primers, scratchpads, memory files and shortcuts, kept so any model can come up to speed with historical access. This draft used: `~/.clauderc_canonical_maths` (PRIME WAVELENGTH, LIMACON_THRESHOLD, constants); `~/.clauderc_user_provenance` §1.18; `scratchpad/2026-08-13_apex_path/half_is_the_equator.py`; `Ainulindale/wiki/73_why_the_half_line.md`; `ValaQuenta/wiki/singularity_null.md`; the `cs-paper-code-conventions` skill.
