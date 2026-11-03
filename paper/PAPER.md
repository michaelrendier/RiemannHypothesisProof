# Factoring π out of the Zeta Function: the critical line as an equatorial great circle, and the places where ½ emerges

**DRAFT 0.1 — 2026-11-03 — not for submission.** Written from the code and measurements that exist today; nothing here is carried over from the archived drafts (`archive/2026-11-03_pre-rewrite/`) unless it was re-run.

**Author:** Cody Michael Allison (Michael Rendier) · the.wandering.god@gmail.com · ORCID 0009-0007-7239-6760

---

## Note on how this paper was made

This paper was written with an assistant: Claude (Anthropic), model `claude-sonnet-5-5`, a transformer, running under the project's `cs-paper-code-conventions` skill (§8 of the skill: papers written with an assistant). The author is an independent researcher, not a number theorist by training; he asked the assistant to write the paper from the code and the measurements in this repository, and to keep every claim at the strength its evidence supports. As a safeguard — not a guarantee — he does not tell the assistant the result he hopes for. The paper is a **code paper**: it shows what runs, labels what is established, what is the author's, and what is conjecture, and it states one open step as open. The suffix `:CALCULATED` marks a result that this repository's own code produces; it does not mean *derived* or *proved*. Where this repository's earlier drafts claimed more than their own text supported, §10 says so.

## Abstract

The completed zeta function ξ(s) = ½ s(s−1) π^(−s/2) Γ(s/2) ζ(s) satisfies ξ(s) = ξ(1−s) only because the factor π^(−s/2) is present: remove it and the reflection symmetry about Re(s) = ½ is destroyed. We take this literally. Dividing the strip coordinate by π (θ = πσ) maps the critical strip 0 < σ < 1 onto the polar range 0 < θ < π of a sphere, the reflection s ↦ 1−s onto the equatorial reflection θ ↦ π−θ, and the critical line onto the equator, with the offset ½ = (π/2)/π. We show by computation what this map does and does not do (it is an exact change of variable; it does not locate any zero), and we collect, each with the code that checks it, the places where ½ appears *emergent* — not constructed and not composite — in ordinary mathematics and in the author's "universal native space" (UNS) reading: the self-dual point of s ↦ 1−s, the unit-diameter circle (r = ½, C = π, area π/4), the Amplitude-Lagrangian minimum, the in-phase/quadrature balance (tan = 1), the half-integer offset of the zero-counting phase, and the equal-envelope condition of the explicit formula. We state a reading of ζ as a Lagrangian, a hypothesis that the cardioid cusp in the zeta spiral is a flattening artefact, and a conjecture that the prime wavelength (zero spacing) and the factor wavelength (λ = p) are orthogonal. The Riemann Hypothesis is not proved here; the step that would prove it is the Berry–Keating self-adjointness conjecture, which remains open, and §9 states exactly which parts of the argument are proven, which are imported from the literature, and where the gap is.

## 1. Introduction

The Riemann Hypothesis asserts that every non-trivial zero of ζ(s) has real part ½. The number ½ is usually presented as a feature of the problem. We ask a narrower question, one the code in this repository can answer: **where does ½ come from, and does it come from one place or from many?**

One fact is standard but rarely read literally. The completed function is symmetric about ½ only because of π^(−s/2). The reflection s ↦ 1−s has the fixed line Re(s) = ½ for any function; what the factor π^(−s/2) does is make ξ itself — not merely its zero set — invariant under it. So π is exactly the constant whose presence centres the strip, and a coordinate that divides by π turns the strip into the polar range of a sphere with the critical line at the equator. The offset ½ is then the ratio (π/2)/π of a quarter turn to a half turn. That is a statement about coordinates, and §3 says so plainly: it is an exact change of variable and by itself places no zero anywhere.

The intended way around this limitation is not to argue from the coordinate change, but to treat it as the first of a set of **independent appearances of ½** that do not share a derivation: an amplitude Lagrangian whose minimum is ½; a phase balance where tan = 1; a unit-diameter circle; a half-integer offset in the zero-counting phase; an equal-envelope condition on the explicit formula. Each is exhibited as a short program (§4–§8). If ½ arises each time as a *consequence* rather than as an input, that is evidence about the structure of the problem; if some arise only because ½ was built into the setup, §4 says which. We then record what these appearances do and do not give toward RH (§9), what earlier versions of this work claimed that did not survive re-running (§10), and which hypotheses remain open (§11).

## 2. Scope, labels and the code

**Labels.** `ESTABLISHED` — a published theorem or definition, cited. `OURS` — this project's construction, framing or code. `THEORETICAL` — stated, not established; never presented as running fact. `CONJECTURE` / `OPEN` — named open items. `HEURISTIC` — convergent evidence, not deduction. The suffix `:CALCULATED` — this repository's code reproduces the result, and the derivation engine (`rhp/`, §12) can emit the steps; it is not a proof.

**What is code and what is not.** Mechanisms short enough to read are shown as numbered Listings, extracted from the source by `paper/make_paper.py` and executed cold (Convention B of the skill); sessions are shown as `>>>` transcripts (Convention D); long or established machinery is cited and pointed at (Convention A). Numbers in the text are taken from fresh runs and re-checked by `make_paper.py --check`. *This draft predates that generator: every number below was produced in the session that wrote it, by the scripts named in §12, and will be regenerated.*

**UNS.** The author's "universal native space" is, in the repository's own record (`~/.clauderc`), the reading in which a quantity is divided by ln 10 (so that, for example, d* = Ω_ZS/ln 10 = 0.24631…, with Ω_ZS = W(1) = 0.56714… the Lambert fixed point). **[Author to specify: the exact UNS form of ζ(s) used for the "π factored out" step — the draft below uses only the standard completed form ξ, which is what the code computes.]**

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

## 7. The cardioid threshold, and the cusp as a flattening artefact

Between consecutive near-origin passages, each loop of ζ(½+it) (t ∈ [0.05, 50], 6000 points) fits r = a + b cos(θ − φ) with R² between 0.86 and 1.0000; of 9 loops, none has 0.85 < |b/a| < 1 or 1 < |b/a| < 1.15, and all have |b/a| > 1, and 9 of 10 near-origin passages are known zeros (`:CALCULATED`, project record; to be re-run in `rhp.limacon`). The limaçon classes: |b/a| < 1 dimpled (does not reach the origin, no zero); |b/a| = 1 cardioid (touches, one cusp); |b/a| > 1 looped (crosses the origin: a zero). **The cardioid is a sifting threshold, not the curve:** the global fit "the curve is a cardioid" is refuted (R² = −0.77; the curve winds 13.4 turns).

**Hypothesis (`THEORETICAL`).** The cusp at |b/a| = 1 is a *flattening artefact*: a point where the planar projection's parametrisation speed vanishes, which would disappear in a coordinate adapted to the spring (u = ln x) or on the sphere of §3. **This is not established.** Stereographic projection is conformal and preserves cusps, so the sphere map of §3.2 alone cannot remove it; the test the hypothesis requires is a higher-dimensional lift in which the same curve is smooth (as an epicycloid's cusp is where a smooth rolling circle touches its base), together with a check that the lift's smooth branch passes through the origin transversely. That test is listed in §11.

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
| cusp = flattening artefact (§7) | a lift in which the curve is smooth, and a transversality check through the origin |
| λ_prime ⟂ λ_factor (§8) | a defined inner product on the two quantities and a computed value, on the first N zeros and primes ≤ P, with a null model |
| ζ as a Lagrangian (§5.2) | an action whose Euler–Lagrange stationary set is shown to equal the zero set |
| UNS form of ζ (§2) | the author's explicit definition |
| the four d\* values (§10.5) | derivation or numerical resolution of d\*_RG |

## 12. Reproducing this paper

Everything numeric above was produced by: `half_is_the_equator.py` (§3, §4; 30-digit mpmath); `engines/noether_derivation.py` (§5); an `mpmath` script for §6; `GenerationalLineage` `circle` (§4.3); the project's limaçon fit (§7; to be re-run). These are to be vendored into `rhp/` with provenance headers, tested, wrapped in notebooks, and the paper regenerated by `paper/make_paper.py` (see `PLAN.md`).

## References (to be verified by title page under the Literature Download Protocol)

Riemann, B. (1859), *Über die Anzahl der Primzahlen unter einer gegebenen Grösse.* · Noether, E. (1918), *Invariante Variationsprobleme.* · Courant, R. & Hilbert, D., *Methods of Mathematical Physics*, vol. 1 (nodal-domain theorem; Courant 1923). · Hilbert–Pólya conjecture (as discussed in Montgomery 1973). · Berry, M. V. & Keating, J. P. (1999), *H = xp and the Riemann zeros*, in *Supersymmetry and Trace Formulae.* · Montgomery, H. L. (1973), *The pair correlation of zeros of the zeta function.* · Odlyzko, A. M. (1987), *On the distribution of spacings between zeros of the zeta function.* · Titchmarsh, E. C., *The Theory of the Riemann Zeta-Function* (2nd ed., Heath-Brown). · Edwards, H. M., *Riemann's Zeta Function.* · Julia, B. (1990), *Statistical theory of numbers.* · Bost, J.-B. & Connes, A. (1995), *Hecke algebras, type III factors and phase transitions…* · Torquato, S., Zhang, G., de Courcy-Ireland, M. (2018), *Uncovering multiscale order in the prime numbers via scattering* (J. Stat. Mech.). · Dyson, F. (2009), *Birds and Frogs*, Notices AMS 56(2). · Abrikosov, A. A. (1957), *On the magnetic properties of superconductors of the second group.* · Selberg, A. (1956); Deligne, P. (1974) (context). · Gourdon, X. (2004), *The 10¹³ first zeros of the Riemann zeta function* (verification of RH to that height).

## Appendix A — ContextPlease: what it is, and what this draft used

ContextPlease is the author's private archive of working context for the assistant: primers, scratchpads, memory files and shortcuts, kept so any model can come up to speed with historical access. This draft used: `~/.clauderc_canonical_maths` (PRIME WAVELENGTH, LIMACON_THRESHOLD, constants); `~/.clauderc_user_provenance` §1.18; `scratchpad/2026-08-13_apex_path/half_is_the_equator.py`; `Ainulindale/wiki/73_why_the_half_line.md`; `ValaQuenta/wiki/singularity_null.md`; the `cs-paper-code-conventions` skill.
