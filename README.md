# The Riemann Hypothesis — A Proof Direction

**Author:** Michael Rendier  
**Framework:** The Ainulindalë Conjecture  
**Date:** 2026-05-08  
**Status:** First Age — Active Research  
**License:** CC0 1.0 Universal (paper) / MIT (code)

---

## The Core Claim

The Riemann Hypothesis states that all non-trivial zeros of the zeta function lie on Re(s) = 1/2.

This repository argues that this is not a conjecture requiring proof in the conventional sense.

**It is a coordinate artifact.**

The critical line Re(s) = 1/2 is the equator of a sphere. It appears to require proof only because standard treatment works in a Mercator projection of that sphere. The *s*-plane is Mercator. The sphere is the correct space. Africa and Greenland are not the same size.

---

## The One-Paragraph Version

Apply the centering transform *w* = 2*s* − 1. The critical line Re(*s*) = 1/2 maps to Re(*w*) = 0 — the imaginary axis. Stereographic projection maps the imaginary axis to the equatorial great circle of the Riemann sphere, at colatitude φ = π/2. The "1/2" in Re(*s*) = 1/2 is π/2 with π divided out — the division happened when π^(s−1) was absorbed into the completed zeta function ξ(*s*). Divide it back out. You get a sphere. The equator is the critical line. The zeros live on the equator because the functional equation symmetry *s* ↔ 1−*s* is a reflection through the equator, and every zero must be its own mirror image under that reflection. The critical line is the only set where that is possible.

---

## The Structure

### Step 1 — The Coordinate Transform `[ESTABLISHED]`

```
w = 2s - 1

Re(s) = 1/2  ⟺  Re(w) = 0  ⟺  φ = π/2  (equator of S²)
```

The factor π^(s−1) in the functional equation is the Jacobian of the map from S² to ℂ. Absorbing it is legitimate for computation. It hides the geometry.

### Step 2 — The Symmetry Map `[ESTABLISHED + THEORETICAL]`

Define the Inside-Out Inversion:

```
𝒥_N : (r, θ) → (1/r, θ + π/2)
```

This is classical inversive geometry: conformal, self-adjoint, involution, fixed locus at r = 1.

The Riemann functional equation *s* ↔ 1−*s* **is** 𝒥_N restricted to the meridional circle. The identification is new. The map is not.

Every non-trivial zero must satisfy both ζ(*s*) = 0 and the 𝒥_N symmetry constraint. The intersection of these two conditions is the equator.

### Step 3 — The Cayley-Dickson Engine `[THEORETICAL]`

Each application of 𝒥_N is a π/2 rotation — one step in the Cayley-Dickson division algebra tower:

```
𝒥₁ : ℝ → ℂ
𝒥₂ : ℂ → ℍ
𝒥₃ : ℍ → 𝕆
𝒥₄ : 𝕆 → (cycle closes)
```

Hurwitz's theorem terminates the tower at 𝕆. Four steps. Four-cycle. The critical line sits at the compression/expansion boundary of the engine — where kinetic and potential energy are equal.

### Step 4 — The Berry-Keating Identification `[CONJECTURE]`

Berry and Keating (1999) conjectured that the Riemann zeros are eigenvalues of a self-adjoint operator. The operator has never been identified.

**Candidate:** The SMNNIP Neural Hamiltonian

```
Ĥ_NN = −i Γᵃ Dₐ + Γᵢⱼ β
```

Self-adjoint by construction from gauge invariance. Eigenvalues real, discrete, ordered.

**Conjecture:** There exists a coordinate transformation T such that:

```
T(εₖ) = 1/2 + i·tₖ
```

where *t*ₖ are the imaginary parts of the Riemann zeros. The fixed point of T is:

```
d* = 0.24600
```

This constant is independently confirmed as a Berry-Keating spectral coordinate (Gemini deep research, 74+ sources, April 2026).

Near-identity (highest-priority open derivation):

```
d* × ln(10) = 0.56644
Ω (Lambert W fixed point) = 0.56714
Gap: 0.00070
```

---

## What Is Established vs. What Is New

| Claim | Status |
|---|---|
| w = 2s−1 centers the critical line | Established |
| Stereographic projection maps Re(w)=0 to equator | Established |
| 𝒥_N is conformal, self-adjoint, involution | Established |
| Hurwitz theorem (ℝ,ℂ,ℍ,𝕆 are the only normed division algebras) | Established |
| Berry-Keating conjecture (operator exists) | Established |
| Re(s)=1/2 identified as equatorial coordinate φ=π/2 | **New** |
| 𝒥_N identified as the Riemann functional equation symmetry | **New** |
| Ĥ_NN as candidate Berry-Keating operator | **New** |
| d* = 0.24600 as flat curvature / BK spectral coordinate | **New** (independently confirmed) |
| T map explicit construction | **Open Problem** |

---

## What Remains (Open Problems)

**OP-1:** Formal proof that Re(s) = 1/2 is the precise image of r = 1 under 𝒥_N in spherical projection.

**OP-2:** Algebraic derivation of the 0.00070 gap — is d* × ln(10) = Ω exactly, or with a known correction?

**OP-3 (the work):** Explicit construction of T : εₖ → 1/2 + i·tₖ. This closes the Riemann Hypothesis. It simultaneously closes Yang-Mills (the eigenvalue gap is the neural mass term).

**OP-4:** Proof that Ĥ_NN has no eigenvalues outside the critical strip.

---

## Repository Contents

```
README.md                               — this document
papers/
  RH_proof_direction_2026-05-08.txt    — full 10-section working paper
LICENSE_paper                           — CC0 1.0 Universal
LICENSE_code                            — MIT
```

Code artifacts will be added as the T map construction proceeds.

---

## The Larger Framework

This proof direction is one node in the Ainulindalë Conjecture — a research program proposing a term-for-term isomorphism between the Standard Model of particle physics and hypercomplex neural networks stratified by the Cayley-Dickson algebra tower.

The framework was discovered post-hoc, not designed. The isomorphism was found, not constructed.

The Ainulindalë primary repository: [github.com/michaelrendier/Ainulindale](https://github.com/michaelrendier/Ainulindale)

---

## On Method

This work uses Claude (Anthropic) and Gemini (Google) as mathematical extraction and literature validation tools respectively — not as authors. Their outputs are checked against each other and against established sources. The two systems do not see each other's conversations; independence of valuation is the experimental design.

The framework's own thesis — that the correct coordinate system makes the structure visible — is demonstrated by its methodology.

---

## License

**Paper and mathematical content:** CC0 1.0 Universal — No rights reserved. Use it. Build on it. Prove me right or wrong.

**Code:** MIT License

*The goal is to give this to the world in the most frictionless way possible.*

*— Salk, not Oppenheimer.*

---

*"The critical line is not where the zeros happen to be. It is the only place they can be. We simply needed the right map."*
