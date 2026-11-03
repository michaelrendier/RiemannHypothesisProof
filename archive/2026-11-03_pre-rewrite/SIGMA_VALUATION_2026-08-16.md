# Sigma Valuation — RiemannHypothesisProof, 2026-08-16

**Assessor:** Claude Opus 5 (Anthropic)
**Scope:** PAPER.md (1425 lines), README.md, `engines/noether_derivation.py`, notebooks 01–08
**Predecessor:** `SIGMA_VALUATION.md` (Claude Sonnet 4.6, 2026-05-10). This pass extends
rather than replaces it, and reaches a **lower** σ on the central argument for a reason
that valuation did not test.

| σ | Meaning |
|---|---------|
| ∞ | Proven theorem — no uncertainty |
| 5+ | Computationally verified beyond reasonable doubt |
| 3–5 | Strong evidence; formally open |
| 2–3 | Plausible, well-motivated; requires formal work |
| 1–2 | Suggestive but underspecified; gap present |
| < 1 | Speculative; connection to the main claim unclear |

---

## Headline

**The central Noether argument does not currently support its conclusion, and the reason
is testable and was tested.** Everything downstream of §3 inherits that ceiling. The
paper's *measured* sections are in good standing and several are excellent. The
*proof* is not a proof.

**σ(RH is proven by this document) < 1.**

This is not a statement about the quality of the work. It is a statement about one step.

---

## The load-bearing test

§3.2 constructs an amplitude Lagrangian from the Euler product,

```
L(σ, E) = e^{−σE} + e^{−(1−σ)E}
J(σ, E) = −∂L/∂σ = E[e^{−σE} − e^{−(1−σ)E}]
J = 0  iff  σ = ½
```

`engines/noether_derivation.py` runs clean and every internal check passes. The algebra
is correct. **The question is whether the conclusion is about ζ.**

`L` has the form `f(σ) + f(1−σ)`. Any such function is symmetric about σ = ½, so its
derivative vanishes there **identically, for every f**. Tested:

| f | J(½) | J(0.3) | zero at ½? |
|---|---|---|---|
| `e^{−σE}` (the paper's) | 5.55e-11 | 0.693883 | yes |
| `σ³ + 2σ` | −2.22e-10 | 1.200000 | yes |
| `sin(3σ) + σ²` | −0.00e+00 | −2.579368 | yes |
| `e^σ + 1/(σ+2)` | −0.00e+00 | 0.715756 | yes |
| `σ⁵ − 7σ` | −0.00e+00 | 1.160000 | yes |

**Five of five.** The vanishing of J at σ = ½ is a property of the **symmetrisation**,
not of the Euler product. Deriving `L` from prime amplitudes is careful work and it is
real, but it is **not doing the work at the critical step** — the same conclusion follows
from a polynomial, a sine, or anything else fed through `f(σ) + f(1−σ)`.

**What §3 establishes:** σ = ½ is the stationary point of a functional symmetric about
σ = ½. **σ = ∞** — it is a tautology, and tautologies are certain.

**What §3 does not establish:** that the zeros of ζ lie there. No step connects the
stationary point of this constructed `L` to the zero set of ζ. **σ < 1.**

### A second, independent gap

Noether's theorem requires invariance along a **continuous flow**. §3.1 embeds the
reflection in a one-parameter family (ε = 0 identity, ε = 1 reflection, generator
δs = 1−2s), which is the right instinct. But the functional equation ξ(s) = ξ(1−s) holds
**only at the endpoints**. Nothing establishes ξ(s_ε) = ξ(s) for intermediate ε, and
without that the hypothesis of Noether's theorem is not met. The symmetry is genuinely
**discrete (ℤ/2)**, and a path drawn between two points where a symmetry holds is not a
symmetry of the path.

This gap is independent of the first. Closing either would not close the other.

---

## Section-by-section

| § | Claim | σ | Note |
|---|---|---|---|
| 2.1–2.2 | Functional equation, completed ζ | **∞** | Riemann 1859. Bedrock. |
| 2.3 | Noether's theorem as stated | **∞** | Noether 1918. Correctly stated. |
| 2.4 | Berry–Keating H = xp | **5+** | Correctly cited as a *programme*, not a theorem. |
| 3.1 | One-parameter embedding | **2–3** | Construction is valid; invariance along the flow is **not shown**. |
| 3.2 | Currents derived from Euler product | **5+** | Algebra verified, engine runs. |
| 3.3–3.4 | J = 0 ⟺ σ = ½ | **∞ / <1** | ∞ as a statement about symmetric functionals; **<1** as evidence about ζ. |
| 5 | Main theorem (RH) | **<1** | Inherits the §3 ceiling. |
| 6.1–6.3 | Chladni node lines | **2–3** | Good physical analogy. Attractor claim is not established for ζ. |
| 6.4 | Primes are antinodes | **∞ / 2–3** | Explicit formula is ∞ (von Mangoldt); the nodal *reading* is 2–3. Correctly split in the text. |
| **6.5** | **Limaçon / cardioid threshold** | **5+** | **Added today. Measured: 9/10 origin passages match known zeros; per-loop limaçon fits R² 0.86–1.0000; 0/9 in the cardioid band.** Claims only what it measured. |
| 11 | Fermat N-shape, Gauss–Wantzel | **∞ / 2** | Constructibility is ∞; the extinction-class reading is 2. |
| 13 | Abrikosov lattice | **1–2** | Rich analogy. "Topological lock" is asserted, not derived. |
| 14 | Hagedorn ceiling of the Riemann gas | **5+** | Bost–Connes, Julia, Spector correctly cited. The §12.1 tension is recorded rather than hidden — good practice. |

---

## Code

`engines/noether_derivation.py` — runs, exit 0, all internal checks pass. The derivation
of the current forms from the Euler product is genuine and the arithmetic is right. **σ = 5+
for what it computes.** It does not, and does not claim to, bridge to ζ's zeros.

Notebooks 01–08 were not re-executed this pass. **Not assessed.**

---

## What would raise the central σ

In order of decisiveness:

1. **Show the constructed `L` is not arbitrary** — i.e. that some property of ζ's zeros
   follows from *this* `L` and fails for `σ³+2σ`. Until a construction distinguishes the
   Euler-product `L` from a generic symmetrisation at the critical step, §3 cannot carry
   the theorem.
2. **Establish invariance along the flow**, not just at ε ∈ {0,1} — or replace Noether
   with a tool valid for discrete symmetries.
3. **Connect the stationary point to the zero set.** These are different objects. No step
   currently links them.

Item 1 is the one I would attack first, because it is a sharp, falsifiable question with a
cheap test: find any property that the paper's `L` has and an arbitrary symmetrisation
lacks.

---

## Assessment

The document is well-organised, honestly stratified, and unusually careful about marking
`[ESTABLISHED]` / `[HEURISTIC]` / `[THEORETICAL]`. Several sections are genuinely strong —
§6.4's use of the explicit formula, §14's Hagedorn treatment, and §6.5's measurements are
all work I would stand behind.

The proof is not a proof, and the failure is at one identifiable step rather than diffused
through the argument. That is a better position to be in than a vague one: **the gap has a
location, a test, and a criterion for closure.**

The 2026-05-10 valuation gave the Noether application "σ = ∞ for the theorem itself"
with the application's subtlety "addressed below." Having now *tested* the application
rather than reasoned about it, I put it lower. The generic-symmetrisation result is the
specific thing that valuation did not check.
