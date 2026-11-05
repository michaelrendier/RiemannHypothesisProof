# Factoring π out of the Zeta Function

### the critical line as an equatorial great circle, and a catalogue of where ½ emerges

**Cody Michael Allison** (Michael Rendier) · the.wandering.god@gmail.com ·
[github.com/michaelrendier](https://github.com/michaelrendier) · ORCID [0009-0007-7239-6760](https://orcid.org/0009-0007-7239-6760)

**Status: draft 0.1 (2026-11-03) — not for submission.** Licence **GNU GPL v3.0** ([`LICENSE`](LICENSE)).
This repository is a code-and-measurement paper in analytic number theory. It **does not prove the Riemann Hypothesis**; it states
one step as open (§9 of the paper) and labels every other claim by what supports it.

---

## An informal introduction

The completed zeta function ξ(s) = ½ s(s−1) π^(−s/2) Γ(s/2) ζ(s) satisfies the functional equation ξ(s) = ξ(1−s). The reflection
s ↦ 1−s fixes the line Re(s) = ½ for *any* function; what the factor π^(−s/2) does is make ξ itself invariant under it. Delete the
factor and the symmetry is lost (the paper shows this at 30 digits). π is therefore the constant that centres the critical strip.

This repository takes that observation literally. Dividing the strip coordinate by π (θ = πσ) identifies the critical strip
0 < σ < 1 with the polar range 0 < θ < π of a sphere, the reflection with the equatorial reflection θ ↦ π−θ, and the critical line
with the equator, with offset ½ = (π/2)/π. The paper states the limit of this construction at once: it is an **exact change of
variable** and, taken alone, places no zero anywhere. Its use is the geometry it gives — the setting in which the remaining
questions are asked.

Those questions are about the number ½ itself: **where does it come from, and does it come from one place or from many?** The
paper collects the places where ½ appears — as the self-dual point of the reflection, as the unit-diameter circle (r = ½, C = π,
A = π/4), as the minimum of an amplitude Lagrangian built from the Euler product, as the in-phase/quadrature balance (tan = 1), as the
half-integer offset of the zero-counting phase, as the equal-envelope condition of the explicit formula — and classifies each as
*emergent* (the output of a computation whose input did not contain ½), *constructed* (imposed by a symmetric choice), or
*conditional* (true only if RH holds). Only the emergent rows count as evidence about the structure of the problem, and none of
them is evidence that a given zero lies on the line. Around that catalogue the paper records three hypotheses (ζ read as a
Lagrangian; a map of what the flat plot of the zeta spiral overlays, read in the UNS log-polar chart; orthogonality of the *prime wavelength* — the zero
spacing — and the *factor wavelength* λ = p), a table of what is proved, what is imported and where the gap is, and a list of
corrections to earlier versions of this work.

**What is claimed, in one line:** ½ is over-determined by independent routes in the surrounding mathematics; that is a statement
about the problem's structure, not a proof that the zeros obey it. The step that would prove RH — the Berry–Keating
self-adjointness conjecture, or equivalently the identification of ζ with a mode whose nodal set is the equator — remains **open**.

## Contents

| | |
|---|---|
| [`paper/PAPER.md`](paper/PAPER.md) | **the paper, draft 0.1** — §1 Introduction · §2 Scope, labels, UNS · §3 The spine: factoring π out of ζ · §4 Where ½ emerges (catalogue) · §5 Zeta as a Lagrangian · §6 The spiral and the half-integer offset · §7 The flat pathway and its UNS map · §8 The prime wavelength and the coil · §9 Proved / imported / gap · §10 Corrections · §11 Open items · §12 Reproducing · References · Appendix A |
| [`PLAN.md`](PLAN.md) | the build plan that follows from the draft: the software to vendor, the code still to write, the documents to create |
| [`archive/`](archive/README.md) | the repository as it stood on 2026-11-03, moved unchanged; **dated work, not the present state** |
| `LICENSE` | GNU GPL v3.0 |

### What exists and what does not yet

The repository is mid-rewrite. **Exists:** the draft paper; the plan; the archive; and the first built piece, the **flattening map**
of §7 — `rhp/uns.py` (Hardy form, θ′, Gram points, UNS cylinder coordinates), `rhp/flatten.py` (artefact classifier and the
viewport sweep), `rhp/render.py` (dependency-free SVG and the CSV a Blender import reads), `rhp/viewport.py` and `rhp/camera.py`
(vendored from the author's nes-viewport and scad-spatial skills), `tests/test_uns.py` and `tests/test_bigpath.py` (29 tests), `notebooks/01_flattening_map.ipynb` and
`02_machinery_flow_long_pathway.ipynb` (prime wheels, the two halves, the ξ′/ξ flow, the argument-principle accounting, a 3065-zero pathway) and
`03_clock_echo_harness_margins_primegame.ipynb` (the θ-clock, negative real-axis crossings, the echo spectrum = explicit formula, an operator test harness with controls, Lehmer margins, the truncated-balance branches, Conway's PRIMEGAME, and a ten-item anomaly log),
`rhp/bigpath.py`, `flow.py`, `epicycles.py`, `render_big.py`, `clock.py`, `gram.py`, `echo.py`, `harness.py`, `lehmer.py`, `balance.py`, `primegame.py`, `figures/{flattening_map,long_pathway,machinery,flow,clock_echo_harness,balance_margins_primegame}.svg`, `data/pathway.csv`, `data/zeros_1_to_3600.csv`, `requirements.txt`. **Does not exist yet** (see [`PLAN.md`](PLAN.md)): the rest of the
vendored code, notebooks 04–07, the sphere-chart picture, the generated Manual, `paper/make_paper.py`. Two items have
**no code anywhere** today: the limaçon sifting fit and any test of the orthogonality conjecture of §8.

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest -q tests                      # 29 tests, ~35 s
cd notebooks && ../.venv/bin/python make_notebook_01.py \
  && ../.venv/bin/jupyter nbconvert --to notebook --execute --inplace 01_flattening_map.ipynb
```

## How to read the paper

* **Labels.** `ESTABLISHED` — a published theorem or definition, cited. `OURS` — this project's construction or framing.
  `THEORETICAL` — stated, not established. `CONJECTURE` / `OPEN` — named open items. `HEURISTIC` — convergent evidence, not
  deduction. The suffix **`:CALCULATED`** means this repository's own code reproduces the result; it does not mean *derived* or *proved*.
* **Code.** Short mechanisms appear as numbered Listings, extracted verbatim and run cold; sessions as `>>>` transcripts; long or
  established machinery is cited and pointed at. (This is the `cs-paper-code-conventions` discipline; the generator that enforces it
  is among the things still to be written.)
* **Statuses are per claim, not per section.** §9 is the place to start if you want the logical structure; §10 is the place to start
  if you want to know what earlier versions claimed that did not survive re-running.

## For mathematicians looking this over

Three places where a specialist's reading would help most:

1. **§3.2 and §9** — is the θ = πσ map useful for anything beyond bookkeeping? The paper says no on its own; where might it do work?
2. **§7** — the flattening map is built (`rhp/uns.py`, `rhp/flatten.py`, notebook 01, `figures/flattening_map.svg`); the σ ≠ ½ sweep as a picture and the sphere chart are not, and no candidate operator passes the harness.
3. **§8** — a natural inner product on (zero spacings, prime periods) for which the orthogonality conjecture is a well-posed
   statement, and a null model against which a computed value would mean something.

## Assistant note

The paper was written with an assistant — Claude (Anthropic), model `claude-sonnet-5-5` — under the project's
`cs-paper-code-conventions` skill. The author is an independent researcher; he does not tell the assistant the result he hopes for,
as a safeguard (not a guarantee). Earlier documents in `archive/` were written in other sessions, some by earlier web-API
assistants, and disagree with each other; they are kept as the record, not as the paper.

## How to cite

A `CITATION.cff` will accompany the first tagged version. For now: Allison, C. M. (2026), *Factoring π out of the Zeta Function*,
draft 0.1, this repository.
