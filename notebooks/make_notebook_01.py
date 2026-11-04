# SPDX-License-Identifier: GPL-3.0-only
"""Builds notebooks/01_flattening_map.ipynb (then execute with nbconvert; see docs/Reproducing.md)."""
import nbformat as nbf

nb = nbf.v4.new_notebook()
M = lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))
C = lambda s: nb.cells.append(nbf.v4.new_code_cell(s.strip("\n")))

M("""
# 01 — Mapping the flattening: ζ(½+it) flat, log-polar and as a spring

**Assistant note.** This notebook was written by an AI assistant (Claude) working with the author, from the code in `rhp/`.
Labels follow the paper: ESTABLISHED (published), OURS, THEORETICAL, CONJECTURE, OPEN; the suffix `:CALCULATED` means the
code here reproduces the statement and the step-by-step is shown. Nothing below places a zero.

**Vantage point** (observer-position). *Orientation:* log-polar in the ζ-plane (radius log₁₀|ζ|, angle arg ζ) with t kept as
an explicit third axis. *Aperture:* the equator σ = ½, continuous in t, 0 ≤ t ≤ 50. σ-off-the-line is not in this notebook.

**Question.** The familiar picture (Wikipedia's *polar graph of ζ(½+it)*) is a flat image of a curve that lives in three
variables. Which features of that picture are the curve, and which are the projection?
""")

C("""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd().parent))          # repo root, so `import rhp` works from notebooks/
import numpy as np
from rhp import uns, flatten, render
p = uns.pathway(50.0, 4001)
zeros = uns.zeros_below(50.0)
grams = uns.gram_points(50.0)
print(len(p["t"]), "samples;", len(zeros), "zeros below 50;", len(grams), "Gram points below 50")
""")

M("""
## 1. The flat curve is a signed-radius polar curve  (ESTABLISHED form, `:CALCULATED` checks)

Hardy: ζ(½+it) = e^{−iθ(t)} Z(t), with θ the Riemann–Siegel phase and Z real. So the flat point has **angle −θ(t)** and
**signed radius Z(t)** — passing through the origin continues on the opposite side, which is why there are petals and no corners.

**Listing 1 — θ′ and the spring coordinates (`rhp/uns.py`, extracted)**
""")
C("""
import inspect
from rhp import uns
print(inspect.getsource(uns.theta_prime))
print(inspect.getsource(uns.cylinder))
""")

M("Check the form on the grid: the imaginary part of e^{iθ}ζ should vanish, and Z should equal mpmath's `siegelz`.")
C("""
from mpmath import mpf, siegelz
w = np.exp(1j * p["theta"]) * p["zeta"]
print("max |Im(e^{iθ} ζ)| =", float(np.max(np.abs(w.imag))))
print("max |Z − siegelz|  =", max(abs(p["Z"][i] - float(siegelz(mpf(p["t"][i])))) for i in range(0, 4001, 400)))
""")

M("""
## 2. The hub: ten zeros, one point

All ten passes through (0, 0) are drawn on top of each other. The direction in which pass n crosses the origin is
arg(i·ζ′(½+iγₙ)) mod π, and Hardy's form says it equals **−θ(γₙ) mod π**.
""")
C("""
rows = flatten.hub_passes(p, zeros)
print(" n    γₙ       direction   −θ(γₙ) mod π   agree   Z sign before→after")
for r in rows:
    print(f"{r['n']:2d} {r['gamma']:8.4f} {r['direction']:10.4f} {r['minus_theta_mod_pi']:12.4f}   {r['agree']}   {r['sign_before']:+d} → {r['sign_after']:+d}")
""")

M("""
## 3. The orbit reading: angular velocity is exactly −θ′(t)

From ζ = e^{−iθ}Z: ζ̄·ζ′ = Z(Z′ − iθ′Z), so **Im(ζ̄ ζ′) = −θ′(t) Z(t)²**, and the angular velocity dφ/dt = −θ′(t)
independent of the radius. θ′(t) = ½ ln(t/2π) + O(1/t²) changes sign at t ≈ 2π: before that the curve turns the other way.
The origin is a point of zero angular momentum: each zero is a radial pass.
""")
C("""
from mpmath import mpc, zeta
print("    t     Im(ζ̄ ζ′)    −θ′Z²")
for t in (3.0, 10.0, 14.0, 20.0, 27.7, 41.3):
    s = mpc(0.5, t); z = complex(zeta(s)); dz = 1j * complex(zeta(s, derivative=1))
    Z = float(siegelz(mpf(t)))
    print(f"{t:6.1f} {(z.conjugate()*dz).imag:11.5f} {-uns.theta_prime(t)*Z*Z:11.5f}")
lc = flatten.lost_clock(p, zeros)
print("numeric dφ/dt vs −θ′: max deviation", lc["max_abs_dev"], "; θ′ changes sign at t =", lc["sign_change_t"], "(2π = 6.283)")
""")

M("""
## 4. The real axis is the Gram points

The flat curve crosses the real axis where θ(t) = nπ. The black dots on the real axis in panel A of the figure are the Gram points gₙ for n ≥ 0 (the early arc also crosses at g₋₁ ≈ 9.667, not drawn).
""")
C("""
for g in grams[:6]:
    z = uns.zeta_line(g)
    print(f"g = {g:9.5f}   θ/π = {uns.theta(g)/np.pi:6.3f}   ζ = {z.real:8.4f} {z.imag:+.1e}i")
""")

M("""
## 5. Artefact classification — what the flattening does

*Crossings.* A self-intersection of the flat polyline is a **false crossing** (the two parameters differ). The ten hub passes
account for C(10,2) = 45 of them; the others are false crossings away from the origin.
""")
C("""
cr = flatten.false_crossings(p)
hub = [c for c in cr if c["hub"]]; away = [c for c in cr if not c["hub"]]
print("total", len(cr), "| hub", len(hub), "= C(10,2) =", 10*9//2, "| away from the origin", len(away))
print("   t1      t2      dt      point")
for c in away:
    print(f"{c['t1']:6.2f}  {c['t2']:6.2f}  {c['dt']:6.2f}   ({c['x']:+.3f}, {c['y']:+.3f})")
""")
M("*Petal size.* The flat plot shows the linear maximum |ζ| of each petal; the UNS (decimal-log) size is far less dramatic.")
C("""
print("petal      t range        max |ζ|    log10 max")
for r in flatten.petals(p, zeros):
    print(f"{r['petal']:3d}   {r['t_from']:6.2f}–{r['t_to']:6.2f}   {r['max_abs']:8.3f}   {r['max_log10']:+8.3f}")
""")
M("""
*Missed petals (Gram's law).* A Gram interval [gₙ, gₙ₊₁) normally holds one zero. Where it does not, the flat curve fails to cross
the real axis between two zeros. None occurs below t = 50; the first is far out.
""")
C("""
g = flatten.gram_law(300.0)
print(g["intervals"], "Gram intervals up to t = 300; Gram's law fails in", len(g["failures"]))
for f in g["failures"]:
    print(f"  interval n={f['n']}  [{f['g_lo']:.3f}, {f['g_hi']:.3f})  zeros = {f['zeros']}")
print("mirror half: max |ζ(½−it) − conj ζ(½+it)| =", flatten.mirror_check())
""")

M("""
## 6. The viewport sweep

A window of 500 samples scrolls over the pathway with the nes-viewport model; the sweep test asserts every sample is shown
(clamp, flush, thumb, sweep). Each row is what one window of the flat plot hides.
""")
C("""
print("  oy      t range       zeros  false  hub")
for r in flatten.window_map(p, zeros, cr, vh=500):
    print(f"{r['oy']:5d}  {r['t_from']:5.1f}–{r['t_to']:5.1f}  {r['zeros']:5d}  {r['false_crossings']:5d}  {r['hub_crossings']:4d}")
""")

M("""
## 7. The figure

A: the flat plot. B: the same curve end-on in UNS (log-radius). C, D: the unrolled radius and clock. E–G: the spring in three
canonical poses (front, iso, top; orthographic). The top pose is panel B: the flat plot is the end-on view with the radius
exponentiated. The spring's zeros sit on its axis.
""")
C("""
from IPython.display import SVG, display
render.flattening_svg(p, zeros, grams, "../figures/flattening_map.svg")
render.write_csv(p, "../data/pathway.csv")
display(SVG("../figures/flattening_map.svg"))
""")

M("""
## 8. What this establishes, and what it does not

| item | label |
|---|---|
| ζ(½+it) = e^{−iθ}Z; flat angle −θ, signed radius Z | ESTABLISHED (Hardy), checks `:CALCULATED` |
| exit direction at the nth zero = −θ(γₙ) mod π | `:CALCULATED` (follows from the form; 10/10 here) |
| Im(ζ̄ζ′) = −θ′Z², dφ/dt = −θ′ | `:CALCULATED` |
| real-axis crossings = Gram points | ESTABLISHED (Gram), `:CALCULATED` here |
| hub = C(10,2) pairwise crossings; 21 false crossings below t = 50 | `:CALCULATED` at this sampling |
| Gram's law fails first in the interval n = 125 (t ≈ 281–282) | ESTABLISHED (known first violation), reproduced |
| the σ ≠ ½ sweep, the sphere (azimuth = t mod 2π) chart as a picture | OPEN (not in this notebook) |

Nothing here is evidence about where the zeros are. It maps what the flat picture overlays, and what it keeps.
""")

nbf.write(nb, "01_flattening_map.ipynb")
