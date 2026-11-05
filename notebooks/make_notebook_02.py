# SPDX-License-Identifier: GPL-3.0-only
"""Builds notebooks/02_machinery_flow_long_pathway.ipynb."""
import nbformat as nbf
nb = nbf.v4.new_notebook()
M = lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))
C = lambda s: nb.cells.append(nbf.v4.new_code_cell(s.strip("\n")))

M("""
# 02 — The machinery, the flow, and the long pathway (3065 zeros)

**Assistant note.** Written by an AI assistant (Claude) with the author, from the code in `rhp/`. Labels: ESTABLISHED, OURS, THEORETICAL,
CONJECTURE, OPEN; `:CALCULATED` = this code reproduces it. Floating-point checks here are not proofs; nothing places a zero.

**Vantage point.** Orientation: log-polar (UNS) next to the flat plane, with t explicit. Aperture: σ = ½ continuous for t ∈ [20, 3600]
(3065 zeros), plus the strip 0 < σ < 1 for t ≤ 60.

**Four questions.** (1) What machinery draws the pathway, and where are the primes in it? (2) What is the pathway's two-hemisphere
structure, concretely? (3) What does the flow of ξ′/ξ do in the two hemispheres, and what is the bookkeeping for "inflow/outflow"?
(4) What does a long pathway look like in each chart?
""")

C("""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd().parent))
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"svg.fonttype": "none", "font.family": "monospace", "font.size": 8, "axes.grid": False})
from rhp import uns, bigpath, flow, epicycles as ep, render, render_big
zeros_all = np.loadtxt("../data/zeros_1_to_3600.csv")
print(len(zeros_all), "zeros below 3600; first", zeros_all[:3].round(6), "last", zeros_all[-1].round(6))
""")

M("""
## 1. The machinery: prime wheels, and the two halves of the pathway

**1a. Prime wheels (ESTABLISHED, σ > 1).** log ζ(s) = Σ_p Σ_k p^{−ks}/k. Each (p, k) is a wheel of radius p^{−kσ}/k turning at k·ln p
per unit t. In the log plane the pathway *is* the sum of the prime wheels; in the flat plane the same wheels multiply. The table
shows the maximum error over t ∈ [0, 60] after the first K primes.
""")
C("""
t = np.linspace(0, 60, 1201)
print("   K    σ=1.5       σ=1.1")
a = dict(ep.wheel_table(1.5, t, [1, 3, 10, 100])); b = dict(ep.wheel_table(1.1, t, [1, 3, 10, 100]))
for K in (1, 3, 10, 100):
    print(f"{K:5d}  {a[K]:9.2e}  {b[K]:9.2e}")
""")
M("""
At σ = 1.5 a handful of wheels already reproduce log ζ. At σ = 1.1 a hundred are not enough, because the tail Σ_{p>P} p^{−σ} decays
like P^{1−σ}. At σ = ½ the series diverges: the prime wheels cannot draw the pathway on the line. The curve there is drawn by the
**approximate functional equation**, below.

**1b. The two halves (ESTABLISHED: Hardy–Littlewood approximate functional equation).**
ζ(s) ≈ A(s) + χ(s)·B(s), A = Σ_{n≤N} n^{−s}, B = Σ_{n≤N} n^{s−1}, N = ⌊√(t/2π)⌋, χ(s) = 2^s π^{s−1} sin(πs/2) Γ(1−s).
A is the forward wheel sum (radius n^{−σ}, frequency −ln n). χB is the functional-equation mirror. On σ = ½: |χ| = 1 and B = conj(A),
so **the two halves have equal modulus** and ζ = A + e^{−2iθ}·conj(A). Off the line one half dominates: |χB|/|A| ≈ (t/2π)^{½−σ}.
""")
C("""
t2 = np.linspace(20, 120, 25)
for sg in (0.4, 0.5, 0.6):
    A_, cB, chi, zt = ep.afe_parts(sg, t2)
    r = np.abs(cB) / np.abs(A_)
    print(f"σ={sg}:  median |χB|/|A| = {np.median(r):.3f}   (t/2π)^(½−σ) at median t = {(np.median(t2)/2/np.pi)**(0.5-sg):.3f}   max ||A|−|χB|| on the line: {np.max(np.abs(np.abs(A_)-np.abs(cB))):.1e}" if sg == 0.5 else
          f"σ={sg}:  median |χB|/|A| = {np.median(r):.3f}   (t/2π)^(½−σ) at median t = {(np.median(t2)/2/np.pi)**(0.5-sg):.3f}")
""")

M("""
**1c. Where the primes are (ESTABLISHED explicit formula).** ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1−x⁻²). The zeros are the spectrum of the
prime-counting residual ψ(x) − x; equivalently the primes are the spectrum of the pathway's log. With a Gaussian factor e^{−(γw)²/2}
on each zero (width w in ln x) the sum is a **continuous sieve**: a smooth function of x that sharpens into the discrete jumps as w → 0.
Reconstruction from the first K zeros (all 3066):
""")
C("""
x = np.linspace(2.0, 60.0, 1400)
psi = ep.psi_exact(x)
rec = {K: ep.psi_zeros(x, zeros_all[:K], 0.0) for K in (10, 100, 3066)}
smooth = {w: ep.psi_zeros(x, zeros_all, w) for w in (0.02, 0.005)}
for K, y in rec.items():
    print(f"K={K:5d}  rms error vs exact ψ(x) on [2,60]: {np.sqrt(np.mean((y-psi)**2)):.3f}")
for w, y in smooth.items():
    print(f"w={w}  (all zeros) rms difference from exact: {np.sqrt(np.mean((y-psi)**2)):.3f}")
""")

C("""
fig, ax = plt.subplots(2, 3, figsize=(16.8, 9.2))
tt = np.linspace(0, 50, 2001)
# (a) the pathway at five σ (flat)
for sg, c in zip((1.5, 1.1, 0.8, 0.65, 0.5), ("#440154", "#3b528b", "#21918c", "#5ec962", "#c0392b")):
    z = np.array([complex(__import__("mpmath").zeta(__import__("mpmath").mpc(sg, x))) for x in tt[::2]])
    ax[0, 0].plot(z.real, z.imag, color=c, lw=0.8, label=f"σ={sg}")
ax[0, 0].plot([0], [0], "o", color="k", ms=3); ax[0, 0].legend(frameon=False); ax[0, 0].set_aspect("equal"); ax[0, 0].set_title("a  pathway ζ(σ+it), t ≤ 50, flat")
# (b) log ζ true vs wheels (σ=1.5)
tl = np.linspace(0, 60, 1201)
tr = ep.logzeta_true(1.5, tl)
ax[0, 1].plot(tr.real, tr.imag, "k", lw=1, label="log ζ (σ=1.5)")
for K, c in ((3, "#c0392b"), (10, "#3b528b")):
    w_ = ep.logzeta_wheels(1.5, tl, K); ax[0, 1].plot(w_.real, w_.imag, color=c, lw=0.8, label=f"first {K} prime wheels")
ax[0, 1].legend(frameon=False); ax[0, 1].set_title("b  the log plane: wheels add"); ax[0, 1].set_aspect("equal")
# (c) the two halves at σ=½, t = 20..120
tc = np.linspace(20, 120, 2400)
A_, cB, chi, zt = ep.afe_parts(0.5, tc)
ax[0, 2].plot(zt.real, zt.imag, color="k", lw=0.9, label="ζ")
ax[0, 2].plot(cB.real, cB.imag, color="#c0392b", lw=0.5, label="χB (mirror)")
ax[0, 2].plot(A_.real, A_.imag, color="#3b528b", lw=0.5, label="A  (forward wheels)"); ax[0, 2].legend(frameon=False, fontsize=7)
ax[0, 2].set_title("c  σ=½: |A| = |χB| exactly"); ax[0, 2].set_aspect("equal")
# (d) explicit formula: exact vs K zeros
ax[1, 0].plot(x, psi - x, "k", lw=1, label="ψ(x) − x exact")
for K, c in ((10, "#c0392b"), (100, "#3b528b"), (3066, "#5ec962")):
    ax[1, 0].plot(x, rec[K] - x, color=c, lw=0.7, label=f"first {K} zeros")
ax[1, 0].legend(frameon=False, fontsize=7); ax[1, 0].set_title("d  primes from zeros")
# (e) continuous sieve
ax[1, 1].plot(x, psi - x, "k", lw=1, label="discrete")
for w, c in ((0.02, "#c0392b"), (0.005, "#3b528b")):
    ax[1, 1].plot(x, smooth[w] - x, color=c, lw=0.9, label=f"smoothed, w={w}")
ax[1, 1].legend(frameon=False, fontsize=7); ax[1, 1].set_title("e  a continuous sieve (Gaussian in ln x)")
# (f) dominance ratio vs t
for sg, c in ((0.4, "#c0392b"), (0.5, "k"), (0.6, "#3b528b")):
    A2, cB2, _, _ = ep.afe_parts(sg, tc[::10]); ax[1, 2].plot(tc[::10], np.abs(cB2) / np.abs(A2), color=c, lw=0.9, label=f"σ={sg}")
ax[1, 2].set_yscale("log"); ax[1, 2].legend(frameon=False); ax[1, 2].set_title("f  |χB|/|A|: balance on the line")
fig.tight_layout(); fig.savefig("../figures/machinery.svg"); plt.close(fig)
print("wrote figures/machinery.svg")
""")

M("""
## 2. The flow of ξ′/ξ and the bookkeeping

**Definition (ESTABLISHED).** Φ = log ξ(s), L = ξ′/ξ = 1/s + 1/(s−1) − ½ln π + ½ψ(s/2) + ζ′/ζ(s). Velocity (u, v) = (Re L, −Im L), u along σ,
v along t; zeros of ξ are sources with vortex. From ξ(s) = ξ(1−s) and ξ(s̄) = conj ξ(s): L(1−s̄) = −conj L(s), so
**u is odd and v is even under σ ↦ 1−σ**, and |L| is even (equal energy density in the two hemispheres). `:CALCULATED` below.
""")
C("""
pts = [0.3+5j, 0.1+14.5j, 0.45+30j, 0.2+21.5j, 0.35+47j, 0.05+55j]
du, dv = flow.symmetry_residuals(pts)
print(f"max |u(s*)+u(s)| = {du:.1e},  max |v(s*)−v(s)| = {dv:.1e}   (s* = 1 − s̄)")
""")
M("""
**The sign map.** If all zeros are on the line, u > 0 for σ > ½ and u < 0 for σ < ½: outflow from the line in both hemispheres, a stagnation
line at σ = ½. The check on a 60 × 80 grid (σ ∈ [0.05, 0.95], t ∈ [0.5, 60], skipping the line):
""")
C("""
sig = np.concatenate([np.linspace(0.05, 0.45, 40), np.linspace(0.55, 0.95, 40)])
ts = np.linspace(0.5, 60, 60)
U, V = flow.flow_grid(sig, ts)
viol = np.sum(U * (sig - 0.5)[None, :] <= 0)
print("grid points:", U.size, " sign violations of u·(σ−½) > 0:", int(viol), " min u/(σ−½):", float((U / (sig - 0.5)[None, :]).min()).__round__(4))
fig, ax = plt.subplots(1, 2, figsize=(16.8, 6))
sg2 = (sig - 0.5)[None, :]
im = ax[0].pcolormesh(sig, ts, np.log10(U / sg2), cmap="viridis", shading="nearest")
fig.colorbar(im, ax=ax[0], label="log₁₀ (u/(σ−½))  > −∞ everywhere ⇒ outflow")
ax[0].plot(np.full(len(zeros_all[zeros_all < 60]), 0.5), zeros_all[zeros_all < 60], "k.", ms=3)
ax[0].set_xlabel("σ"); ax[0].set_ylabel("t"); ax[0].set_title("margin of outflow: log₁₀ u/(σ−½) (finite = same sign everywhere); black = zeros")
sub = (slice(None, None, 3), slice(None, None, 4))
SS, TT = np.meshgrid(sig, ts)
ax[1].quiver(SS[sub], TT[sub], np.sign(U[sub]) * 1.0, np.sign(V[sub]) * 1.0, np.hypot(U[sub], V[sub]), cmap="viridis", scale=40)
ax[1].axvline(0.5, color="k", lw=0.5); ax[1].set_title("direction of (u, v): u mirrors across σ=½, v does not")
fig.tight_layout(); fig.savefig("../figures/flow.svg"); plt.close(fig)
""")
M("""
**Gauss's law for the flow: the argument principle.** The number of zeros inside the box 0 < σ < 1, 0 < t < T equals θ(T)/π + 1 + S(T),
with S(T) = (1/π) arg ζ(½+iT) followed continuously along 2 → 2+iT → ½+iT. The same count on the line is the number of sign changes of Hardy's Z.
Equality at every T means **no source lies off the equator below T**. This is the check Turing's method makes rigorous (Platt–Trudgian
use it with interval arithmetic); the version here is floating-point, `:CALCULATED`, not a proof.
""")
C("""
rng = np.random.default_rng(3)
mids = (zeros_all[:-1] + zeros_all[1:]) / 2
pick = np.sort(rng.choice(len(mids), 24, replace=False))
Ts = np.concatenate([mids[:6], mids[pick]])
rows = flow.argument_accounting(Ts, zeros_all)
print("      T        S(T)     N_arg=θ/π+1+S   zeros on line   gap")
for r in rows:
    print(f"{r['T']:9.2f} {r['S']:9.4f} {r['N_arg']:14.6f} {r['N_line']:13d} {r['gap']:+10.1e}")
print("max |gap| =", max(abs(r['gap']) for r in rows), " max Δarg per step =", max(r['step_max'] for r in rows))
""")

M("""
## 3. The long pathway

3065 zeros in t ∈ [20, 3600]: about 1500 full turns of the angle −θ(t). Panel A is the flat plot: the ten-pass hub of notebook 01 has become a dense core;
the **envelope** of the overlay is cardioid-like. That is the envelope of a statistical overlay, not the curve (the curve has no cusp).
Panel B is the same curve in the UNS chart (radius log₁₀|ζ|): the zeros are rays into the centre. C and D: the hub fan, S(γₙ) mod 1 — the direction
at the origin of the nth pass is π/2 + π·S(γₙ) (mod π), where S is the zero-counting error. E: S(t) between zeros. F, G: the spring.
""")
C("""
p = bigpath.pathway_big(3600.0, 0.02)
zz = zeros_all[zeros_all > 20]
check = bigpath.zeros_from_Z(p["t"], p["Z"])
print("grid samples:", len(p["t"]), "| zeros found on the grid by sign change:", len(check), "| polished zeros:", len(zz),
      "| max grid-vs-polished:", float(np.max(np.abs(check - zz))))
th = bigpath.theta_fast(zz)
S = ((-th / np.pi) % 1.0) - 0.5
print(f"S(γₙ) mod 1: mean {S.mean():+.4f}, sd {S.std():.4f}, fraction within ±0.25: {np.mean(np.abs(S) < 0.25):.3f}")
mid = (zz[:-1] + zz[1:]) / 2
St = (np.arange(1, len(mid) + 1) + 1) - (bigpath.theta_fast(mid) / np.pi + 1)
print(f"S(t) between zeros: mean {St.mean():+.4f}, sd {St.std():.3f}, max |S| {np.abs(St).max():.3f}")
print("turns of the angle −θ over [20, 3600]:", round(float((bigpath.theta_fast([3600.0])[0] - bigpath.theta_fast([20.0])[0]) / (2*np.pi)), 1))
""")
C("""
render_big.long_pathway_svg(p, zz, "../figures/long_pathway.svg")
from IPython.display import SVG, display
display(SVG("../figures/long_pathway.svg"))
""")
M("""
## 4. What is established here, and what is not

| item | label |
|---|---|
| log ζ is a sum of prime wheels for σ > 1; the table of K | ESTABLISHED, `:CALCULATED` |
| approximate functional equation; |A| = |χB| on σ = ½ (error 1e-14), ratio ≈ (t/2π)^{½−σ} off it | ESTABLISHED, `:CALCULATED` |
| explicit formula; smoothed explicit formula as a continuous sieve (its "difference from exact" is the intended blur) | ESTABLISHED, `:CALCULATED` |
| u odd, v even under σ ↦ 1−σ (residuals ~1e-16) | `:CALCULATED` (follows from the two symmetries) |
| sign of u·(σ−½): 0 violations, 4800 grid points | `:CALCULATED` on a grid; equivalent to RH (Lagarias), **not a proof** |
| argument-principle accounting to T = 3600: gap ~1e-14 at all sampled T | `:CALCULATED`; floating point, not rigorous |
| hub fan = S(γₙ) | `:CALCULATED` (exact consequence of the Hardy form and N(T) = θ/π + 1 + S) |
| the zeros are sources of a flow whose outflow is "factor/prime exchange" of the Two Trees | **THEORETICAL**: the bookkeeping that carries it (J_fwd + J_bwd = 0) is a different conserved quantity and is not computed here |
| the cardioid-like envelope of the flat overlay | observation; no fit made |

Nothing here is evidence about where the zeros are beyond the known verification to 3·10¹².
""")
nbf.write(nb, "02_machinery_flow_long_pathway.ipynb")
