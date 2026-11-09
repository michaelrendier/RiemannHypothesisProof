# SPDX-License-Identifier: GPL-3.0-only
"""polarisation.svg — the two counter-rotating waves of a term, and where they are a pure standing wave.

Each term n of the approximate functional equation is a pair:  forward n^{−s}  and mirror χ(s)·n^{−(1−s)}.
On σ = ½ both have modulus n^{−½}: equal, counter-rotating — a LINEAR polarisation, a pure standing wave, the nodal state.
Off the equator the amplitudes differ by r = |χ|·n^{2σ−1}; the two phasors trace an ELLIPSE (axial ratio |1−r|/(1+r)).
r = 1 off the equator only at n ≈ √(t/2π): the Riemann–Siegel cut-off.  So the standing-wave set is a CROSS: the equator
(every term) and the horizontal line n = √(t/2π) (every latitude).

Code drawn:  rhp.sphere (abs_chi, amplitude_ratio, polarisation, ellipse_axes).   Run: python -m rhp.atlas_riemann_fermat.polarisation
"""
import numpy as np
from matplotlib.patches import Ellipse
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import sphere as S


def build():
    fig = plt.figure(figsize=(11.5, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1, 1])
    # (a) the ellipses: term n = 3 at t = 100 across latitude
    ax = fig.add_subplot(gs[0]); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("a  term n = 1, t = 100", loc="left", fontsize=9)
    n, t = 1, 100.0
    for k, sg in enumerate([0.5, 0.4, 0.3, 0.15, 0.0]):
        a1 = n ** (-sg); a2 = S.abs_chi(sg, t) * n ** (sg - 1)
        A, B = S.ellipse_axes(a1, a2); sc = 0.46 / (a1 + a2)
        col = HOT if sg == 0.5 else (NORTH if sg > .5 else SOUTH)
        ax.add_patch(Ellipse((k * 1.05 + .5, 0), 2 * A * sc, 2 * B * sc, fill=False, ec=col, lw=1.6))
        ax.text(k * 1.05 + .5, -.62, f"σ={sg}\nr={a2 / a1:.2f}\nb/a={S.polarisation(a2 / a1):.2f}", ha="center", fontsize=7)
    ax.text(.5, .62, "σ=½: a line\n(standing wave)", ha="center", fontsize=7, color=HOT)
    ax.set_xlim(-.2, 5.6); ax.set_ylim(-1.05, 1); ax.text(2.6, .85, "pair polarisation by latitude", ha="center", fontsize=7)
    # (b, c) the cross
    for j, T in enumerate((100.0, 1000.0)):
        ax = fig.add_subplot(gs[1 + j])
        sig = np.linspace(0.02, 0.98, 97); ns = np.arange(1, 41)
        M = np.array([[S.polarisation(S.amplitude_ratio(m, s, T)) for s in sig] for m in ns])
        im = ax.imshow(M, origin="lower", extent=(sig[0], sig[-1], .5, 40.5), aspect="auto", cmap="magma_r", vmin=0, vmax=1)
        N = np.sqrt(T / (2 * np.pi))
        ax.axvline(.5, color=HOT, lw=1.2); ax.axhline(N, color=INK, lw=1, ls="--")
        ax.text(.51, 38.5, "equator", color=HOT, fontsize=7); ax.text(.03, N + .8, f"n = √(t/2π) = {N:.2f}", fontsize=7)
        ax.set_title(f"{'bc'[j]}  b/a at t = {T:g} (pale = standing)", loc="left", fontsize=9)
        ax.set_xlabel("σ   (north →)"); ax.set_ylabel("term n")
        if j == 1: fig.colorbar(im, ax=ax, fraction=.046, label="b/a  (0 linear … 1 circular)")
    fig.suptitle("Two counter-rotating waves per term: the equator and the RS cut-off form the nodal cross  (CALCULATED: rhp/sphere.py)",
                 fontsize=9, x=.01, ha="left")
    fig.tight_layout()
    return save(fig, "polarisation.svg")


if __name__ == "__main__":
    print(build())
