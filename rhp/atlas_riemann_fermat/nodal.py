# SPDX-License-Identifier: GPL-3.0-only
"""nodal.svg — the nodal map of ξ on both hemispheres: Re ξ = 0 (red) and Im ξ = 0 (blue).

a  the strip σ ∈ [−1.5, 2.5], t ∈ [0, 50]:  the blue set contains the whole equator σ = ½ (Im ξ = 0 there); the red 'rungs'
   cross the equator only at the zeros, level with t = γ.  The picture is mirror-symmetric about σ = ½ (ξ(s) = ξ(1−s)).
b  zoom on the first five zeros, with the heartbeat's extrema ◆ (F′ = 0): saddles of the analytic landscape on the equator.
c  the heartbeat F(t) = ξ(½+it) with its zeros (rungs) and extrema (saddles).
d  log10 |χ(σ+it)| (ζ(s) = χ(s)ζ(1−s)): zero on the equator, the weighting that turns the north's picture into the south's.
e  ζ itself on the real side: Re ζ = 0 and Im ζ = 0 for σ ∈ [−9, 2], t ∈ [0, 14] — the trivial zeros at −2, −4, −6, −8 are on the real axis
   in the SOUTH; ξ has none there (Γ(s/2) cancels them), so they belong to χ, not to the equatorial structure.
Code drawn:  rhp.nodal, rhp.sphere.abs_chi.   Run: python -m rhp.atlas_riemann_fermat.nodal
"""
import numpy as np
from mpmath import mp, mpc, zeta
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import nodal as N
from ..sphere import abs_chi

ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719, 43.327073, 48.005151, 49.773832]


def build():
    fig = plt.figure(figsize=(11.5, 8.6))
    # a
    ax = fig.add_subplot(2, 3, (1, 4))
    sg = np.linspace(-1.5, 2.5, 161); tv = np.linspace(0.05, 50, 330)
    Z = N.grid(sg, tv)
    ax.contour(sg, tv, Z.imag, levels=[0], colors=NORTH, linewidths=.8)
    ax.contour(sg, tv, Z.real, levels=[0], colors=HOT, linewidths=.8)
    ax.axvline(.5, color=INK, lw=.6, ls=":")
    ax.plot([.5] * len(ZEROS), ZEROS, "k.", ms=4)
    ax.plot([], [], color=NORTH, label="Im ξ = 0 (contains the equator)"); ax.plot([], [], color=HOT, label="Re ξ = 0 (rungs)")
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    ax.set_xlabel("σ   (south ← ½ → north)"); ax.set_ylabel("t"); ax.set_title("a  nodal map of ξ, both hemispheres", loc="left", fontsize=9)
    # b
    ax = fig.add_subplot(2, 3, 2)
    sg2 = np.linspace(-0.2, 1.2, 141); tv2 = np.linspace(12, 34, 220)
    Z2 = N.grid(sg2, tv2)
    ax.contour(sg2, tv2, Z2.imag, levels=[0], colors=NORTH, linewidths=.9)
    ax.contour(sg2, tv2, Z2.real, levels=[0], colors=HOT, linewidths=.9)
    ex = [e for e in N.heartbeat_extrema(34, 1400) if e > 12]
    ax.plot([.5] * len(ex), ex, "D", color=INK, ms=4, label="extrema of F (saddles)")
    ax.plot([.5] * 5, ZEROS[:5], "o", mfc="none", color=INK, ms=5, label="zeros")
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    ax.set_xlabel("σ"); ax.set_ylabel("t"); ax.set_title("b  rungs and saddles", loc="left", fontsize=9)
    # c
    ax = fig.add_subplot(2, 3, 3)
    t = np.linspace(0.5, 50, 1200); F = np.array([N.F(x) for x in t]) * np.exp(np.pi * t / 4)
    ax.plot(t, F, color=INK, lw=1)
    ax.axhline(0, color=GREY, lw=.5)
    ax.plot(ZEROS, [0] * len(ZEROS), "o", mfc="none", color=HOT, ms=5)
    ee = N.heartbeat_extrema(50, 1400); ax.plot(ee, [N.F(x) * np.exp(np.pi * x / 4) for x in ee], "D", color=SOUTH, ms=3.5)
    ax.set_xlabel("t"); ax.set_ylabel("F(t)·e^{πt/4},  F = ξ(½+it)"); ax.set_title("c  the heartbeat on the equator", loc="left", fontsize=9)
    # d
    ax = fig.add_subplot(2, 3, 5)
    sg3 = np.linspace(-1.5, 2.5, 100); t3 = np.linspace(10, 300, 100)
    L = np.log10(np.array([[abs_chi(s, tt) for s in sg3] for tt in t3]))
    im = ax.imshow(L, origin="lower", aspect="auto", extent=[sg3[0], sg3[-1], t3[0], t3[-1]], cmap="PuOr_r", vmin=-3, vmax=3)
    ax.axvline(.5, color="k", lw=.7); fig.colorbar(im, ax=ax, shrink=.75, label="log10|χ|")
    ax.set_xlabel("σ"); ax.set_ylabel("t"); ax.set_title("d  the weight that mirrors N into S", loc="left", fontsize=9)
    # e
    ax = fig.add_subplot(2, 3, 6)
    sg4 = np.linspace(-9, 2, 220); t4 = np.linspace(-0.0, 14, 200)
    with mp.workdps(15):
        Zz = np.array([[complex(zeta(mpc(s, tt))) if abs(mpc(s, tt) - 1) > 1e-9 else 1e9 for s in sg4] for tt in t4])
    ax.contour(sg4, t4, Zz.imag, levels=[0], colors=NORTH, linewidths=.8)
    ax.contour(sg4, t4, Zz.real, levels=[0], colors=HOT, linewidths=.8)
    ax.plot([-2, -4, -6, -8], [0] * 4, "ks", ms=4, label="trivial zeros"); ax.plot([1], [0], "k*", ms=7, label="pole s = 1")
    ax.axvline(.5, color=INK, lw=.6, ls=":"); ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    ax.set_xlabel("σ"); ax.set_ylabel("t"); ax.set_title("e  ζ on the real side: trivial zeros", loc="left", fontsize=9)
    fig.tight_layout()
    return save(fig, "nodal.svg")


if __name__ == "__main__":
    print(build())
