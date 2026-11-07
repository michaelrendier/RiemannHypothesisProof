# SPDX-License-Identifier: GPL-3.0-only
"""spirals.svg — the spiral is not flat and not static: the pathway and its prime-side twin as 3-D conical helices.

a  ζ(½+it) in the cylindrical UNS chart (angle arg ζ, radius (2/π)arctan|ζ|, height t): the pathway THREADS the t-axis (ζ = 0) at each zero
b  the same at σ = ½, ¾, 1: off the equator the axis is a hole that widens (min |ζ| on the window printed); the pitch tightens with t in all three
c  the partial sums Σ n^{−s} about ζ(s), height u = ln N: a CONE at σ = ½ (flare e^{u/2}), a CYLINDER at σ = 1, a FUNNEL at σ = 3/2 — logarithmic spirals
d  the prime side: the conical helices of the explicit formula (first zeros, flare e^{u/2}) and their sum; prime heights u = ln p marked on the axis
e  a flat pancake coil (POE: 17 turns, 0.8 mm) with the zeros as beads, one turn per two zeros: what a flat spiral keeps (angle, index) and loses (height, radius)
f  reading
Code drawn: rhp.coil, rhp.unsmap.   Run: python -m rhp.atlas.spirals
"""
import numpy as np
from matplotlib.collections import LineCollection
from ._common import plt, save, view, INK, NORTH, SOUTH, HOT, GREY
from .. import coil as K


def _line3(ax, P, c=None, cmap="viridis", lw=1.0, color=None, elev=18, azim=-60, label=None, scale=1.0):
    x, up, depth = view(P * scale, elev=elev, azim=azim)
    if c is None:
        ax.plot(x, up, color=color, lw=lw, label=label); return
    pts = np.column_stack([x, up]).reshape(-1, 1, 2)
    lc = LineCollection(np.concatenate([pts[:-1], pts[1:]], axis=1), cmap=cmap, linewidths=lw)
    lc.set_array(np.asarray(c)[:-1]); ax.add_collection(lc); ax.autoscale_view()


def build():
    zeros = np.loadtxt(K.__file__.replace("rhp/coil.py", "data/zeros_1_to_3600.csv"))
    fig = plt.figure(figsize=(11.5, 8.4))
    ts = np.linspace(10, 60, 2400)
    # a
    ax = fig.add_subplot(2, 3, 1)
    P, z = K.zeta_cylinder(ts, 0.5)
    sc = np.array([1, 1, 0.045])
    ax.plot(*view(np.array([[0, 0, 10], [0, 0, 60]]) * sc)[:2], color=GREY, lw=.8, ls="--")
    _line3(ax, P, c=ts, lw=1.1, scale=sc)
    zz = np.array([[0, 0, g] for g in zeros[:10]]) * sc
    ax.plot(*view(zz)[:2], "o", color=HOT, ms=3.5)
    ax.set_title("a  ζ(½+it) around the t-axis (zeros ● on it)", loc="left", fontsize=8.5); ax.axis("off")
    # b
    ax = fig.add_subplot(2, 3, 2)
    d = K.axis_distance(np.linspace(10, 60, 3000))
    for sg, col in [(0.5, HOT), (0.75, SOUTH), (1.0, NORTH)]:
        Pq, _ = K.zeta_cylinder(ts[::2], sg)
        _line3(ax, Pq, color=col, lw=.9, scale=sc, label=f"σ = {sg}:  min|ζ| = {d[sg]:.3f}")
    ax.plot(*view(np.array([[0, 0, 10], [0, 0, 60]]) * sc)[:2], color=GREY, lw=.8, ls="--")
    ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(0, 0.02)); ax.set_title("b  off the line the axis is a widening hole", loc="left", fontsize=8.5); ax.axis("off")
    # c
    ax = fig.add_subplot(2, 3, 3)
    t0 = 40.0
    for sg, col in [(0.5, HOT), (1.0, GREY), (1.5, NORTH)]:
        C = K.partial_sum_cone(t0, sg, 3000, 60)
        C[:, :2] /= max(np.abs(C[:, :2]).max(), 1e-9)
        _line3(ax, C * np.array([1, 1, 0.25]), color=col, lw=.8, label=f"σ = {sg}")
    ax.legend(fontsize=7, frameon=False, loc="upper right"); ax.set_title("c  partial sums: cone / cylinder / funnel", loc="left", fontsize=8.5); ax.axis("off")
    # d
    ax = fig.add_subplot(2, 3, 4)
    u = np.linspace(0, np.log(40), 5000)
    for k, col in zip(range(3), [NORTH, SOUTH, GREY]):
        H = K.helix(zeros[k], u); _line3(ax, H * np.array([.35, .35, 1]), color=col, lw=.7)
    Cn = K.explicit_cone(zeros[:200], u)
    _line3(ax, Cn * np.array([.35, .35, 1]), color=INK, lw=1.0)
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        ax.plot(*view(np.array([[0, 0, np.log(p)]]) * np.array([.35, .35, 1]))[:2], "|", color=HOT, ms=6)
    ax.set_title("d  explicit formula: helices (colour), sum (black), u=ln p ticks", loc="left", fontsize=8); ax.axis("off")
    # e
    ax = fig.add_subplot(2, 3, 5)
    x, y, phi, r = K.pancake()
    ax.plot(x, y, color=GREY, lw=.7)
    xb, yb, th, rb = K.zeros_on_coil(zeros[:34])
    ax.plot(xb, yb, "o", color=HOT, ms=3)
    for name, (kidx, g) in K.tap_zeros(zeros).items():
        ax.plot([xb[kidx - 1]], [yb[kidx - 1]], "s", color=INK, ms=5); ax.annotate(name, (xb[kidx - 1], yb[kidx - 1]), fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    rl = K.radius_law(zeros[:34])
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"e  flat coil, 34 zeros (r² idx {rl['index_vs_angle']:.3f})", loc="left", fontsize=8)
    # f
    ax = fig.add_subplot(2, 3, 6); ax.axis("off")
    ax.text(0, 1, "Reading\n\n"
            "• Flat is not enough: angle and index survive,\n  height and radius are lost.\n"
            "• a: the pathway is a helix about the t-axis whose\n  radius is Z(t); it passes THROUGH the axis only at\n  zeros. The axis is the 'attractor' (ζ = 0).\n"
            "• b: off the line the axis is a hole that widens\n  with σ−½ (printed): the pull is strongest on the line.\n"
            "• c,d: the same e^{u/2} flare on both sides. The\n  partial-sum spiral about ζ(s) is a log spiral:\n  d ln r/dφ = −(1−σ)/t.\n"
            "• The primes are where all the conical helices of\n  d come into phase: no winding rule is chosen.\n"
            "• e: radius = zero index makes the Hardy angle an\n  exact Archimedean (pancake) spiral up to S(T).\n\n"
            "Not static: pitch 2π/θ′ tightens as t grows; the\nradius is Z(t), not a constant.", va="top", fontsize=7.4, family="monospace")
    fig.tight_layout()
    return save(fig, "spirals.svg")


if __name__ == "__main__":
    print(build())
