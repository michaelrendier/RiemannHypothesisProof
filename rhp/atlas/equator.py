# SPDX-License-Identifier: GPL-3.0-only
"""equator.svg — the critical line as the equator; π factored out; the circle versus the winding.

Code drawn:  rhp.sphere (cayley, to_xyz, theta_split, longitude_cayley, winding_hardy, mirror_residual, pole_residues).
Run:         python -m rhp.atlas.equator
"""
import numpy as np
from ._common import plt, save, view, INK, NORTH, SOUTH, HOT, GREY
from .. import sphere as S, bigpath


def build():
    z = np.loadtxt(__file__.rsplit("/rhp/", 1)[0] + "/data/zeros_1_to_3600.csv")
    fig = plt.figure(figsize=(11, 8.6))
    # (a) the Cayley sphere
    ax = fig.add_subplot(2, 2, 1)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("a  the critical line is the equator of w = (s−1)/s", loc="left", fontsize=9)
    u = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(u), np.sin(u), color=GREY, lw=.8)
    for lat in (-60, -30, 30, 60):                                   # parallels, faint
        r, h = np.cos(np.radians(lat)), np.sin(np.radians(lat))
        X, Y, D = view(np.column_stack([r * np.cos(u), r * np.sin(u), np.full_like(u, h)]))
        ax.plot(np.where(D > 0, X, np.nan), np.where(D > 0, Y, np.nan), color=GREY, lw=.3)
    for az in range(0, 180, 30):                                     # meridians, faint
        v = np.linspace(-np.pi / 2, np.pi / 2, 200)
        for sgn in (1, -1):
            P = np.column_stack([sgn * np.cos(v) * np.cos(np.radians(az)), sgn * np.cos(v) * np.sin(np.radians(az)), np.sin(v)])
            X, Y, D = view(P)
            ax.plot(np.where(D > 0, X, np.nan), np.where(D > 0, Y, np.nan), color=GREY, lw=.3)
    E = np.column_stack([np.cos(u), np.sin(u), np.zeros_like(u)])    # the equator
    X, Y, D = view(E)
    ax.plot(np.where(D > 0, X, np.nan), np.where(D > 0, Y, np.nan), color=HOT, lw=2)
    ax.plot(np.where(D <= 0, X, np.nan), np.where(D <= 0, Y, np.nan), color=HOT, lw=1, ls=(0, (3, 3)))
    g = z[:80]
    lon = S.longitude_hardy(g)
    X, Y, D = view(S.xyz_from(np.zeros(len(g)), lon))
    ax.scatter(X[D > 0], Y[D > 0], s=9, color=INK, zorder=5)
    ax.scatter(X[D <= 0], Y[D <= 0], s=7, facecolor="none", edgecolor=INK, lw=.5, zorder=5)
    for name, p3, col in (("s = 1  (N pole of ζ, Λ residue +1)", [0, 0, 1], NORTH), ("s = 0  (S pole, residue −1)", [0, 0, -1], SOUTH)):
        x, y, d = view(np.array([p3]))
        ax.scatter(x, y, s=40, color=col, zorder=6)
        ax.annotate(name, (x[0], y[0]), (x[0] + .12, y[0] + (.1 if p3[2] > 0 else -.12)), fontsize=8, color=col)
    lo = float(S.longitude_hardy(20.0)[0])
    pair = S.xyz_from([S.latitude_cayley(0.75 + 20j), -S.latitude_cayley(0.75 + 20j)], [lo, lo])
    X, Y, D = view(pair)
    ax.plot(X, Y, color=INK, lw=.8); ax.scatter(X, Y, s=16, color=[NORTH, SOUTH], zorder=7)
    ax.annotate("s ↦ 1−s̄ (t = 20): same longitude,\nopposite latitude", (X[0], Y[0]), (.45, .78), fontsize=7.2, arrowprops=dict(arrowstyle="-", lw=.5))
    ax.text(0, -1.3, "longitude = Hardy phase θ(t) mod 2π (the spring end-on); dots = first 80 zeros, open = far side:\n"
            "they scatter about two antipodal carrier points (sd 0.25 half-turn; uniform = 0.29), they do not sit on them.  north |w|<1: σ>½.", ha="center", fontsize=7.2)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.4, 1.25)
    # (b) the flat strip it came from
    ax = fig.add_subplot(2, 2, 2)
    ax.set_title("b  the same thing flat: the critical strip", loc="left", fontsize=9)
    ax.axvspan(0.5, 1, color=NORTH, alpha=.10); ax.axvspan(0, 0.5, color=SOUTH, alpha=.10)
    ax.axvline(.5, color=HOT, lw=2)
    ax.scatter(np.full(40, .5), z[:40], s=6, color=INK, zorder=4)
    ax.text(.75, 3, "north\n(Euler product\nconverges for σ>1)", ha="center", fontsize=7.5, color=NORTH)
    ax.text(.25, 3, "south\n(its mirror)", ha="center", fontsize=7.5, color=SOUTH)
    ax.set_xlim(0, 1); ax.set_ylim(0, 80)
    ax.set_xlabel("σ   (σ = 0 and σ = 1 are the two poles; the edges collapse to points on the sphere)")
    ax.set_ylabel("t")
    ax2 = ax.secondary_xaxis("top", functions=(lambda s: np.degrees(S.latitude_linear(s)), lambda l: np.radians(l) / np.pi + .5))
    ax2.set_xlabel("latitude π(σ−½)  — π re-introduced (degrees)", fontsize=7.5)
    # (c) π factored out of θ
    ax = fig.add_subplot(2, 2, 3)
    t = np.linspace(0.5, 100, 400)
    gp, pp, tot = S.theta_split(t[::8])
    ax.plot(t[::8], gp, color=NORTH, lw=1.2, label="Im lnΓ(¼+it/2): carries the t·ln t growth")
    ax.plot(t[::8], pp, color=SOUTH, lw=1.2, label="−(t/2)·ln π: the π part, a constant-rate rotation")
    ax.plot(t[::8], tot, color=INK, lw=1.6, label="θ(t) = sum  (mpmath siegeltheta)")
    ax.set_title("c  factoring π out of the phase: θ(t) = Γ-part + π-part", loc="left", fontsize=9)
    ax.set_xlabel("t"); ax.set_ylabel("rad"); ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax.text(100, tot[-1] * .15, f"max |sum − θ| = {np.max(np.abs(gp + pp - tot)):.1e}", ha="right", fontsize=7.2)
    # (d) what the flat circle forgets
    ax = fig.add_subplot(2, 2, 4)
    tt = np.geomspace(0.2, 3600, 500)
    lon = S.longitude_cayley(tt) / np.pi
    wind = bigpath.theta_fast(tt) / np.pi
    ax.semilogx(tt, lon, color=HOT, lw=1.6, label="Cayley longitude/π (left): wraps once, then crowds at w = 1")
    ax.set_ylabel("Cayley longitude / π", color=HOT); ax.set_ylim(-1.05, 1.05)
    ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
    ax2.semilogx(tt, wind, color=INK, lw=1.6, label="Hardy winding θ/π (right): the helix, one half-turn per zero")
    ax2.set_ylabel("θ(t)/π", color=INK)
    ax.set_title("d  the equator forgets the winding; the winding is the zero index", loc="left", fontsize=9)
    ax.set_xlabel("t")
    h = ax.get_legend_handles_labels(); h2 = ax2.get_legend_handles_labels()
    ax.legend(h[0] + h2[0], h[1] + h2[1], fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(.5, -.2))
    fig.suptitle("The equatorial construction (CALCULATED: rhp/sphere.py; a change of coordinates, not evidence about RH)", fontsize=9.5, x=.02, ha="left")
    fig.tight_layout()
    return save(fig, "equator.svg")


if __name__ == "__main__":
    print(build())
