# SPDX-License-Identifier: GPL-3.0-only
"""flow.svg — the extinction flow on the integers, and its Laplace dual.

a  survivor fraction of the sieve against the cutoff time τ = ln y, with Mertens' e^{−γ}/τ
b  Telperion: Φ(N,y)·ln y / N  against U = ln N / ln y, for four N, on the Buchstab function ω(U)   (delay-differential, delay 1)
c  Laurelin:  Ψ(N,y) / N       against U, for four N, on the Dickman function ρ(U)
d  the gain of the sieve book in the dual plane:  ln|1+ω̂(s′)| = Re E₁(s′),  s′ = σ′ + it  (log branch at s′ = 0 = ζ's pole, north)
e  reciprocity  s′ρ̂(s′)·(1+ω̂(s′)) = 1 checked from the NUMERIC ω and ρ along σ′ = 1  (error is the grid, first order in h)
Code drawn:  rhp.extinction, rhp.wrap.  Run: python -m rhp.atlas_riemann_fermat.flow
"""
import math
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import wrap as W
from .. import extinction as X


def build(N=2_000_000):
    spf, gpf = W.spf_gpf(N)
    ur, rho = X.dickman(); uo, om = X.buchstab()
    fig = plt.figure(figsize=(11.5, 8.2))
    # a
    ax = fig.add_subplot(2, 3, 1)
    ys = np.unique(np.geomspace(2, N, 160).astype(int))
    cs = np.cumsum(np.bincount(np.minimum(spf[2:N + 1], N), minlength=N + 1))   # counts by spf value
    S = np.array([(N - 1 - cs[y]) + 1 for y in ys]) / N
    ax.loglog(np.log(ys), S, color=INK, lw=1.6, label="Φ(N,y)/N  counted, N = 2·10⁶")
    mer = np.array([X.mertens(int(y)) for y in ys[ys < 20000]])
    ax.loglog(np.log(ys[ys < 20000]), mer, color=SOUTH, lw=1, ls="--", label="Π_{p≤y}(1−1/p)")
    tt = np.log(ys[1:]); ax.loglog(tt, np.exp(-X.GAMMA) / tt, color=HOT, lw=.8, ls=":", label="e^{−γ}/τ")
    ax.set_xlabel("τ = ln y   (time of the flow)"); ax.set_ylabel("survivors / N")
    ax.set_title("a  the sieve as a flow", loc="left", fontsize=9); ax.legend(fontsize=6.5, frameon=False)
    # b, c
    cols = [GREY, NORTH, SOUTH, INK]
    Ns = [10 ** 4, 10 ** 5, 10 ** 6, N]
    for k, (name, curve_u, curve) in enumerate([("b  Telperion: survivors → Buchstab ω(U)", uo, om), ("c  Laurelin: smooth → Dickman ρ(U)", ur, rho)]):
        ax = fig.add_subplot(2, 3, 2 + k)
        m = (curve_u >= 1) & (curve_u <= 8) if k == 0 else (curve_u <= 7)
        ax.plot(curve_u[m], curve[m], color=HOT, lw=1.4, label="ω(U)" if k == 0 else "ρ(U)", zorder=5)
        for Nn, c in zip(Ns, cols):
            sp, gp = spf[:Nn + 1], gpf[:Nn + 1]
            Us = np.linspace(1.15, 7 if k else 6, 40)
            ysx = np.unique(np.round(Nn ** (1 / Us)).astype(int))
            ysx = ysx[ysx >= 2]
            Uv = math.log(Nn) / np.log(ysx)
            if k == 0:
                v = np.array([X.survivors(sp, Nn, y) * math.log(y) / Nn for y in ysx])
            else:
                v = np.array([X.smooth(gp, Nn, y) / Nn for y in ysx])
            ax.plot(Uv, v, "o", ms=2.2, color=c, label=f"N = 10^{math.log10(Nn):.1f}")
        if k == 1: ax.set_yscale("log"); ax.set_ylim(1e-5, 1.2)
        ax.set_xlabel("U = ln N / ln y"); ax.set_title(name, loc="left", fontsize=9); ax.legend(fontsize=6.5, frameon=False)
    # d
    ax = fig.add_subplot(2, 3, 4)
    sg = np.linspace(-1.5, 2.0, 140); tv = np.linspace(-12, 12, 160)
    G = np.array([[X.log_gain(complex(a, b)) for a in sg] for b in tv])
    im = ax.contourf(sg, tv, np.clip(G, -3, 3), levels=21, cmap="RdBu_r")
    ax.contour(sg, tv, G, levels=[0], colors="k", linewidths=.7)
    ax.plot([0], [0], "k*", ms=8); ax.axvline(0, color=GREY, lw=.5)
    fig.colorbar(im, ax=ax, shrink=.8, label="ln|1+ω̂| = Re E₁(s′)")
    ax.set_xlabel("σ′ = σ − 1   (damping per ln n)"); ax.set_ylabel("t   (frequency in ln n)")
    ax.set_title("d  the dual plane: gain of the sieve", loc="left", fontsize=9)
    # e
    ax = fig.add_subplot(2, 3, 5)
    ts = np.linspace(0, 30, 61)
    for h, c in [(1 / 500, GREY), (1 / 2000, NORTH), (1 / 8000, INK)]:
        u_r, r_ = X.dickman(14, h); u_o, o_ = X.buchstab(14, h)
        ax.semilogy(ts, [X.reciprocal_residual(1 + 1j * t, u_r, r_, u_o, o_) for t in ts], color=c, lw=1.2, label=f"h = 1/{int(round(1 / h))}")
    ax.set_xlabel("t   on σ′ = 1"); ax.set_ylabel("|s′ρ̂(1+ω̂) − 1|")
    ax.set_title("e  reciprocity from the numeric ω, ρ", loc="left", fontsize=9); ax.legend(fontsize=7, frameon=False)
    # reading
    ax = fig.add_subplot(2, 3, 6); ax.axis("off")
    ax.text(0, 1, "Reading\n\n"
            "• The sieve and the un-sieve are CONTINUOUS\n  flows in τ = ln y; the primes are the clicks.\n"
            "• Their delay equations (delay 1 in U) are the\n  'something it avoids': ω rings down to e^{−γ}.\n"
            "• In the dual plane the two books are exact\n  reciprocals:  s′ρ̂ · (1+ω̂) = 1,\n  the Laplace twin of ζ · (1/ζ) = 1.\n"
            "• The dual variable of u = ln x is the s of ζ.\n  σ′ = damping, t = frequency in ln n.\n\n"
            "ESTABLISHED: ω, ρ, Mertens.\nOURS-CALCULATED: the reciprocity (checked\nabove; may be known).",
            va="top", fontsize=7.8, family="monospace")
    fig.tight_layout()
    return save(fig, "flow.svg")


if __name__ == "__main__":
    print(build())
