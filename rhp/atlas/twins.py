# SPDX-License-Identifier: GPL-3.0-only
"""twins.svg — the two books (sieve / un-sieve), the circles' twin constant, the clock–wheel resonance.

a  Telperion × Laurelin: each n ≤ 10⁵ at (π(spf), π(gpf)) — when the sieve first kills it, when the wheels first complete it.
b  the circles' ratio Π p(p−2)/(p−1)² against the largest circle P: it converges to the Hardy–Littlewood constant C₂ = 0.66016.
c  twin pairs below x: counted against 2C₂Σ1/ln²x (a heuristic prediction, not a proof).
d  gaps between consecutive twin pairs, in units of the local mean (ln p)²/(2C₂), against an exponential law.
e  wheel n against the clock: θ(t) − t ln n is stationary at t = 2πn² — the one wheel at rest relative to the clock at each height.
Code drawn:  rhp.wrap, rhp.uns.   Run: python -m rhp.atlas.twins
"""
import math
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import wrap as W


def build():
    fig = plt.figure(figsize=(11.5, 8.4))
    # a
    ax = fig.add_subplot(2, 3, 1)
    n = 100000
    spf, gpf = W.spf_gpf(n)
    pi_ = np.cumsum(W.sieve(n))
    m = np.arange(4, n + 1)
    comp = m[spf[m] != m]
    ax.scatter(pi_[spf[comp]], pi_[gpf[comp]], s=.3, color=GREY, alpha=.35, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.plot([1, 1e4], [1, 1e4], color=INK, lw=.7)
    ax.set_title("a  sieve × un-sieve: composites ≤ 10⁵", loc="left", fontsize=9)
    ax.set_xlabel("π(spf)  Telperion: when it dies"); ax.set_ylabel("π(gpf)  Laurelin: when it is born")
    ax.text(1.2, 3000, "prime powers sit on the\ndiagonal (not drawn: primes)", fontsize=7)
    # b
    ax = fig.add_subplot(2, 3, 2)
    Ps = np.unique(np.geomspace(5, 100000, 60).astype(int))
    ax.semilogx(Ps, [W.singular_series_from_circles(int(P)) for P in Ps], color=INK, lw=1.5)
    ax.axhline(W.TWIN_C2, color=HOT, lw=1, ls="--"); ax.text(6, W.TWIN_C2 + .006, "C₂ = 0.66016", color=HOT, fontsize=7.5)
    ax.set_title("b  the circles' twin constant", loc="left", fontsize=9)
    ax.set_xlabel("largest circle P"); ax.set_ylabel("Π p(p−2)/(p−1)²")
    # c
    ax = fig.add_subplot(2, 3, 3)
    N = 10 ** 7
    tw = W.twin_pairs(N)
    xs = np.geomspace(100, N, 80)
    ax.loglog(xs, [np.searchsorted(tw, x, side="right") for x in xs], color=INK, lw=1.6, label="counted")
    ax.loglog(xs, [W.hl_twin_count(int(x)) for x in xs], color=SOUTH, lw=1, ls="--", label="2C₂ Σ 1/ln²x")
    ax.set_title("c  twin pairs below x", loc="left", fontsize=9); ax.set_xlabel("x"); ax.legend(fontsize=7, frameon=False)
    # d
    ax = fig.add_subplot(2, 3, 4)
    z = W.twin_gaps_normalised(N)
    ax.hist(z, bins=60, range=(0, 6), density=True, color=NORTH, alpha=.6, label=f"twin gaps (mean {z.mean():.3f}, CV {z.std() / z.mean():.2f})")
    x = np.linspace(0, 6, 200); ax.plot(x, np.exp(-x), color=HOT, lw=1.4, label="exponential (memoryless)")
    ax.set_title("d  lengths between twin pairs", loc="left", fontsize=9); ax.set_xlabel("gap / local mean (ln p)²/2C₂"); ax.legend(fontsize=7, frameon=False)
    # e
    ax = fig.add_subplot(2, 3, 5)
    for k, (nn, col) in enumerate(zip([2, 3, 5, 7], [NORTH, SOUTH, HOT, INK])):
        t0 = W.resonance_height(nn)
        t = np.linspace(max(2.0, t0 * .3), t0 * 1.7, 160)
        ph = W.relative_phase(nn, t); ph0 = float(W.relative_phase(nn, [t0])[0])
        ax.plot(t / t0, ph - ph0, color=col, lw=1.3, label=f"n = {nn}, rests at t = 2π·{nn}² = {t0:.1f}")
        ax.scatter([1], [0], color=col, s=14, zorder=5)
    ax.set_ylim(-1.5, 3); ax.axvline(1, color=GREY, lw=.5)
    ax.set_title("e  wheel n against the Hardy clock", loc="left", fontsize=9)
    ax.set_xlabel("t / 2πn²"); ax.set_ylabel("θ(t) − t ln n  (shifted to 0 at rest)"); ax.legend(fontsize=6.5, frameon=False, loc="lower center")
    ax = fig.add_subplot(2, 3, 6); ax.axis("off")
    ax.text(0, .95, "reading", fontsize=9, fontweight="bold")
    ax.text(0, .45, "a  every composite has two coordinates; primes close both\n"
                    "   books at once.\n"
                    "b  the 0-spots of the circles are not independent: the\n"
                    "   correction is exactly C₂.\n"
                    "c  the circles predict the count to 0.4%: a prediction,\n"
                    "   not a proof.\n"
                    "d  gaps between twin pairs are nearly memoryless: only\n"
                    "   their mean is a smooth spring.\n"
                    "e  at height t the one wheel at rest against the clock\n"
                    "   is n = √(t/2π): the clock visits the wheels in order.",
            fontsize=7.6, va="center")
    fig.suptitle("Two books, the twin circles, and the clock–wheel resonance (CALCULATED: rhp/wrap.py)", fontsize=9.5, x=.01, ha="left")
    fig.tight_layout()
    return save(fig, "twins.svg")


if __name__ == "__main__":
    print(build())
