# SPDX-License-Identifier: GPL-3.0-only
"""wrap.svg — the primes on the two springs.

a,b  the sieve spring: integers wrapped on circles of circumference 6 and 30; primes marked; every twin pair straddles the 0 spot.
c    the allowed-centre fraction of the primorial circle at each level (the 0 spot thins).
d    the wheel spring: twin wheels drift apart at t·ln(1 + 2/p); the beat carries the gap.
e    the tension: the Hardy spring tightens like 1/ln t; the prime spikes (t = ln n) crowd like eᵗ/t.
Code drawn:  rhp.wrap.   Run: python -m rhp.atlas.wrap
"""
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import wrap as W


def _circle(ax, m, title, mark_twin=True):
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, loc="left", fontsize=9)
    ang = lambda n: np.pi / 2 - 2 * np.pi * n / m                    # 0 spot at the top, clockwise
    u = np.linspace(0, 2 * np.pi, 300); ax.plot(np.cos(u), np.sin(u), color=GREY, lw=.6)
    P = set(int(p) for p in W.primes_to(10 * m))
    for r in range(m):
        # residue r is "prime-like" if some prime > m sits on it (coprime to m) or r is the 0 spot
        import math
        coprime = math.gcd(r, m) == 1
        x, y = np.cos(ang(r)), np.sin(ang(r))
        if r == 0: ax.scatter(x, y, s=70, color=HOT, zorder=5); ax.text(x * .5, y * .5, "0 spot\n(multiples of %d)" % m, ha="center", va="center", fontsize=7, color=HOT)
        elif coprime: ax.scatter(x, y, s=42, color=INK, zorder=4); ax.text(x * 1.16, y * 1.16, str(r), ha="center", va="center", fontsize=7)
        else: ax.scatter(x, y, s=8, color=GREY, zorder=3)
    if mark_twin:
        # centres of twin pairs mod m: arcs joining c−1 and c+1 around the centre
        ok = sorted(set(int(c) for c in W.allowed_centres(5)[1])) if m == 30 else [0]
        for c in ok:
            t = np.linspace(ang(c - 1), ang(c + 1), 20)
            ax.plot(.82 * np.cos(t), .82 * np.sin(t), color=NORTH, lw=2.2, solid_capstyle="round")


def build():
    fig = plt.figure(figsize=(11.5, 8.6))
    _circle(fig.add_subplot(2, 3, 1), 6, "a  integers wrapped on a circle of 6")
    _circle(fig.add_subplot(2, 3, 2), 30, "b  wrapped on 30: arcs = allowed twin centres")
    # c — the 0 spot thins
    ax = fig.add_subplot(2, 3, 3)
    Ps = [5, 7, 11, 13, 17, 19, 23]
    fr = [W.allowed_centre_count(P) / W.primorial(P) for P in Ps]
    ax.semilogy(range(len(Ps)), fr, "o-", color=INK)
    for i, P in enumerate(Ps): ax.text(i, fr[i] * (1.18 if i % 2 else .8), f"{W.allowed_centre_count(P)}/{W.primorial(P)}", fontsize=6, ha="center", va="bottom" if i % 2 else "top")
    ax.set_xticks(range(len(Ps))); ax.set_xticklabels([f"{P}#" for P in Ps])
    ax.set_title("c  allowed centres on the P# circle: the 0 spot thins", loc="left", fontsize=9)
    ax.set_ylabel("fraction of the circle"); ax.set_ylim(0.015, 0.2)
    n, ok = W.straddle_check(10 ** 7)
    ax.text(.02, .03, f"twin pairs ≤ 10⁷ with p>3: {n}; centre ≡ 0 (mod 6): {ok}", transform=ax.transAxes, fontsize=7)
    # d — twin wheels
    ax = fig.add_subplot(2, 3, 4)
    t = np.linspace(0, 60, 300)
    for (p, q), col in zip([(101, 103), (1019, 1021), (10007, 10009), (101, 107)], [NORTH, SOUTH, HOT, GREY]):
        ax.plot(t, [W.wheel_phase_lag(p, q, x) for x in t], color=col, lw=1.4, label=f"({p},{q}): antiphase at t = {W.antiphase_height(p, q):.0f}")
    ax.axhline(np.pi, color=INK, lw=.6, ls="--"); ax.text(1, np.pi + .05, "π: antiphase", fontsize=7)
    ax.set_title("d  twin wheels part at t·ln(q/p)", loc="left", fontsize=9)
    ax.set_xlabel("height t"); ax.set_ylabel("phase lag (rad)"); ax.set_ylim(0, 4); ax.legend(fontsize=6.8, frameon=False, loc="center left")
    # e — the tension
    ax = fig.add_subplot(2, 3, 5)
    tt = np.geomspace(14, 3600, 80)
    ax.semilogx(tt, W.hardy_pitch(tt), color=INK, lw=1.6)
    ax.set_title("e  the zeros' spring tightens", loc="left", fontsize=9)
    ax.set_xlabel("height t"); ax.set_ylabel("height of one 2π turn  = 2π/θ′(t)")
    ax = fig.add_subplot(2, 3, 6)
    th = np.linspace(1, 20, 200)
    ax.semilogy(th, W.prime_density_in_t(th), color=HOT, lw=1.6)
    ax.set_title("f  the primes' spikes (t = ln n) crowd like eᵗ/t", loc="left", fontsize=9)
    ax.set_xlabel("t = ln n"); ax.set_ylabel("prime spikes per unit t")
    fig.suptitle("Primes on the springs (CALCULATED: rhp/wrap.py).  Two different springs: the sieve (integers) and the wheel (Euler product)",
                 fontsize=9, x=.01, ha="left")
    fig.tight_layout()
    return save(fig, "wrap.svg")


if __name__ == "__main__":
    print(build())
