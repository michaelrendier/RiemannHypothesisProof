# SPDX-License-Identifier: GPL-3.0-only
"""net.svg — the four coordinates of the whole process, unfolded as the SIX squares of a hypercube's faces (a 'net').

Axes:  u = ln n  (size of the number)   τ = ln y  (time of the flow, cutoff)   σ′ = σ−1 (damping)   t (frequency in ln n).
Each square is one PAIR of axes; no single square shows the process, the six together do (observer-position: change one
parameter at a time).  Read the figure as a cube net:   top row  (u,τ) (u,σ′) (u,t) — the three squares that contain u;
                                                         bottom row (τ,σ′) (τ,t) (σ′,t) — the three that do not.
 (u,τ)  each integer n ≤ 3000 at its two clocks:  death  ln spf(n)  (Telperion, blue) and birth  ln gpf(n)  (Laurelin, orange); primes on the diagonal.
 (u,σ′) energy of the first n ≤ e^u integers  Σ n^{−2σ}  (log10): finite for σ>½ (σ′>−½), logarithmic at ½, divergent below.
 (u,t)  the dual lattices: primes are VERTICAL lines u = ln p (the leaf-drops of the Screw), the first zeros are HORIZONTAL lines t = γ.
 (τ,σ′) sieve survivors Π_{p≤e^τ}(1 − p^{−1−σ′}) (log10): settles for σ′>0, falls like 1/τ at σ′ = 0 (Mertens), faster below.
 (τ,t)  phase of the Euler partial product P(½+it) with the wheels up to e^τ switched on — the spring seen from the side.
 (σ′,t) gain of the sieve book Re E₁(s′) (the dual plane of flow.svg).
Code drawn: rhp.wrap, rhp.extinction.   Run: python -m rhp.atlas_riemann_fermat.net
"""
import math
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import wrap as W
from .. import extinction as X

ZEROS = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187, 43.3271, 48.0052, 49.7738]


def build():
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 8.0))
    # (u,τ)
    ax = axs[0, 0]
    N = 3000
    spf, gpf = W.spf_gpf(N)
    n = np.arange(4, N + 1)
    comp = n[spf[n] != n]
    ax.scatter(np.log(comp), np.log(spf[comp]), s=1.2, color=NORTH, label="death ln spf (Telperion)", rasterized=True)
    ax.scatter(np.log(comp), np.log(gpf[comp]), s=1.2, color=SOUTH, label="birth ln gpf (Laurelin)", rasterized=True)
    pr = n[spf[n] == n]; ax.scatter(np.log(pr), np.log(pr), s=3, color=INK, label="primes (both clocks at once)")
    ax.set_xlabel("u = ln n"); ax.set_ylabel("τ = ln y"); ax.legend(fontsize=6, frameon=False, loc="upper left")
    ax.set_title("(u,τ)  two clocks per integer", loc="left", fontsize=9)
    # (u,σ')
    ax = axs[0, 1]
    us = np.linspace(0.5, 11, 120); sg = np.linspace(-1.0, 2.0, 120)
    E = np.log10(np.array([[np.sum(np.arange(1, int(math.exp(u)) + 1, dtype=float) ** (-2 * (s + 1))) if math.exp(u) < 6e4 else
                            (np.log(math.exp(u)) if abs(s + .5) < 1e-9 else (math.exp(u) ** (-2 * s - 1) - 1) / (-2 * s - 1) + 1) for u in us] for s in sg]))
    im = ax.imshow(E, origin="lower", aspect="auto", extent=[us[0], us[-1], sg[0], sg[-1]], cmap="magma", vmin=0, vmax=2.5)
    ax.axhline(-0.5, color="w", lw=.8, ls="--"); ax.text(0.8, -0.45, "σ = ½ : energy diverges like ln N", color="w", fontsize=6.5)
    ax.set_xlabel("u = ln N"); ax.set_ylabel("σ′ = σ − 1"); ax.set_title("(u,σ′)  energy of the first N integers", loc="left", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.75, label="log10 Σ n^{−2σ}")
    # (u,t)
    ax = axs[0, 2]
    for p in W.primes_to(150):
        ax.axvline(math.log(p), color=NORTH, lw=.5, alpha=.7)
    for g in ZEROS:
        ax.axhline(g, color=SOUTH, lw=.6, alpha=.8)
    ax.set_xlim(0, math.log(150)); ax.set_ylim(0, 52)
    ax.set_xlabel("u = ln n   (primes: vertical, u = ln p)"); ax.set_ylabel("t   (zeros: horizontal, t = γ)")
    ax.set_title("(u,t)  the dual lattices", loc="left", fontsize=9)
    # (τ,σ')
    ax = axs[1, 0]
    P = W.primes_to(20000).astype(float)
    taus = np.linspace(0.7, math.log(20000), 140); sg2 = np.linspace(-0.6, 2.0, 120)
    lp = np.log(P)
    cum = {}
    M = np.empty((len(sg2), len(taus)))
    for i, s in enumerate(sg2):
        terms = np.log10(np.abs(1 - P ** (-(1 + s))) + 1e-300)
        c = np.cumsum(terms)
        M[i] = np.interp(taus, lp, c)
    im = ax.imshow(M, origin="lower", aspect="auto", extent=[taus[0], taus[-1], sg2[0], sg2[-1]], cmap="viridis", vmin=-4, vmax=0)
    ax.axhline(0, color="w", lw=.8, ls="--"); ax.text(0.9, 0.05, "σ′ = 0: ~ e^{−γ}/τ", color="w", fontsize=6.5)
    ax.set_xlabel("τ = ln y"); ax.set_ylabel("σ′"); ax.set_title("(τ,σ′)  survivors of the sieve", loc="left", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.75, label="log10 Π(1 − p^{−1−σ′})")
    # (τ,t)
    ax = axs[1, 1]
    tv = np.linspace(0.5, 50, 400)
    Pp = W.primes_to(3000).astype(float)
    ph = np.zeros((len(Pp), len(tv)))
    acc = np.zeros(len(tv), complex)
    for i, p in enumerate(Pp):
        acc += -np.log(1 - p ** (-(0.5 + 1j * tv)))
        ph[i] = np.angle(np.exp(acc))
    taus2 = np.log(Pp)
    im = ax.imshow(ph.T, origin="lower", aspect="auto", extent=[taus2[0], taus2[-1], tv[0], tv[-1]], cmap="twilight")
    for g in ZEROS: ax.axhline(g, color="w", lw=.3, alpha=.6)
    ax.set_xlabel("τ = ln y   (wheels up to e^τ on)"); ax.set_ylabel("t"); ax.set_title("(τ,t)  Euler phase on the equator", loc="left", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.75, label="arg P_y(½+it)")
    # (σ',t)
    ax = axs[1, 2]
    sg3 = np.linspace(-1.5, 2.0, 120); t3 = np.linspace(-12, 12, 140)
    G = np.array([[X.log_gain(complex(a, b)) for a in sg3] for b in t3])
    im = ax.contourf(sg3, t3, np.clip(G, -3, 3), levels=21, cmap="RdBu_r"); ax.plot([0], [0], "k*", ms=7)
    ax.set_xlabel("σ′"); ax.set_ylabel("t"); ax.set_title("(σ′,t)  gain Re E₁ of the sieve book", loc="left", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.75, label="ln|1+ω̂|")
    fig.suptitle("net of the process: the six faces of (u, τ, σ′, t)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    return save(fig, "net.svg")


if __name__ == "__main__":
    print(build())
