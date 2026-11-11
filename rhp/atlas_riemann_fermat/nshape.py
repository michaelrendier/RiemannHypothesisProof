# SPDX-License-Identifier: GPL-3.0-only
"""nshape.svg — The Shape: the Extinguishable Bulk (Laurelin) around the holes (Telperion).

a  the Ω-onion: share of n ≤ x in each layer Ω(n) = k (k = 1 the primes, the core; k = 2 semiprimes wrap them; …), x = 10² … 10⁷. The core thins like 1/ln x; the mode moves from layer 2 to layer 3 at ≈ 10⁵.
b  the two factoring conics: the fraction of odd composites N ≤ x reached by Euler's circle (two primitive representations as a sum of two squares) against Fermat's hyperbola a² − b² = N (every odd composite, no prime). The circle's reach is the composites built from primes ≡ 1 mod 4 only.
c  the tower: Z, Z[i] (D2), the Hurwitz quaternions (D4), the integral octonions (E8): irreducible classes per rational prime p (log axis) and the share of rational primes that stay holes (1, ½, 0, 0).
Code drawn: rhp.nshape.   Data: experiments/nshape_shape.json.   Run: python -m rhp.atlas_riemann_fermat.nshape
"""
import json, math
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY, HERE
from .. import nshape as S

ROOT = HERE.parents[1]


def build():
    d = json.load(open(ROOT / "experiments" / "nshape_shape.json"))
    fig = plt.figure(figsize=(11.8, 4.4))
    ax = fig.add_subplot(1, 3, 1)
    X = 10 ** 7; om = S.big_omega(X); xs = np.unique(np.round(np.logspace(2, 7, 40)).astype(int)); K = 8
    cols = plt.cm.viridis(np.linspace(0.05, 0.95, K + 1))
    sh = np.array([[v / (x - 1) for v in S.layer_counts(om, [x], kmax=K)[x]] for x in xs])
    ax.stackplot(np.log10(xs), sh.T, colors=cols, labels=[f"Ω = {k}" + ("+" if k == K else "") for k in range(1, K + 1)] + [])
    ax.set_xlabel("log₁₀ x"); ax.set_ylabel("share of 2 ≤ n ≤ x"); ax.legend(fontsize=5.8, frameon=False, loc="upper right", ncol=2)
    ax.set_title("a  the Ω-onion (core = primes)", loc="left", fontsize=8.5)
    ax = fig.add_subplot(1, 3, 2)
    NB = 10 ** 6; P = S.primes_upto(NB); isp = np.zeros(NB + 1, bool); isp[P] = True; n = np.arange(NB + 1)
    cnt = S.circle_primitive_reps(NB); oddc = (n % 2 == 1) & ~isp & (n > 1); two = oddc & (cnt >= 2)
    grid = np.unique(np.round(np.logspace(2, 6, 40)).astype(int)); frac = [two[:x + 1].sum() / oddc[:x + 1].sum() for x in grid]
    ax.semilogx(grid, frac, color=HOT, lw=1.4, label="Euler: circle x² + y² = N"); ax.semilogx(grid, [1] * len(grid), color=NORTH, lw=1.4, label="Fermat: hyperbola a² − b² = N")
    ax.set_xlabel("x"); ax.set_ylabel("share of odd composites ≤ x reached"); ax.set_ylim(0, 1.05); ax.legend(fontsize=6.5, frameon=False, loc="center right")
    ax.set_title("b  two factoring conics", loc="left", fontsize=8.5)
    ax = fig.add_subplot(1, 3, 3)
    p = np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47])
    ax.semilogy(p, np.ones(len(p)), "o-", color=GREY, ms=3, label="Z  (units 2; holes 100 %)")
    ax.semilogy(p, np.where(p % 4 == 1, 2, 1), "s-", color=SOUTH, ms=3, label="Z[i]  (units 4; holes 50 %)")
    ax.semilogy(p, np.where(p == 2, 1, p + 1), "^-", color=NORTH, ms=3, label="Hurwitz  (units 24; holes 0)")
    ax.semilogy(p, 1 + p ** 3, "d-", color=HOT, ms=3, label="E8 octonions  (units 240; holes 0)")
    ax.set_xlabel("rational prime p"); ax.set_ylabel("irreducible classes of norm p"); ax.legend(fontsize=6, frameon=False, loc="upper left")
    ax.set_title("c  the tower: primes absorbed, new holes born", loc="left", fontsize=8.5)
    fig.tight_layout()
    return save(fig, "nshape.svg")


if __name__ == "__main__":
    print(build())
