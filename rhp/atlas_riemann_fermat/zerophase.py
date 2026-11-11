# SPDX-License-Identifier: GPL-3.0-only
"""zerophase.svg — zero-phase windings against Archimedean (polygonal) windings: which arrangement do the primes follow?

a  the critical line is the equator: each zero-tone x^ρ = e^{(½+iγ)u} lifted by inverse stereographic projection (equator at |z| = 1) is a loxodrome (constant bearing arctan 2γ, pitch arctan 1/2γ); zeros 1, 2, 3.
b  bearing along u for zero 1: stereographic lift constant (87.974°), the linear-in-ln-n lift coded in lattice_planes.py wobbles ±1.4°.
c  the zero-phase scan A(t) = |Σ (Λ(n) − 1) n^{−½ − it}| / ln x at x = 5·10⁴ against the zero ordinates (ticks): 29 lines, height ≈ 0.92, background 0.15.
d  the Archimedean scan A(h), θ = 2π h √n (h = 1 Sacks/Ulam; m-gon at h = √(2/(m−2))): nothing above the random-phase level 0.53 (dashed); polygon values ticked.
Code drawn: rhp.zerophase.   Data: experiments/zerophase_probe.json, zerophase_scans.npz.   Run: python -m rhp.atlas_riemann_fermat.zerophase
"""
import json
import numpy as np
from ._common import plt, save, view, INK, NORTH, SOUTH, HOT, GREY, HERE
from .. import zerophase as Z

ROOT = HERE.parents[1]


def build():
    d = json.load(open(ROOT / "experiments" / "zerophase_probe.json")); sc = np.load(ROOT / "experiments" / "zerophase_scans.npz")
    g = np.array([float(l.split()[-1]) for l in open(ROOT / "data" / "zeros_25000_hp.tsv") if l[0].isdigit()])
    fig = plt.figure(figsize=(11.8, 8.2))
    ax = fig.add_subplot(2, 2, 1)
    t = np.linspace(0, 2 * np.pi, 200)
    for lat in np.radians([-60, -30, 30, 60]):
        x, up, dp = view(np.c_[np.cos(lat) * np.cos(t), np.cos(lat) * np.sin(t), np.full_like(t, np.sin(lat))]); ax.plot(x, up, color="#d6d6de", lw=0.5)
    x, up, dp = view(np.c_[np.cos(t), np.sin(t), np.zeros_like(t)]); ax.plot(x, up, color="#c9a227", lw=1.4, label="critical line = equator (|z| = 1)")
    ax.plot(np.cos(t), np.sin(t), color=GREY, lw=0.8)
    u = np.linspace(0.02, 6.5, 6000)
    for k, col in ((1, NORTH), (2, SOUTH), (3, HOT)):
        x, up, dp = view(Z.lift_stereo(Z.zero_spiral(u, g[k - 1]))); front = dp > 0
        ax.plot(np.where(front, x, np.nan), np.where(front, up, np.nan), color=col, lw=0.6, label=f"zero {k}: γ = {g[k-1]:.2f}, pitch {Z.pitch_deg(g[k-1]):.2f}°")
    ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=6.2, frameon=False, loc="lower left")
    ax.set_title("a  x^ρ is a loxodrome on the sphere; the critical line is the equator", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 2, 2); u = np.linspace(0.5, 6, 4000); g1 = g[0]
    ax.plot(u[:-1], np.degrees(Z.bearing(Z.lift_stereo(Z.zero_spiral(u, g1)))), color=NORTH, label="stereographic lift (conformal)")
    ax.plot(u[:-1], np.degrees(Z.bearing(Z.lift_linear(u, g1 * u, 6.0))), color=SOUTH, label="linear in ln n (lattice_planes._lift)")
    ax.set_xlabel("u = ln x"); ax.set_ylabel("bearing from the meridian (deg), zero 1"); ax.legend(fontsize=6.5, frameon=False)
    ax.set_title("b  one lift keeps the zero-tone's bearing, the other does not", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 2, 3); m = sc["ts"] < 60
    ax.plot(sc["ts"][m], sc["Az"][m], color=INK, lw=0.7); ax.plot(g[(g > 10) & (g < 60)], np.full(((g > 10) & (g < 60)).sum(), -0.06), "|", color=HOT, ms=7)
    ax.set_xlabel("t"); ax.set_ylabel("A(t)"); ax.set_title("c  zero-phase scan: a line at every zero ordinate (ticks)", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(sc["hs"], sc["Aa"], color=INK, lw=0.5); ax.axhline(0.53, color=GREY, ls="--", lw=0.8)
    for mm, c in ((3, HOT), (4, NORTH), (5, SOUTH), (6, "#2e8b57")):
        ax.axvline(Z.polygon_h(mm), color=c, lw=0.8, alpha=0.8, label=f"{mm}-gon  h = {Z.polygon_h(mm):.3f}")
    ax.set_xlabel("h  (θ = 2π h √n)"); ax.set_ylabel("A(h)"); ax.legend(fontsize=6, frameon=False, ncol=2)
    ax.set_title("d  Archimedean scan: random-phase noise (dashed = 0.53)", loc="left", fontsize=8.5)
    fig.tight_layout()
    return save(fig, "zerophase.svg")


if __name__ == "__main__":
    print(build())
