# SPDX-License-Identifier: GPL-3.0-only
"""nodal_extract.svg — the nodal-line extractor on both fields.   (cymatics: sand collects where the field is zero)

a  ξ on the strip σ ∈ [−0.5, 1.5], t ∈ [10, 50]: cymatic sand (density ∝ exp(−Re ξ²/2ε²) with Re ξ scaled by its rms on each row t, ε = 0.12; heaviest on Re ξ = 0), the extracted Re ξ = 0 'rungs' (red) and the
   Im ξ = 0 lines off the equator (blue), the equator σ = ½ (Im ξ ≡ 0).  Rungs cross the equator at the zeros ●; the blue lines leave the equator at the extrema of the heartbeat F′(t) = 0 ◆ and tilt toward larger σ.
b  the prime wave W(u) = (ψ(x) − x)/√x on x ∈ [1000, 1100]: exact nodes (ψ(x) = x, ticks) against the nodes of the wave rebuilt from K zero-tones (K = 100, 3066, 25000).
c  node count in x ∈ [1000, 1500] against K: found nodes, Rice's law for a Gaussian band-limited wave, the exact count 39.
Code drawn: rhp.nodal_extract (zero_set_2d, seg_intersections, sign_nodes, exact_prime_nodes, zero_tone_wave, rice_rate, sand), data experiments/nodal_extract/r1.json.   Run: python -m rhp.atlas_riemann_fermat.nodal_extract
"""
import json, math, os
import numpy as np
from mpmath import mpc
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY, HERE
from .. import nodal as N, nodal_extract as NX

ROOT = HERE.parents[1]


def _grid(t, sig):
    Zh = np.array([[complex(N.xi(mpc(s, tt))) for s in sig] for tt in t])            # σ ≥ ½ only; mirror: ξ(σ+it) = conj ξ(1−σ+it)
    sigfull = np.concatenate([1 - sig[:0:-1], sig]); return sigfull, np.concatenate([np.conj(Zh[:, :0:-1]), Zh], axis=1), Zh


def build():
    G = np.array([float(l.split("\t")[1]) for l in open(ROOT / "data" / "zeros_25000_hp.tsv").read().split("\n")[1:] if l.strip()])
    r1 = json.load(open(ROOT / "experiments" / "nodal_extract" / "r1.json"))
    fig = plt.figure(figsize=(11.8, 7.4))
    ax = fig.add_subplot(1, 3, 1)
    t = np.arange(10, 50.0001, 0.04); sig = np.round(np.arange(0.5, 1.5001, 0.02), 10)
    sf, Zf, Zh = _grid(t, sig)
    ax.pcolormesh(sf, t, NX.sand(Zf.real / np.sqrt(np.mean(Zf.real ** 2, axis=1, keepdims=True)), 0.12), cmap="Greys", shading="auto", vmin=0, vmax=1, rasterized=True)
    Rl, _ = NX.zero_set_2d(sf, t, Zf.real); Il, _ = NX.zero_set_2d(sig[1:], t, Zh[:, 1:].imag)
    for L in Rl: ax.plot(L[:, 0], L[:, 1], color=HOT, lw=1)
    for L in Il: ax.plot(L[:, 0], L[:, 1], color=NORTH, lw=1)
    ax.axvline(.5, color=INK, lw=.6, ls=":")
    zz = G[(G > 10) & (G < 50)]; ax.plot([.5] * len(zz), zz, "o", color=HOT, ms=3.5)
    ex = np.array(N.heartbeat_extrema(50, n=4000)); ex = ex[ex > 10]; ax.plot([.5] * len(ex), ex, "D", color=NORTH, ms=3)
    ax.set_xlabel("σ"); ax.set_ylabel("t"); ax.set_title(f"a  ξ: sand, {len(Rl)} rungs, {len(Il)} Im-lines", loc="left", fontsize=8.5)
    ax = fig.add_subplot(1, 3, 2)
    lo, hi = 1000, 1100; u = np.arange(math.log(lo), math.log(hi), 1e-5); psi, _ = NX.psi_table(hi + 2)
    ex_pos, ex_dir, _k = NX.exact_prime_nodes(lo, hi, psi)
    xs = np.exp(u); floor = np.floor(xs).astype(int); W = (psi[floor] - xs) / np.sqrt(xs)
    ax.plot(xs, W, color=INK, lw=.6, label="exact W (ψ(x) − x)/√x"); ax.axhline(0, color=GREY, lw=.4)
    for K, c in ((100, GREY), (3066, SOUTH), (25000, HOT)):
        w = NX.zero_tone_wave(u, G[:K]); pos, _d = NX.sign_nodes(u, w)
        ax.plot(xs, w, color=c, lw=.6, label=f"K = {K}  ({len(pos)} nodes)")
    ax.plot(ex_pos, np.zeros(len(ex_pos)), "|", color=INK, ms=14, mew=1.4, label=f"exact nodes ({len(ex_pos)})")
    ax.set_xlim(lo, hi); ax.set_ylim(-.5, .5); ax.legend(fontsize=6, frameon=False, loc="lower left"); ax.set_xlabel("x")
    ax.set_title("b  ψ(x) = x: exact vs K zero-tones", loc="left", fontsize=8.5)
    ax = fig.add_subplot(1, 3, 3)
    row = r1["1000-1500"]; Ks = sorted(int(k) for k in row["per_K"])
    ax.loglog(Ks, [row["per_K"][str(k)]["nodes"] for k in Ks], "o-", color=HOT, label="nodes found")
    ax.loglog(Ks, [row["per_K"][str(k)]["rice_predicted"] for k in Ks], "s--", color=GREY, label="Rice's law")
    ax.axhline(row["exact_nodes"], color=INK, lw=1, label=f"exact: {row['exact_nodes']} (jump {row['exact_jump']}, slide {row['exact_slide']})")
    ax.set_xlabel("zero-tones K"); ax.set_ylabel("nodes in x ∈ [1000, 1500]"); ax.legend(fontsize=6.5, frameon=False)
    ax.set_title("c  nodes vs K: saturate; Rice does not", loc="left", fontsize=8.5)
    fig.tight_layout()
    return save(fig, "nodal_extract.svg")


if __name__ == "__main__":
    print(build())
