# SPDX-License-Identifier: GPL-3.0-only
"""wave.svg — the Heartbeat waveform derived from the two sets, and the 'wind' it makes.

a  Telperion's wave W_T (from the primes: (ψ(x)−x)/√x) and Laurelin's W_L (from the composites: (⌊x⌋−1−ψ)/√x): exact mirror images, no zeros used
b  the same W_T against the sum of the first 10 / 100 / 3066 zero tones H(u) = −2Σcos(γu − arg ρ)/|ρ| (ESTABLISHED explicit formula)
c  the wind: |W_T| and its envelope (Hilbert) — gusts and near-calms; rms against √(2+γ_E−ln4π)
d  one factor's step in W units, ln p/√p, against the rms wind: it peaks at p = e² (= 7.389, where the engine's Mingling sits) at 2/e = 0.736
e  fraction of u-time the wind is at least that factor's step, against p
f  how many zero tones are needed to resolve the jump of a prime p (fraction of the true jump ln p recovered)
Code drawn: rhp.wave, rhp.diverge, rhp.heartbeat.   Run: python -m rhp.atlas_riemann_fermat.wave
"""
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import wave as Wv, diverge as D


def build():
    Z = D.load_zeros()
    lam, psi = Wv.lam_table(9_000_000)
    u = np.linspace(6, 15.9, 30000)
    wt, wl = Wv.trees(u, psi)
    st = Wv.wind_stats(wt)
    fig = plt.figure(figsize=(11.5, 7.8))
    ax = fig.add_subplot(2, 3, 1)
    uu = np.linspace(8, 8.4, 3000); a1, a2 = Wv.trees(uu, psi)
    ax.plot(uu, a1, color=SOUTH, lw=1.6, label="Telperion  W_T (primes)")
    ax.plot(uu, -a2, color=HOT, lw=.8, label="−Laurelin  −W_L (composites)")
    ax.axhline(0, color=GREY, lw=.5); ax.legend(fontsize=6.5, frameon=False, loc="lower left"); ax.set_xlabel("u = ln x")
    d, b = Wv.antisymmetry_defect(u, psi)
    ax.set_title(f"a  two sets, exact mirror (|W_T+W_L| ≤ {d:.3f} ≤ 2/√x)", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 3, 2)
    uu = np.linspace(8, 9.2, 2400); w2, _ = Wv.trees(uu, psi)
    ax.plot(uu, w2, color=INK, lw=1.1, label="from the primes")
    for K, c in ((10, GREY), (100, NORTH), (3066, HOT)):
        ax.plot(uu, Wv.heartbeat(uu, Z[:K]) - np.log(2 * np.pi) * np.exp(-uu / 2), color=c, lw=.8, label=f"{K} zero tones")
    ax.legend(fontsize=6.5, frameon=False, loc="lower left"); ax.set_xlabel("u = ln x")
    c, r = Wv.reconstruction_error(u, psi, Z)
    ax.set_title(f"b  zeros rebuild it (corr {c:.3f}, 3066 tones)", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 3, 3)
    hs = Wv.heartbeat(u, Z[:500]); env = Wv.envelope(hs)
    ax.plot(u, np.abs(hs), color=GREY, lw=.4); ax.plot(u, env, color=HOT, lw=1.1, label="envelope (500 tones)")
    ax.axhline(st["rms"], color=INK, ls="--", lw=.9, label=f"rms {st['rms']:.3f} (√Σ1/|ρ|² = {Wv.MEAN_POWER ** .5:.3f})")
    ax.legend(fontsize=6.5, frameon=False, loc="upper left"); ax.set_xlabel("u = ln x")
    ax.set_title(f"c  the wind: gust factor {st['gust']:.2f}, envelope sd/mean {env.std() / env.mean():.2f}", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 3, 4)
    ps = np.logspace(0.3, 6, 400)
    ax.semilogx(ps, np.log(ps) / np.sqrt(ps), color=INK, lw=1.2, label="ln p/√p (one factor's step)")
    ax.axhline(st["rms"], color=HOT, ls="--", lw=.9, label="rms wind")
    ax.axvline(np.e ** 2, color=NORTH, ls=":", lw=1); ax.text(np.e ** 2 * 1.2, .05, "p = e²", fontsize=7, color=NORTH)
    ax.legend(fontsize=6.5, frameon=False); ax.set_xlabel("p"); ax.set_title("d  a factor's step peaks at p = e², 2/e", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 3, 5)
    fr = [Wv.exceedance(wt, np.log(p) / np.sqrt(p)) for p in ps]
    ax.semilogx(ps, fr, color=SOUTH, lw=1.3); ax.set_xlabel("p"); ax.set_ylabel("fraction of u with |W| ≥ step")
    ax.set_title("e  how often the wind equals a factor", loc="left", fontsize=8.5)
    ax = fig.add_subplot(2, 3, 6)
    Ks = [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3066]
    for p, c in ((11, GREY), (101, NORTH), (307, HOT)):
        ax.semilogx(Ks, [Wv.jump_fraction(p, Z, K) for K in Ks], "o-", color=c, ms=3, lw=1, label=f"p = {p}")
    ax.axhline(1, color=INK, lw=.5); ax.set_ylim(-.5, 1.6); ax.legend(fontsize=6.5, frameon=False); ax.set_xlabel("zero tones K")
    ax.set_title("f  tones needed to resolve a prime's jump", loc="left", fontsize=8.5)
    fig.tight_layout()
    return save(fig, "wave.svg")


if __name__ == "__main__":
    print(build())
