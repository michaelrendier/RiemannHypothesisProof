# SPDX-License-Identifier: GPL-3.0-only
"""boundary_uclock.svg — Gemini's constant-angle clock read as u = ln p, against the ζ Euler path.

a  loop boundaries of the ζ path at t = γ₁: loop number against u = ln p (straight line, slope t/2π) and against prime index n (log axis) — the worm's loops are uniform in θ_b ≡ t·u.
b  radius growth factor per loop: ζ at t = 3, 5, 10 (measured, dots) with the smooth-spiral prediction (lines), ζ at t = γ₁ (near a zero: suppressed); the worm (hyperbolic: 1.25, 1.33, 1.49, 1.95, 19).
c  tangent-angle wobble of ζ against the analytic ½p^{-½}: the worm's index is the constant 5.
Code drawn: rhp.boundary.   Run: python -m rhp.atlas_riemann_fermat.boundary_uclock
"""
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import boundary as B, primepath as PP


def build():
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
    L = B.zeta_loops(PP.GAMMA1, 200000); j = np.arange(len(L["u"])); ax = axs[0]
    ax.plot(L["u"], j, "o-", color=NORTH, ms=3, lw=.8, label="loop number vs u = ln p")
    ax.plot(L["u"], L["u"] * PP.GAMMA1 / (2 * np.pi) - (L["u"][0] * PP.GAMMA1 / (2 * np.pi)), "--", color=GREY, lw=.8, label="slope t/2π")
    ax2 = ax.twiny(); ax2.set_xscale("log"); ax2.plot(L["idx"], j, ".", color=HOT, ms=3); ax2.set_xlabel("prime index n (red dots, log axis)", fontsize=7)
    ax.set_xlabel("u = ln p", fontsize=7); ax.legend(fontsize=6, loc="upper left")
    ax.set_title(f"a  ζ loops are uniform in u (spacing {L['spacing_u'].mean():.4f} vs 2π/t = {L['expected_spacing']:.4f})", fontsize=7.5, loc="left", y=1.18)
    ax = axs[1]
    for t, c in ((3.0, NORTH), (5.0, SOUTH), (10.0, HOT)):
        Lt = B.zeta_loops(t, 200000); k = np.arange(1, len(Lt["growth"]) + 1)
        ax.plot(k, Lt["growth"], "o", color=c, ms=3, label=f"ζ t = {t:g}"); ax.plot(k, Lt["smooth_growth"], "-", color=c, lw=.8)
    Lg = B.zeta_loops(PP.GAMMA1, 200000); ax.plot(np.arange(1, len(Lg["growth"]) + 1), Lg["growth"], "x", color=INK, ms=3, label="ζ t = γ₁ (on a zero)")
    ax.plot(np.arange(1, 6), B.worm_loop_growth(), "s-", color="k", ms=4, label="worm (hyperbolic)"); ax.set_yscale("log"); ax.set_xlim(0, 20)
    ax.set_xlabel("loop number", fontsize=7); ax.set_ylabel("radius growth per loop", fontsize=7); ax.legend(fontsize=6)
    ax.set_title("b  radius growth per loop (dots measured, lines smooth law)", fontsize=7.5, loc="left")
    ax = axs[2]; P = PP.first_primes(9592); z = PP.log_product_path(PP.GAMMA1, P, .5); dz = np.diff(np.concatenate([[0], z]))
    resid = np.unwrap(np.angle(dz)) + PP.GAMMA1 * np.log(P); resid -= np.round(np.median(resid[P > 1000]) / (2 * np.pi)) * 2 * np.pi
    m = P > 100; ax.loglog(P[m], np.abs(resid[m] - np.median(resid[m])), ".", ms=1.2, color=NORTH, alpha=.5, label="|step angle + t ln p|")
    ax.loglog(P[m], .5 / np.sqrt(P[m]), color=HOT, lw=1.2, label="½ p^{-½}"); ax.axhline(5, color="k", lw=.8, ls="--", label="worm β = 5")
    ax.set_xlabel("p", fontsize=7); ax.legend(fontsize=6); ax.set_title("c  ζ wobble index decays; worm's is constant", fontsize=7.5, loc="left")
    fig.tight_layout(); return save(fig, "boundary_uclock.svg")


if __name__ == "__main__":
    print(build())
