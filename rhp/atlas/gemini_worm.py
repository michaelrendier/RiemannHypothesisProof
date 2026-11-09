# SPDX-License-Identifier: GPL-3.0-only
"""gemini_worm.svg — the Gemini event-horizon sims' curves (re-drawn from their formulas) beside the ζ pathway.

Row 1: the oldest sim's φ-flow z ↦ 1 + 1/z (13 points, polar, θ ∈ [0, π]); Revised_2 Inertial; Metric Swap; the Transition sim (path crosses r = 1 and inverts).
Row 2: Revised_2 Heartbeat at t = 1, 3, 7 (600 points, polyline through consecutive indices — rosettes whose arm count is ≈ π²/t); ζ path log P_K at s = ½ + iγ₁ for K = 65 and K = 2000.
Code drawn: rhp.geminisim, rhp.primepath.   Run: python -m rhp.atlas.gemini_worm
"""
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import geminisim as G, primepath as PP


def build():
    fig, axs = plt.subplots(2, 4, figsize=(13, 6.6))
    z = G.phi_flow(); th = np.linspace(0, np.pi, len(z)); ax = axs[0, 0]
    ax.plot(z * np.cos(th), z * np.sin(th), ".--", color=HOT, lw=.8)
    ax.set_title(f"a  oldest sim: z ↦ 1 + 1/z  (→ φ, err {G.fib_ratio_error():.1e})", fontsize=7.5, loc="left")
    for ax, (ph, t, lab) in zip((axs[0, 1], axs[0, 2]), (("Inertial", 1, "b  Revised_2 Inertial, t = 1"), ("Metric Swap", 1, "c  Revised_2 Metric Swap, t = 1"))):
        th, r = G.revised2(ph, t); ax.plot(r * np.cos(th), r * np.sin(th), "-", lw=.35, color=INK); ax.set_title(lab, fontsize=7.5, loc="left")
    (to, ro), (ti, ri) = G.transition(); ax = axs[0, 3]
    ax.plot(to * np.cos(to), ro * np.sin(to), color=NORTH, lw=1.2, label="inertial r ≥ 1")
    ax.plot(ti * np.cos(ti), ri * np.sin(ti), color=SOUTH, lw=1.2, label="inverted r → (1+gap)/r")
    ax.add_patch(plt.Circle((0, 0), 1, fill=False, ls="--", color=GREY)); ax.legend(fontsize=6)
    ax.set_title("d  Transition sim: the path crosses r = 1", fontsize=7.5, loc="left")
    for ax, t, lab in zip(axs[1, :3], (1, 3, 7), "efg"):
        x, y = G.heartbeat_xy(t); ax.plot(x, y, "-", lw=.3, color=INK)
        ax.set_title(f"{lab}  Heartbeat t = {t}  ({G.turns_per_index(t):.3f} turn/step)", fontsize=7.5, loc="left")
    ax = axs[1, 3]; P = PP.first_primes(2000); L = PP.log_product_path(PP.GAMMA1, P)
    ax.plot(L.real[:65], L.imag[:65], color=HOT, lw=.9); ax.plot(L.real, L.imag, color=NORTH, lw=.35)
    ax.set_title("h  ζ path log P_K, s = ½ + iγ₁; red = first 65 primes, K = 2000", fontsize=7.5, loc="left")
    for a in axs.flat: a.set_aspect("equal", "datalim")
    fig.tight_layout()
    return save(fig, "gemini_worm.svg")


if __name__ == "__main__":
    print(build())
