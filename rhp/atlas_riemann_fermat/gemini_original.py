# SPDX-License-Identifier: GPL-3.0-only
"""gemini_original.svg — the ORIGINAL event_horizon_sim.py (Gemini, 143 lines) phases 1-3, and what they are as spectra.

Row 1: phase 2 (the worm, metric swap) at cycle time t = 0.1, 0.3, 0.6, 1.0.   Row 2: phase 3 (the rosette) at t = 0.2, 0.524 (k = 1/3), 0.785 (k = 1/2), 1.0;
col 5 of row 1: line spectrum of the worm's unit phasor against Bessel |J_n(5)|;  col 5 of row 2: rosette sidebands (carrier 1, 1 ± k at 0.1).
Code drawn: rhp.geminisim.   Run: python -m rhp.atlas_riemann_fermat.gemini_original
"""
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import geminisim as G


def build():
    fig, axs = plt.subplots(2, 5, figsize=(15, 6.2), gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 1.1]))
    for ax, t in zip(axs[0, :4], (0.1, 0.3, 0.6, 1.0)):
        th, r = G.original_phase(1, t); ax.plot(r * np.cos(th), r * np.sin(th), lw=.5, color=SOUTH)
        ax.set_title(f"phase 2 (worm)  t = {t}", fontsize=7.5, loc="left")
    for ax, t in zip(axs[1, :4], (0.2, np.pi / 6, np.pi / 4, 1.0)):
        th, r = G.original_phase(2, t); ax.plot(r * np.cos(th), r * np.sin(th), lw=.3, color=HOT)
        ax.set_title(f"phase 3 rosette  t={t:.2f}, k={2 / np.pi * t:.3f}", fontsize=7.5, loc="left")
    for a in axs[:, :4].flat: a.set_aspect("equal", "datalim")
    sb = G.worm_sidebands(0.6); ax = axs[0, 4]
    ax.bar([k - .17 for k, j, m in sb], [j for k, j, m in sb], .34, color=GREY, label="|J_n(5)|")
    ax.bar([k + .17 for k, j, m in sb], [m for k, j, m in sb], .34, color=SOUTH, label="measured line at 1 + n·t")
    ax.set_title("worm = phase modulation, β = 5", fontsize=7.5, loc="left"); ax.legend(fontsize=6); ax.set_xlabel("sideband n", fontsize=7)
    ax = axs[1, 4]
    for i, k in enumerate((1 / 3, 1 / 2, 2 / np.pi)):
        c, up, dn = G.rosette_sidebands(k); ax.bar([-1 + i * .0 + (i - 1) * .25, 0 + (i - 1) * .25, 1 + (i - 1) * .25], [dn, c, up], .24, label=f"k = {k:.3f}")
    ax.set_xticks([-1, 0, 1]); ax.set_xticklabels(["1 − k", "carrier 1", "1 + k"], fontsize=7)
    ax.set_title("rosette = amplitude modulation, depth 0.2", fontsize=7.5, loc="left"); ax.legend(fontsize=6)
    fig.tight_layout()
    return save(fig, "gemini_original.svg")


if __name__ == "__main__":
    print(build())
