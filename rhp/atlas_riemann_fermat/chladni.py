# SPDX-License-Identifier: GPL-3.0-only
"""chladni.svg — the two hemispheres as one oscillation: parities, the flux ('motion'), the helix of nodal lines, the lobe.

a  nodal lines on the sphere in the chart (latitude π(σ−½), longitude θ(t) mod 2π): Re ξ = 0 (red) and Im ξ = 0 (blue) are meridian arcs, a ladder wrapped around the axis
b  parity at t = 30: Re ξ even (free, ∂σ = 0 at the equator), Im ξ odd (held, = 0 there): two polarisations with different boundary conditions
c  the flux, scale-free ∂σ ln|ξ| = J/|ξ|² (J = ½∂σ|ξ|² = Im(conj ξ ∂tξ)): zero on the equator, opposite in the two hemispheres, positive everywhere in the north in this window
d  inject an off-line zero pair (δ = 0.08 at γ = 25) and the flux reverses in a pocket: the footprint of an off-line zero (a test of the detector, not a search)
e  the lobe: angle swept by the Hardy phase between consecutive zeros: mean π, spread = the GUE zero-spacing distribution × π
Code drawn: rhp.chladni, rhp.nodal.   Run: python -m rhp.atlas_riemann_fermat.chladni   (≈ 1 min)
"""
import numpy as np
from mpmath import mpc
from ._common import plt, save, view, INK, NORTH, SOUTH, HOT, GREY
from .. import chladni as C, nodal as N


def build():
    zeros = np.loadtxt(__file__.replace("rhp/atlas_riemann_fermat/chladni.py", "data/zeros_1_to_3600.csv"))
    fig = plt.figure(figsize=(11.5, 8.2))
    # a
    ax = fig.add_subplot(2, 3, 1)
    sig = np.linspace(0.0, 1.0, 51); t = np.linspace(0.5, 50, 330)
    Z = N.grid(sig, t)
    segs = C.nodal_segments(sig, t, Z)
    th = np.linspace(0, 2 * np.pi, 200); ax.plot(np.cos(th), np.sin(th), color=GREY, lw=.5)
    for name, col in (("re", HOT), ("im", NORTH)):
        for sg in segs[name]:
            if name == "im" and np.all(np.abs(sg[:, 0] - .5) < 0.004): continue
            P = C.chart_xyz(sg[:, 0], sg[:, 1])
            x, up, depth = view(P, elev=22, azim=-50)
            front = depth > 0
            ax.plot(np.where(front, x, np.nan), np.where(front, up, np.nan), color=col, lw=.7)
            ax.plot(np.where(~front, x, np.nan), np.where(~front, up, np.nan), color=col, lw=.4, alpha=.25)
    eq = np.column_stack([np.cos(th), np.sin(th), 0 * th]); x, up, depth = view(eq, elev=22, azim=-50)
    ax.plot(np.where(depth > 0, x, np.nan), np.where(depth > 0, up, np.nan), color=INK, lw=1.2)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title("a  nodal lines wrapped by the Hardy phase", loc="left", fontsize=8.5)
    # b
    ax = fig.add_subplot(2, 3, 2)
    ss = np.linspace(0.02, 0.98, 97)
    v = np.array([complex(N.xi(mpc(s, 30.0))) for s in ss])
    sc = np.abs(v).max()
    ax.plot(ss, v.real / sc, color=HOT, label="Re ξ  (even, free)"); ax.plot(ss, v.imag / sc, color=NORTH, label="Im ξ  (odd, held)")
    ax.axvline(.5, color=INK, lw=.7, ls=":"); ax.axhline(0, color=GREY, lw=.5)
    wre, wim = C.parity_defect([0.05, 0.2, 0.4], [14.0, 25.0, 40.0])
    ax.text(.02, .03, f"parity defects: Re {wre:.0e}, Im {wim:.0e}", transform=ax.transAxes, fontsize=7)
    ax.set_xlabel("σ  at t = 30"); ax.legend(fontsize=7, frameon=False); ax.set_title("b  two polarisations at the fret", loc="left", fontsize=8.5)
    # c
    ax = fig.add_subplot(2, 3, 3)
    sg2 = np.linspace(0.02, 0.98, 49); t2 = np.linspace(12, 40, 100)
    J = C.rel_flux_grid(sg2, t2)
    im = ax.imshow(np.sign(J) * np.log10(1 + np.abs(J)), origin="lower", aspect="auto", extent=[0.02, 0.98, 12, 40], cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.axvline(.5, color="k", lw=.8)
    v_ = C.north_violations(J, sg2)
    ax.text(.03, .03, f"north J ≤ 0: {v_} points; |J| on equator ≤ {np.abs(J[:, 24]).max():.0e}", transform=ax.transAxes, fontsize=6.5, color="k")
    fig.colorbar(im, ax=ax, shrink=.75, label="signed log J"); ax.set_xlabel("σ"); ax.set_ylabel("t"); ax.set_title("c  flux: nothing moves on the equator", loc="left", fontsize=8.5)
    # d
    ax = fig.add_subplot(2, 3, 4)
    f = C.with_offline_quadruple(0.08, 25.0)
    sg3 = np.linspace(0.5, 0.72, 23)[1:]; t3 = np.linspace(23, 27, 61)
    J3 = C.rel_flux_grid(sg3, t3, f)
    ax.imshow(np.sign(J3) * np.log10(1 + np.abs(J3)), origin="lower", aspect="auto", extent=[sg3[0], sg3[-1], 23, 27], cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.contour(sg3, t3, J3, levels=[0], colors="k", linewidths=.8)
    ax.set_title("d  injected δ = 0.08 at γ = 25: flux reverses", loc="left", fontsize=8.5); ax.set_xlabel("σ"); ax.set_ylabel("t")
    ax.text(.03, .03, f"J ≤ 0 at {int((J3 <= 0).sum())} of {J3.size} points", transform=ax.transAxes, fontsize=7, color="k")
    # e
    ax = fig.add_subplot(2, 3, 5)
    la = C.lobe_angles(zeros[:400])
    ax.hist(la, bins=50, color=NORTH, alpha=.65); ax.axvline(np.pi, color=HOT, lw=1.2)
    ax.set_title(f"e  angle per lobe: mean {la.mean():.3f}, π = 3.142", loc="left", fontsize=8.5); ax.set_xlabel("θ(γ_{k+1}) − θ(γ_k)")
    ax = fig.add_subplot(2, 3, 6); ax.axis("off")
    ax.text(0, 1, "Reading\n\n"
            "• ξ(1−σ+it) = conj ξ(σ+it): the mirror is complex\n  conjugation. Re is even, Im is odd.\n"
            "• 'Moves' defined as the flux J = Im(conj ξ ∂tξ)\n  = ½∂σ|ξ|²: zero on the equator, opposite in the\n  hemispheres. The equator is where nothing flows.\n"
            "• Each rung is a meridian arc; rungs step in longitude\n  by the lobe angle: mean π, spread = GUE spacing × π.\n"
            "• An off-line zero would reverse J in a pocket\n  (d). In this window the north has no reversal.\n\n"
            "READING, not a theorem: the Chladni sand is the\nnodal set, the 'plate' is the analytic ξ.",
            va="top", fontsize=7.4, family="monospace")
    fig.tight_layout()
    return save(fig, "chladni.svg")


if __name__ == "__main__":
    print(build())
