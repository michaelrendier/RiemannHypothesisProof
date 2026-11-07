# SPDX-License-Identifier: GPL-3.0-only
"""uns_map.svg — flat Cartesian → UNS, with no numbers in it: what shapes become, and how to tell an artefact from a shape.

a  lines σ = c and a circle in the plane (the chart s = σ + it)
b  their images on the sphere (Cayley, then stereographic): all circles; σ = ½ (thick) is the great circle = the equator, the others are small circles
c  a cusp in the plane   d  its image: still a cusp (conformal map), the flat and sphere speeds vanish together
e  a back-and-forth trace (reversals) that is a UNIFORM ROTATION seen edge-on: the Hilbert companion supplies the missing coordinate
f  three spiral laws on test spirals (R² of the law's straight line): Archimedean (pancake), logarithmic (loxodrome), Euler/Cornu
Code drawn: rhp.unsmap.   Run: python -m rhp.atlas.uns_map
"""
import numpy as np
from scipy.special import fresnel
from ._common import plt, save, view, INK, NORTH, SOUTH, HOT, GREY
from .. import unsmap as U


def _sphere(ax, P, color, lw=1.2, label=None, alpha_back=.25):
    x, up, depth = view(P)
    front = depth > 0
    ax.plot(np.where(front, x, np.nan), np.where(front, up, np.nan), color=color, lw=lw, label=label)
    ax.plot(np.where(~front, x, np.nan), np.where(~front, up, np.nan), color=color, lw=lw * .8, alpha=alpha_back)


def _outline(ax):
    th = np.linspace(0, 2 * np.pi, 200); ax.plot(np.cos(th), np.sin(th), color=GREY, lw=.6); ax.set_aspect("equal"); ax.axis("off")


def build():
    fig = plt.figure(figsize=(11.5, 7.8))
    t = np.linspace(-60, 60, 6001)
    cols = {0.2: NORTH, 0.35: "#5a8fb8", 0.5: HOT, 0.7: SOUTH, 0.9: "#a05a1a"}
    ax = fig.add_subplot(2, 3, 1)
    for c, col in cols.items(): ax.plot([c] * 2, [-3, 3], color=col, lw=2 if c == .5 else 1)
    th = np.linspace(0, 2 * np.pi, 300); ax.plot(1.6 + .5 * np.cos(th), .8 * np.sin(th) * 0 + .5 * np.sin(th) + 1.5, color=INK, lw=1)
    ax.set_xlim(-.3, 2.3); ax.set_ylim(-3, 3); ax.set_xlabel("σ"); ax.set_ylabel("t"); ax.set_title("a  lines and a circle in the plane", loc="left", fontsize=9)
    ax = fig.add_subplot(2, 3, 2); _outline(ax)
    for c, col in cols.items(): _sphere(ax, U.lift(c + 1j * t), col, 2 if c == .5 else 1)
    _sphere(ax, U.lift(1.6 + .5 * np.exp(1j * th) + 1.5j), INK, 1)
    ax.set_title("b  circles; σ = ½ = equator", loc="left", fontsize=9)
    # cusp
    c = np.linspace(-1, 1, 2001)
    z = 0.6 + 0.5j + 0.4 * (c ** 2 + 1j * c ** 3)
    ax = fig.add_subplot(2, 3, 3); ax.plot(z.real, z.imag, color=INK); ax.plot([z[1000].real], [z[1000].imag], "o", color=HOT)
    f0, g0 = U.cusp_survives_lift(z, c, 0.0)
    ax.set_title("c  a cusp in the plane", loc="left", fontsize=9); ax.set_xlabel("σ"); ax.set_ylabel("t"); ax.set_aspect("equal")
    ax.text(.02, .02, f"tip speed: flat {f0:.0e}, sphere {g0:.0e}\n(× median): the lift keeps a true cusp", transform=ax.transAxes, fontsize=7)
    ax = fig.add_subplot(2, 3, 4); _outline(ax); _sphere(ax, U.lift(z), INK, 1.4)
    x, up, d = view(U.lift(z[1000:1001])); ax.plot(x, up, "o", color=HOT)
    ax.set_title("d  its image on the sphere", loc="left", fontsize=9)
    # edge-on
    ax = fig.add_subplot(2, 3, 5)
    tt = np.linspace(0, 12, 3000); xx = np.cos(2 * np.pi * 0.6 * tt)
    amp, ph = U.analytic_signal(xx); h = amp * np.sin(ph)
    ax.plot(tt, xx, color=INK, lw=1, label="flat trace: reverses")
    ax.plot(tt, h, color=HOT, lw=.8, alpha=.8, label="Hilbert companion")
    res = U.edge_on_rotation_test(tt, xx)
    ax.set_title("e  a reversal that is a rotation edge-on", loc="left", fontsize=9)
    ax.text(.02, .02, f"amplitude CV {res['amp_cv']:.1e}, phase linear r² {res['phase_linear_r2']:.6f}", transform=ax.transAxes, fontsize=7)
    ax.legend(fontsize=7, frameon=False, loc="upper right"); ax.set_xlabel("t")
    # spirals
    ax = fig.add_subplot(2, 3, 6); ax.axis("off")
    ph_ = np.linspace(0, 12 * np.pi, 3000)
    S, C = fresnel(np.linspace(0.05, 6, 3000))
    rows = [("Archimedean (pancake)", (1 + .3 * ph_) * np.exp(1j * ph_), 0j), ("logarithmic (loxodrome)", np.exp(.12 * ph_) * np.exp(1j * ph_), 0j),
            ("Euler / Cornu", C + 1j * S, C[-1] + 1j * S[-1])]
    txt = "f  which law is it?   r² of each law\n\n" + f"{'':26s}{'arch':>7s}{'log':>7s}{'euler':>7s}\n"
    for name, zz, ctr in rows:
        L = U.spiral_laws(zz, ctr)
        txt += f"{name:26s}{L['archimedean']:7.3f}{L['logarithmic'] if L['logarithmic'] == L['logarithmic'] else float('nan'):7.3f}{L['euler']:7.3f}\n"
    txt += ("\nReading\n• Möbius: lines and circles → circles; the\n  map is conformal, so a true cusp survives.\n"
            "• A reversal in a 1-D trace may be a rotation\n  edge-on: constant amplitude + linear phase.\n"
            "• Number-free: these detectors take any trace.")
    ax.text(0, 1, txt, va="top", fontsize=7.8, family="monospace")
    fig.tight_layout()
    return save(fig, "uns_map.svg")


if __name__ == "__main__":
    print(build())
