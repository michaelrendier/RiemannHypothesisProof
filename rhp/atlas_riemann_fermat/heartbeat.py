# SPDX-License-Identifier: GPL-3.0-only
"""heartbeat.svg — what the standing wave contains under its curve, and ADD / SCALE / SIGN acting on it.

a  the unit circle's area π split: the zero circles (π/|ρ|², all zeros) hold 2.31 %, the Γ/π ground state the rest (ESTABLISHED sum Σ1/ρ = 1+γ_E/2−½ln4π)
b  how the zero-circle area is shared: the first zero alone holds 21.6 %, the first 10 hold 58.6 %, the 34 on the pancake coil 76 %, the first 100 86.5 %
c  mean power of the real wave on the exact ψ, windows of u, against the Parseval value 2+γ_E−ln4π
d  one standing wave as two counter-rotating circles: equal radii, the sum is real; off the line (δ = 0.15) both radii grow as e^{δu}
e  ADD, SCALE, SIGN on a tone: a step in ln x is a spiral similarity (scale e^{Δu/2}, rotation γΔu); conjugation is the reflection
f  reading
Code drawn: rhp.heartbeat, rhp.diverge.   Run: python -m rhp.atlas_riemann_fermat.heartbeat
"""
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import heartbeat as H, diverge as D


def build():
    z = np.loadtxt(__file__.replace("rhp/atlas_riemann_fermat/heartbeat.py", "data/zeros_1_to_3600.csv"))
    fig = plt.figure(figsize=(11.5, 7.8))
    # a
    ax = fig.add_subplot(2, 3, 1)
    ax.barh([0], [H.SUM_INV_RHO], color=HOT, label=f"zero circles: {H.SUM_INV_RHO:.5f} π")
    ax.barh([0], [H.GROUND_SHARE], left=[H.SUM_INV_RHO], color=NORTH, alpha=.6, label=f"Γ/π ground state: {H.GROUND_SHARE:.5f} π")
    ax.set_xlim(0, 1); ax.set_yticks([]); ax.set_xlabel("fraction of the unit circle's area π"); ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(.5, 1.0))
    ax.set_title("a  area budget of the unit circle", loc="left", fontsize=9)
    # b
    ax = fig.add_subplot(2, 3, 2)
    cum = np.cumsum(H.tone_radius(z) ** 2) / H.SUM_INV_RHO
    k = np.arange(1, len(z) + 1)
    ax.semilogx(k, H.SUM_INV_RHO * cum, color=INK, lw=1.4, label="partial sum Σ_{k≤K} 1/|ρ|²")
    ax.axhline(H.SUM_INV_RHO, color=HOT, ls="--", lw=1, label="exact 1+γ_E/2−½ln4π")
    for kk in (1, 10, 34, 100):
        ax.plot([kk], [H.SUM_INV_RHO * cum[kk - 1]], "o", color=SOUTH, ms=4); ax.annotate(f"{cum[kk - 1] * 100:.1f} %", (kk, H.SUM_INV_RHO * cum[kk - 1]), fontsize=6.5, xytext=(2, -9), textcoords="offset points")
    s, st, ex = H.budget(z)
    ax.text(.03, .62, f"3066 zeros: {s:.6f}\n+ tail estimate: {st:.8f}\nexact: {ex:.8f}", transform=ax.transAxes, fontsize=6.8)
    ax.set_xlabel("number of zeros K"); ax.legend(fontsize=6.5, frameon=False, loc="lower right"); ax.set_title("b  who holds the zero-circle area", loc="left", fontsize=9)
    # c
    ax = fig.add_subplot(2, 3, 3)
    psi = D.psi_table(10 ** 7)
    wins = [(6, 10), (6, 12), (8, 12), (8, 14), (10, 14), (10, 16.1), (6, 16.1)]
    vals = [H.mean_power_measured(psi, a, b)[0] for a, b in wins]
    ax.plot(range(len(wins)), vals, "o", color=INK); ax.axhline(H.MEAN_POWER, color=HOT, ls="--", lw=1.2)
    ax.set_xticks(range(len(wins))); ax.set_xticklabels([f"{a:g}–{b:g}" for a, b in wins], rotation=40, fontsize=6.5)
    ax.text(.03, .05, f"Parseval 2+γ_E−ln4π = {H.MEAN_POWER:.5f}", transform=ax.transAxes, fontsize=7, color=HOT)
    ax.set_xlabel("window in u = ln x"); ax.set_ylabel("⟨E²⟩"); ax.set_title("c  mean power, real ψ", loc="left", fontsize=9)
    # d
    ax = fig.add_subplot(2, 3, 4)
    u = np.linspace(0, 4 * 2 * np.pi / z[0], 2500)
    for dlt, col, ls in [(0.0, INK, "-"), (0.15, SOUTH, "--")]:
        f, m, s_ = H.counter_rotating_pair(u, z[0], dlt)
        ax.plot(f.real, f.imag, color=col, lw=.9, ls=ls); ax.plot(m.real, m.imag, color=col, lw=.5, ls=ls, alpha=.6)
    ax.plot([0], [0], "k+"); ax.set_aspect("equal")
    ax.set_title("d  two circles (γ₁); dashed δ = 0.15", loc="left", fontsize=8.5)
    ax.set_xlabel("Re"); ax.set_ylabel("Im")
    # e
    ax = fig.add_subplot(2, 3, 5)
    u = np.linspace(0, 3.0, 800); f = H.tone(u, z[0]) / np.exp(0)
    ax.plot(f.real, f.imag, color=INK, lw=1); ax.plot(f.real, -f.imag, color=GREY, lw=.8, ls="--")
    for du in (0, 1, 2, 3):
        p = H.tone(np.array([du]), z[0])[0]; ax.plot([p.real], [p.imag], "o", color=HOT, ms=4); ax.plot([p.real], [-p.imag], "o", color=GREY, ms=3)
    sc, rot = H.spiral_similarity(z[0], 1.0)
    ax.text(.03, .03, f"Δu = 1: scale {sc:.4f} = e^{{1/2}}, rotation {rot:.4f} rad\n(γ₁·Δu = {z[0]:.3f} ≡ {z[0] % (2 * np.pi):.3f} mod 2π); dashed: SIGN", transform=ax.transAxes, fontsize=6.5)
    ax.set_aspect("equal"); ax.set_title("e  ADD = spiral similarity", loc="left", fontsize=8)
    ax = fig.add_subplot(2, 3, 6); ax.axis("off")
    ax.text(0, 1, "Reading\n\n"
            "• One zero pair = one standing wave = two circles\n  of radius 1/|ρ| turning opposite ways. Area π/|ρ|².\n"
            "• Σ over all zeros of 1/|ρ|² = 0.0230957 (Hadamard\n  constant −B): the zeros fill 2.31 % of π.\n"
            "• First zero alone: 21.6 % of that; first 10: 58.6 %.\n"
            "• Mean power = 2+γ_E−ln4π = 0.04619; measured\n  0.0453–0.0459 on exact ψ (u ≤ 16, γ ≤ 3066 only).\n"
            "• ADD (step in ln x) = scale e^{Δu/2} ∘ rotation γΔu.\n  SCALE is the e^{u/2} flare. SIGN is conjugation.\n"
            "• Unit circle of area π: only as the normalisation;\n  no single tone holds it (largest holds 0.5 %).",
            va="top", fontsize=7.4, family="monospace")
    fig.tight_layout()
    return save(fig, "heartbeat.svg")


if __name__ == "__main__":
    print(build())
