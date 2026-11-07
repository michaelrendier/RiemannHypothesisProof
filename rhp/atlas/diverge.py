# SPDX-License-Identifier: GPL-3.0-only
"""diverge.svg — the divergence detector: a growing mode (a DILATE) among rotations.

a  CONTROL on the real ψ (x ≤ 10⁷): slope of ln A against u at the first 20 true zeros; scatter about 0 is the false-alarm floor (±3 sd band).
b  an injected off-line quadruple at γ₁ for δ = 0, 0.1, 0.2, 0.3: the lock-in amplitude along u; the mirror partner makes it a cosh, not a line.
c  recovered slope against δ (dots) with the 3-sd floor, and the smallest δ the detector can see at reach u_max = 11, 13, 16 (labels).
d  frequency scan: amplitude against γ for the real ψ (the true zeros are the peaks) and with δ = 0.3 injected at γ = 17.5 (not a zero): it is found there.
Code drawn: rhp.diverge.   Run: python -m rhp.atlas.diverge
"""
import numpy as np
from ._common import plt, save, INK, NORTH, SOUTH, HOT, GREY
from .. import diverge as D


def build():
    psi = D.psi_table(10 ** 7)
    fig = plt.figure(figsize=(11.5, 7.6))
    # a
    ax = fig.add_subplot(2, 3, 1)
    s, m, sd = D.control(psi, 20)
    g = D.load_zeros(20)
    ax.axhspan(-3 * sd, 3 * sd, color=GREY, alpha=.25, label=f"±3 sd = ±{3 * sd:.3f}")
    ax.plot(g, s, "o", color=INK, ms=4); ax.axhline(0, color=GREY, lw=.5)
    ax.set_xlabel("zero height γ"); ax.set_ylabel("d ln A / du"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("a  control: real ψ, true zeros", loc="left", fontsize=9)
    # b
    ax = fig.add_subplot(2, 3, 2)
    u = np.linspace(6, 16, 60000); E = D.residual(u, psi)
    for d, c in [(0, INK), (0.1, NORTH), (0.2, SOUTH), (0.3, HOT)]:
        Ed = E + (D.quadruple(u, d, 14.134725) if d > 0 else 0)
        cs, A = D.amplitude_track(u, Ed, 14.134725, 6, 16, 3.0, 0.5)
        ax.semilogy(cs, A, "o-", ms=3, color=c, label=f"δ = {d}")
    ax.set_xlabel("u = ln x"); ax.set_ylabel("lock-in amplitude at γ₁"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("b  injected growth at γ₁", loc="left", fontsize=9)
    # c
    ax = fig.add_subplot(2, 3, 3)
    ds = [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
    rec = D.inject_and_recover(psi, ds)
    ax.plot([d for d, _ in rec], [r for _, r in rec], "o-", color=INK, label="recovered slope")
    ax.plot([0, 0.5], [0, 0.5], color=GREY, ls=":", label="slope = δ")
    ax.axhline(3 * sd, color=HOT, lw=1, ls="--", label="3-sd floor (u ≤ 16)")
    txt = []
    for u1 in (11, 13, 16):
        dm, _ = D.detection_limit(psi, u1); txt.append(f"u_max {u1}:  δ_min ≈ {dm:.2f}")
    ax.text(0.02, 0.36, "\n".join(txt), fontsize=7.5, family="monospace")
    ax.set_xlabel("injected δ"); ax.set_ylabel("slope"); ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.set_title("c  what it can see", loc="left", fontsize=9)
    # d
    ax = fig.add_subplot(2, 3, (4, 5))
    gam = np.arange(10, 40, 0.1)
    _, s0, a0 = D.scan(psi, gam)
    _, s1, a1 = D.scan(psi, gam, extra=lambda uu: D.quadruple(uu, 0.3, 17.5))
    ax.semilogy(gam, a0, color=INK, lw=1, label="real ψ (peaks = true zeros)")
    ax.semilogy(gam, a1, color=HOT, lw=1, label="+ δ = 0.3 quadruple at γ = 17.5")
    for z in D.load_zeros(8): ax.axvline(z, color=GREY, lw=.4, ls=":")
    pk = D.candidate_peaks(gam, s1, a1, 0.1, 0.06)
    for (gp, a, sl) in pk: ax.plot([gp], [a], "v", color=HOT, ms=6)
    ax.set_xlabel("frequency γ"); ax.set_ylabel("mean lock-in amplitude"); ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.set_title("d  frequency scan: found at 17.5, not a zero", loc="left", fontsize=9)
    ax = fig.add_subplot(2, 3, 6); ax.axis("off")
    ax.text(0, 1, "Reading\n\n"
            "• On the line every zero is a pure rotation in\n  u (gain 1). Off the line a zero adds a growth\n  e^{δu}: a DILATE, graded by δ.\n"
            "• Its mirror partner decays as e^{−δu}, so the\n  amplitude is a cosh: slope ~ δ²u when small.\n"
            "• Control first: on the real ψ the slopes sit\n  at 0 ± 0.033 (u ≤ 16).\n"
            "• Reach: δ_min ≈ 0.13 at x ≤ 10⁷. Real zeros\n  are on the line to far beyond this; the detector\n  is machinery, and its null result is a\n  statement about this window only.\n\n"
            "Lineage: tier 1 (ln of an amplitude = ADD∘\nSCALE); emergence signature 1 (§5) if it fires.",
            va="top", fontsize=7.6, family="monospace")
    fig.tight_layout()
    return save(fig, "diverge.svg")


if __name__ == "__main__":
    print(build())
