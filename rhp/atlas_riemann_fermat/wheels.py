# SPDX-License-Identifier: GPL-3.0-only
"""wheels.svg (static) and wheels_animated.svg (SMIL) — the wheel spring in motion, and how far the primes' grip reaches.

a  error of K prime wheels against |ζ|, by latitude σ: far north the primes win; at the equator they stop helping after ~10² wheels.
b  K*(σ): the number of wheels that does best, against latitude.
c  the first zero emerging from the wheels: ln|P_K| against ln|ζ| near γ₁ = 14.1347.
animated: the first six prime wheels (radius p^{−½}, rate ln p), tip-to-tail, with the trace of their sum.
Code drawn:  rhp.wrap (log_abs_euler, euler_grip).   Run: python -m rhp.atlas_riemann_fermat.wheels   (≈ 1 min)
"""
import math
import numpy as np
from mpmath import zeta, mpc
from ._common import plt, save, HERE, INK, NORTH, SOUTH, HOT, GREY
from .. import wrap as W

KS = [3, 10, 30, 100, 300, 1000, 3000, 10000]


def static():
    t = np.linspace(10, 60, 101)
    sig = [0.5, 0.6, 0.75, 1.0, 1.5]
    G = {s: W.euler_grip(s, t, KS) for s in sig}
    fig = plt.figure(figsize=(11.5, 4.3))
    ax = fig.add_subplot(1, 3, 1)
    for s, col in zip(sig, [HOT, "#b5573a", SOUTH, NORTH, INK]):
        ax.loglog(KS, G[s], "o-", color=col, ms=3, lw=1.3, label=f"σ = {s}")
    ax.set_title("a  error of K wheels vs |ζ|, t ∈ [10,60]", loc="left", fontsize=9)
    ax.set_xlabel("K  (prime wheels)"); ax.set_ylabel("median | ln|ζ| − ln|P_K| |"); ax.legend(fontsize=7, frameon=False)
    ax = fig.add_subplot(1, 3, 2)
    sg = [0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5]
    Ks = [KS[int(np.argmin(W.euler_grip(s, t, KS)))] for s in sg]
    ax.semilogy(sg, Ks, "o-", color=INK)
    ax.axvline(.5, color=HOT, lw=1.2); ax.text(.52, 1.3e4, "equator", color=HOT, fontsize=7); ax.set_ylim(60, 3e4)
    ax.set_title("b  the K that does best, by latitude", loc="left", fontsize=9)
    ax.set_xlabel("σ"); ax.set_ylabel("K*  (grid 3 … 10⁴)")
    ax = fig.add_subplot(1, 3, 3)
    tt = np.linspace(13.0, 15.4, 241)
    true = np.array([math.log(abs(complex(zeta(mpc(.5, x))))) for x in tt])
    L = W.log_abs_euler(0.5, tt, [10, 100, 1000])
    ax.plot(tt, true, color=INK, lw=1.8, label="ln|ζ(½+it)|")
    for row, K, col in zip(L, [10, 100, 1000], [GREY, SOUTH, NORTH]):
        ax.plot(tt, row, color=col, lw=1, label=f"K = {K} wheels")
    ax.axvline(14.134725, color=HOT, lw=.8, ls="--"); ax.set_ylim(-5, 2.5)
    ax.set_title("c  the first zero out of the wheels", loc="left", fontsize=9)
    ax.set_xlabel("t"); ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("How far the primes' grip reaches (CALCULATED: rhp/wrap.py).  At the equator more wheels stop helping.",
                 fontsize=9, x=.01, ha="left")
    fig.tight_layout()
    return save(fig, "wheels.svg"), dict(zip(sg, Ks))


def animated(K=6, speed=1.5, R=78.0, T=40.0):
    """SMIL epicycle machine: wheel p has radius R·p^{−½}, turns clockwise at ln p rad per unit height; one unit of height = 1/speed s."""
    P = W.primes_to(50)[:K]
    cx, cy = 400, 400
    # the trace of Σ p^{-1/2} e^{-it ln p}  (y flipped so the picture is the usual y-up plane)
    t = np.linspace(0, T, 2400)
    z = sum(R * p ** -0.5 * np.exp(-1j * t * math.log(p)) for p in P)
    pts = " ".join(f"{cx + a.real:.1f},{cy - a.imag:.1f}" for a in z)
    arms, close = [], ""
    for p in P:
        r = R * p ** -0.5
        dur = 2 * math.pi / (math.log(p) * speed)
        arms.append(f'<g><animateTransform attributeName="transform" type="rotate" from="0" to="360" dur="{dur:.3f}s" repeatCount="indefinite"/>'
                    f'<circle r="{r:.1f}" class="wh"/><line x1="0" y1="0" x2="{r:.1f}" y2="0" class="arm"/><circle cx="{r:.1f}" r="2.6" class="tip"/>'
                    f'<g transform="translate({r:.1f},0)">')
        close += "</g></g>"
    machine = "".join(arms) + close
    labels = "".join(f'<text x="14" y="{44 + 15 * i}" class="lab">p = {int(p)}: radius p^(−½) = {p ** -0.5:.3f}, rate ln p = {math.log(p):.3f} rad per unit height</text>' for i, p in enumerate(P))
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 740" width="800" height="740">
<style>.bg{{fill:#fff}}.wh{{fill:none;stroke:#b9b9c2;stroke-width:.6}}.arm{{stroke:#1b1b1f;stroke-width:1}}.tip{{fill:#c0392b}}
.tr{{fill:none;stroke:#2a6fb0;stroke-width:.9;opacity:.8}}text{{font-family:ui-monospace,Menlo,monospace;fill:#1b1b1f}}.lab{{font-size:10.5px}}.ttl{{font-size:13px;font-weight:600}}
@media (prefers-color-scheme: dark){{.bg{{fill:#15151a}}.arm{{stroke:#e8e8ee}}text{{fill:#e8e8ee}}.wh{{stroke:#4a4a55}}}}</style>
<rect width="800" height="740" class="bg"/>
<text x="14" y="22" class="ttl">The first {K} prime wheels of the equator, tip to tail: the sum is Σ p^(−½−it), trace to t = {T:g}</text>
{labels}
<polyline points="{pts}" class="tr"/>
<g transform="translate({cx},{cy})">{machine}</g>
<text x="14" y="722" class="lab">clockwise = e^(−it ln p); the mirror wheels p^(−½+it) turn the other way: together, a standing wave.</text>
</svg>'''
    path = HERE / "wheels_animated.svg"
    path.write_text(svg, encoding="utf-8")
    return path


if __name__ == "__main__":
    print(static()); print(animated())
