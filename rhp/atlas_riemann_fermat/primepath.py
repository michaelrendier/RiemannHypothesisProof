# SPDX-License-Identifier: GPL-3.0-only
"""primepath_flat.svg, primepath_3axis.svg, primepath_uns.svg — the pathway of the first K primes at s = ½ + it.

flat    log P_K(s) in the complex plane for K = 65, 10⁴, 10⁵, 10⁶ primes; row 1 at t = γ₁ = 14.1347 (first zero), row 2 at t = 3 (off the zeros). First 65 primes drawn in red, the rest in blue.
3axis   (Re log P, Im log P, ln p) for 10⁴ primes at t = γ₁: isometric, top, front, side; |log P_K| against K with the PNT smooth part; the torus radius R_K = √Σ1/p against K.
uns     the partial product P_K on the Riemann sphere (|P| = 1 is the equator): isometric and polar views for 10⁴ primes, and ln|P_K| against K with the equator crossings counted.
Code drawn: rhp.primepath, rhp.unsmap.   Run: python -m rhp.atlas_riemann_fermat.primepath
"""
import numpy as np
from ._common import plt, save, view, INK, NORTH, SOUTH, HOT, GREY
from .. import primepath as PP

KS = (65, 10 ** 4, 10 ** 5, 10 ** 6)


def _data():
    P = PP.first_primes(10 ** 6)
    return P, {t: PP.log_product_path(t, P) for t in (PP.GAMMA1, 3.0)}


def flat(P, L):
    fig, axs = plt.subplots(2, 4, figsize=(12, 6.4))
    for r, t in enumerate((PP.GAMMA1, 3.0)):
        for c, K in enumerate(KS):
            ax = axs[r, c]; z = L[t][:K]; ix = PP.downsample(K)
            ax.plot(z.real[ix], z.imag[ix], color=NORTH, lw=.5)
            ax.plot(z.real[:65], z.imag[:65], color=HOT, lw=1.3)
            ax.plot([0], [0], "k+", ms=5); ax.plot(z.real[-1], z.imag[-1], "o", color=SOUTH, ms=3)
            ax.set_aspect("equal", "datalim")
            ax.set_title(f"{'t = γ₁' if r == 0 else 't = 3'}   K = {K:,}   |log P| = {abs(z[-1]):.1f}", fontsize=7.5, loc="left")
    fig.suptitle("log of the partial Euler product at s = ½ + it: red = first 65 primes, blue = the rest, orange dot = end", fontsize=8.5, x=.01, ha="left")
    fig.tight_layout()
    return save(fig, "primepath_flat.svg")


def three_axis(P, L):
    K = 10 ** 4; z = L[PP.GAMMA1][:K]; tau = np.log(P[:K])
    xyz = np.column_stack([z.real, z.imag, (tau - tau.mean()) * (np.ptp(z.real) / np.ptp(tau)) * .6])
    fig = plt.figure(figsize=(12, 7.2))
    def panel(pos, title, a, b, pts=None):
        ax = fig.add_subplot(2, 3, pos); ax.plot(a[65:], b[65:], color=NORTH, lw=.5); ax.plot(a[:65], b[:65], color=HOT, lw=1.3)
        ax.set_aspect("equal", "datalim"); ax.set_title(title, loc="left", fontsize=8); return ax
    x, up, _ = view(xyz, elev=35.264, azim=45)
    panel(1, "isometric  (Re, Im, ln p)", x, up)
    panel(2, "top  (Re, Im)", xyz[:, 0], xyz[:, 1])
    panel(3, "front  (Re, ln p)", xyz[:, 0], xyz[:, 2])
    panel(4, "side  (Im, ln p)", xyz[:, 1], xyz[:, 2])
    ax = fig.add_subplot(2, 3, 5)
    Kk = np.arange(1, 10 ** 6 + 1); ix = PP.downsample(10 ** 6)
    for t, c, lab in ((PP.GAMMA1, HOT, "t = γ₁"), (3.0, NORTH, "t = 3")):
        ax.loglog(Kk[ix], np.abs(PP.prime_sum_path(t, P))[ix], color=c, lw=1, label=f"|Σ p^(−s)|  {lab}")
        ax.loglog(Kk[ix], np.abs(PP.smooth_prime_sum(t, P))[ix], color=c, lw=.8, ls="--")
    ax.legend(fontsize=6.5, frameon=False); ax.set_xlabel("number of primes K"); ax.set_title("e  the path keeps growing (dashed = PNT smooth part)", loc="left", fontsize=8)
    ax = fig.add_subplot(2, 3, 6)
    R = np.sqrt(np.cumsum(1 / P)); ax.semilogx(Kk[ix], R[ix], color=INK)
    for k in KS:
        ax.plot([k], [R[k - 1]], "o", color=SOUTH, ms=3); ax.annotate(f"{R[k - 1]:.2f}", (k, R[k - 1]), fontsize=7, xytext=(2, -10), textcoords="offset points")
    ax.set_xlabel("number of wheels K"); ax.set_title("f  sphere radius R_K = √Σ1/p holding the torus", loc="left", fontsize=8)
    fig.tight_layout()
    return save(fig, "primepath_3axis.svg")


def uns(P, L):
    K = 10 ** 4; z = L[PP.GAMMA1][:K]; S = PP.sphere_path(z)
    fig = plt.figure(figsize=(12, 4.6))
    th = np.linspace(0, 2 * np.pi, 400); eq = np.column_stack([np.cos(th), np.sin(th), 0 * th])
    for pos, (el, az, title) in enumerate(((24, -35, "isometric"), (90, 0, "from the north pole"), (0, 0, "side, equator edge-on")), 1):
        ax = fig.add_subplot(1, 4, pos)
        a, b, d = view(np.column_stack([S[:, 0], S[:, 1], S[:, 2]]), el, az); e, f, _ = view(eq, el, az)
        ax.plot(e, f, color=GREY, lw=.8); ax.plot(a[65:], b[65:], color=NORTH, lw=.5); ax.plot(a[:65], b[:65], color=HOT, lw=1.3)
        ax.set_aspect("equal"); ax.set_title(title, loc="left", fontsize=8); ax.axis("off")
    ax = fig.add_subplot(1, 4, 4)
    Kk = np.arange(1, 10 ** 6 + 1); ix = PP.downsample(10 ** 6); zz = L[PP.GAMMA1]
    ax.semilogx(Kk[ix], zz.real[ix], color=INK, lw=.7); ax.axhline(0, color=GREY, lw=.7); ax.axvline(65, color=HOT, ls=":", lw=1)
    n = PP.equator_crossings(zz); ax.set_xlabel("number of primes K"); ax.set_ylabel("ln|P_K|   (0 = equator)")
    ax.set_title(f"|P_K| = 1 crossed {n}× in 10⁶ primes", loc="left", fontsize=8)
    fig.suptitle("partial Euler product on the Riemann sphere (t = γ₁): red = first 65 primes", fontsize=8.5, x=.01, ha="left")
    fig.tight_layout()
    return save(fig, "primepath_uns.svg")


def build():
    P, L = _data()
    return flat(P, L), three_axis(P, L), uns(P, L)


if __name__ == "__main__":
    for p in build():
        print(p)
