# SPDX-License-Identifier: GPL-3.0-only
"""SVG figures for notebook 02: the long pathway (thousands of zeros) and the machinery/flow panels."""
from __future__ import annotations
import numpy as np

from . import bigpath
from .render import Panel, colour, CSS, project


def bands(P, X, Y, T, tmax, nb=40, stride=2, width=0.45, opacity=0.4):
    """Draw X,Y as nb contiguous polylines coloured by t (keeps SVG size bounded)."""
    n = len(X)
    edges = np.linspace(0, n - 1, nb + 1).astype(int)
    for a, b in zip(edges[:-1], edges[1:]):
        idx = np.arange(a, b + 1, stride)
        px = [P.px(X[k], Y[k]) for k in idx]
        pts = " ".join("%.1f,%.1f" % q for q in px)
        P.out.append(f'<polyline points="{pts}" fill="none" stroke="{colour(T[(a + b) // 2] / tmax)}" '
                     f'stroke-width="{width}" stroke-opacity="{opacity}" stroke-linejoin="round"/>')


def _svg(W, H, panels, legend):
    body = "".join("".join(p.out) for p in panels)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">'
            f'<style>{CSS}</style><rect width="{W}" height="{H}" class="bg"/>{body}{legend}</svg>')


def long_pathway_svg(p, zeros, path, floor=-2.5):
    t, z = p["t"], p["zeta"]
    tmax, tmin = float(t[-1]), float(t[0])
    n = len(zeros)
    th = bigpath.theta_fast(zeros)
    S = ((-th / np.pi - 0.5 + 0.5) % 1.0) - 0.5            # (direction − π/2)/π wrapped to (−½, ½]  = S(γₙ) mod 1
    panels = []

    R = float(np.abs(z).max()) * 1.03
    A = Panel(40, 50, 520, 520, (-R, R), (-R, R), f"A  flat, t = {tmin:.0f}…{tmax:.0f}  ({n} zeros)", equal=True)
    A.frame(); A.axes(); bands(A, z.real, z.imag, t, tmax); A.dot(0, 0, 3.5, "fh"); panels.append(A)

    rho = np.clip(p["u"] - floor, 0, None)
    ang = np.angle(z)
    Rm = float(rho.max()) * 1.04
    B = Panel(600, 50, 520, 520, (-Rm, Rm), (-Rm, Rm), "B  UNS end-on: ρ = log₁₀|ζ| + 2.5, angle = arg ζ", equal=True)
    B.frame(); B.axes(); bands(B, rho * np.cos(ang), rho * np.sin(ang), t, tmax); B.dot(0, 0, 3, "fh"); panels.append(B)

    C = Panel(1160, 50, 480, 250, (0, n), (-0.5, 0.5), "C  the hub fan: S(γₙ) mod 1, one dot per zero")
    C.frame(); C.axes(0, 0)
    for k in range(n):
        a = C.px(k + 1, float(S[k])); C.out.append(f'<circle cx="{a[0]:.1f}" cy="{a[1]:.1f}" r="0.9" class="fi" fill-opacity=".55"/>')
    C.text(20, 0.43, "direction at the origin = π/2 + π·S(γₙ) (mod π)", "small"); panels.append(C)

    D = Panel(1160, 340, 480, 230, (-0.5, 0.5), (0, 1), "D  histogram of that fan")
    D.frame()
    h, e = np.histogram(S, bins=40, range=(-0.5, 0.5)); h = h / h.max()
    for i in range(40):
        a, b = D.px(e[i], 0), D.px(e[i + 1], h[i])
        D.out.append(f'<rect x="{a[0]:.1f}" y="{b[1]:.1f}" width="{b[0]-a[0]-.6:.1f}" height="{a[1]-b[1]:.1f}" class="fi" fill-opacity=".6"/>')
    D.text(-0.48, 0.93, f"mean {S.mean():+.3f}, sd {S.std():.3f}", "small"); panels.append(D)

    # E: S(t) = N(t) − θ(t)/π − 1 at the midpoints between consecutive zeros (where N is unambiguous)
    mid = (zeros[:-1] + zeros[1:]) / 2
    St = (np.arange(1, len(mid) + 1) + 1) - (bigpath.theta_fast(mid) / np.pi + 1)     # +1: the zero below t = 20
    E = Panel(40, 640, 520, 330, (tmin, tmax), (-2.0, 2.0), "E  S(t) = N(t) − θ(t)/π − 1, between zeros")
    E.frame(); E.axes(0, 0)
    E.line(mid, St, mid, tmax, 0.7)
    E.text(tmin + 40, 1.75, f"sd {St.std():.2f}, max |S| {np.abs(St).max():.2f}", "small")
    panels.append(E)

    # F: the spring (front pose), height = t
    zs = 1.7 * Rm / (tmax - tmin)
    X, Y, Zh = rho * np.cos(ang), rho * np.sin(ang), (t - tmin) * zs
    F = Panel(600, 640, 520, 330, (-Rm * 1.05, Rm * 1.05), (0, 1.7 * Rm), "F  spring, front: height = t (x = ρ cos φ)", equal=False)
    F.frame(); bands(F, X, Zh, t, tmax, nb=60, stride=3, width=0.3, opacity=0.5); panels.append(F)

    # G: spring, iso, t ≤ 300
    m = t <= 300
    Xi, Yi, Zi = rho[m] * np.cos(ang[m]), rho[m] * np.sin(ang[m]), (t[m] - tmin) * (1.7 * Rm / 280)
    pts = np.column_stack([Xi, Yi, Zi])
    sx, sy, _ = project(pts, (55, 0, 25))
    s = max(np.ptp(sx), np.ptp(sy) * 480 / 330) / 2 * 1.08
    cx, cy = (sx.min() + sx.max()) / 2, (sy.min() + sy.max()) / 2
    G = Panel(1160, 640, 480, 330, (cx - s, cx + s), (cy - s * 330 / 480, cy + s * 330 / 480), "G  spring, iso [55,0,25], t = 20…300")
    G.frame(); bands(G, sx, sy, t[m], 300.0, nb=24, stride=1, width=0.5, opacity=0.6); panels.append(G)

    legend = (f'<text x="40" y="1010" class="small">colour = t (violet {tmin:.0f} → yellow {tmax:.0f}); Riemann–Siegel (first correction), zeros polished with mpmath. '
              f'Orthographic. Panels A, B overlay {int((bigpath.theta_fast([tmax])[0] - bigpath.theta_fast([tmin])[0]) / (2 * np.pi))} full turns.</text>')
    open(path, "w", encoding="utf-8").write(_svg(1680, 1040, panels, legend))
    return path
