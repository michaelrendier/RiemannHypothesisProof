# SPDX-License-Identifier: GPL-3.0-only
"""Dependency-free SVG renderer for the flattening map, plus the CSV a Blender import reads.

Orthographic (scad-spatial: orthographic to measure). Framing of the 3D spring is gated by
rhp.camera.assert_model_framed. Maths is drawn in Unicode text labels (no LaTeX).
"""
from __future__ import annotations

import math
import numpy as np

from . import uns
from .camera import assert_model_framed

STOPS = [(68, 1, 84), (59, 82, 139), (33, 145, 140), (94, 201, 98), (253, 231, 37)]


def colour(x):
    """x in [0,1] → viridis-like hex."""
    x = min(max(x, 0.0), 1.0) * (len(STOPS) - 1)
    i = min(int(x), len(STOPS) - 2)
    f = x - i
    c = [round(STOPS[i][k] + f * (STOPS[i + 1][k] - STOPS[i][k])) for k in range(3)]
    return "#%02x%02x%02x" % tuple(c)


class Panel:
    """A framed plot area mapping data coordinates to SVG pixels."""

    def __init__(self, x, y, w, h, xr, yr, title, equal=False):
        self.x, self.y, self.w, self.h, self.title = x, y, w, h, title
        (x0, x1), (y0, y1) = xr, yr
        if equal:
            s = min(w / (x1 - x0), h / (y1 - y0))
            cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            x0, x1, y0, y1 = cx - w / 2 / s, cx + w / 2 / s, cy - h / 2 / s, cy + h / 2 / s
        self.xr, self.yr = (x0, x1), (y0, y1)
        self.out = []

    def px(self, X, Y):
        return (self.x + (X - self.xr[0]) / (self.xr[1] - self.xr[0]) * self.w,
                self.y + self.h - (Y - self.yr[0]) / (self.yr[1] - self.yr[0]) * self.h)

    def frame(self):
        self.out.append(f'<rect x="{self.x}" y="{self.y}" width="{self.w}" height="{self.h}" fill="none" '
                        f'class="sr"/><text x="{self.x}" y="{self.y - 8}" class="ttl">{self.title}</text>')

    def axes(self, ox=0.0, oy=0.0):
        if self.xr[0] < ox < self.xr[1]:
            a, b = self.px(ox, self.yr[0]), self.px(ox, self.yr[1])
            self.out.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" class="ax"/>')
        if self.yr[0] < oy < self.yr[1]:
            a, b = self.px(self.xr[0], oy), self.px(self.xr[1], oy)
            self.out.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" class="ax"/>')

    def line(self, X, Y, T, tmax, width=1.3):
        """Polyline coloured by t, drawn as short coloured runs."""
        step = max(1, len(X) // 400)
        for i in range(0, len(X) - 1, step):
            j = min(i + step + 1, len(X))
            pts = " ".join("%.1f,%.1f" % self.px(X[k], Y[k]) for k in range(i, j))
            self.out.append(f'<polyline points="{pts}" fill="none" stroke="{colour(T[i] / tmax)}" '
                            f'stroke-width="{width}" stroke-linecap="round"/>')

    def dot(self, X, Y, r=3, fill="fi", label=None):
        a = self.px(X, Y)
        self.out.append(f'<circle cx="{a[0]:.1f}" cy="{a[1]:.1f}" r="{r}" class="{fill}"/>')
        if label:
            self.out.append(f'<text x="{a[0] + 5:.1f}" y="{a[1] - 4:.1f}" class="lab">{label}</text>')

    def text(self, X, Y, s, cls="lab"):
        a = self.px(X, Y)
        self.out.append(f'<text x="{a[0]:.1f}" y="{a[1]:.1f}" class="{cls}">{s}</text>')


def _view_matrix(pose):
    """$vpr-style pose → matrix M; scene points map to view space by p @ M (the camera is rotated by
    (rx,ry,rz), so the scene is rotated by the inverse).  Pose [0,0,0] looks down −z at the xy-plane."""
    rx, ry, rz = (math.radians(a) for a in pose)
    Rx = np.array([[1, 0, 0], [0, math.cos(rx), -math.sin(rx)], [0, math.sin(rx), math.cos(rx)]])
    Ry = np.array([[math.cos(ry), 0, math.sin(ry)], [0, 1, 0], [-math.sin(ry), 0, math.cos(ry)]])
    Rz = np.array([[math.cos(rz), -math.sin(rz), 0], [math.sin(rz), math.cos(rz), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def spring_points(p, floor=-3.0, zscale=0.06):
    """Cartesian xyz of the spring: ρ cos φ, ρ sin φ, zscale·t   (ρ = log10|ζ| − floor, clipped)."""
    rho, phi, z = uns.cylinder(p, floor)
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), zscale * z])


def project(pts, pose):
    """Orthographic: (screen x, screen y, depth). Front [90,0,0] puts z (= t) up."""
    q = pts @ _view_matrix(pose)
    return q[:, 0], q[:, 1], q[:, 2]


def write_csv(p, path, floor=-3.0):
    xyz = spring_points(p, floor)
    rho, phi, _ = uns.cylinder(p, floor)
    cols = ["t", "re_zeta", "im_zeta", "log10_abs_zeta", "theta", "Z", "areal_velocity_L",
            "rho", "phi", "x", "y", "z"]
    data = np.column_stack([p["t"], p["zeta"].real, p["zeta"].imag, np.maximum(p["u"], -9),
                            p["theta"], p["Z"], p["L"], rho, phi, xyz])
    np.savetxt(path, data, delimiter=",", header=",".join(cols), comments="", fmt="%.10g")


# plain class CSS (no custom properties) so non-browser rasterisers render it too; dark mode via media query
CSS = """
.bg{fill:#fff}.fi{fill:#1b1b1f}.fh{fill:#c0392b}.sr{stroke:#b9b9c2}
text{fill:#1b1b1f;font-family:ui-monospace,Menlo,monospace}
.ttl{font-size:13px;font-weight:600}.lab{font-size:10px}.small{font-size:9.5px;fill:#6d6d7a}
.ax{stroke:#8a8a96;stroke-width:.8}
@media (prefers-color-scheme: dark){.bg{fill:#15151a}.fi{fill:#e8e8ee}.fh{fill:#ff7b6b}.sr{stroke:#4a4a55}
text{fill:#e8e8ee}.small{fill:#9a9aa8}.ax{stroke:#7a7a88}}
"""


def flattening_svg(p, zeros, grams, path, floor=-3.0):
    tmax = float(p["t"][-1])
    W, H = 1680, 1180
    pn = []

    # A — the flat pathway (what the web shows)
    A = Panel(40, 50, 520, 520, (-1.9, 3.3), (-2.4, 2.4), "A  flat: Re ζ, Im ζ   (t not drawn)", equal=True)
    A.frame(); A.axes()
    A.line(p["zeta"].real, p["zeta"].imag, p["t"], tmax)
    A.dot(0, 0, 4, "fh")
    A.text(0.35, 0.6, f"hub: {len(zeros)} zeros, one point", "small")
    z0 = p["zeta"][0]
    A.dot(z0.real, z0.imag, 3, label="ζ(½) = −1.4604")
    for g in grams:
        zz = uns.zeta_line(g)
        A.dot(zz.real, zz.imag, 2.4, "fi")
    A.text(-1.8, -2.25, "dots on the real axis = Gram points gₙ, n ≥ 0 (θ = nπ)", "small")
    pn.append(A)

    # B — UNS: log-polar, end-on view of the spring
    rho, phi, _ = uns.cylinder(p, floor)
    Rm = float(rho.max()) * 1.08
    B = Panel(600, 50, 520, 520, (-Rm, Rm), (-Rm, Rm), "B  UNS end-on: ρ = log₁₀|ζ| + 3, angle = arg ζ", equal=True)
    B.frame(); B.axes()
    B.line(rho * np.cos(phi), rho * np.sin(phi), p["t"], tmax)
    B.dot(0, 0, 3, "fh")
    B.text(-Rm * .95, -Rm * .93, "axis = zero set; ρ = 0 at the floor (rendering choice)", "small")
    pn.append(B)

    # C — the unrolled strip: log₁₀|Z| against t, zeros as spikes
    t = p["t"]
    C = Panel(1160, 50, 480, 250, (0, tmax), (floor, 1.0), "C  log₁₀|ζ| against t   (zeros = spikes to the floor)")
    C.frame()
    C.line(t, np.clip(p["u"], floor, 1.0), t, tmax, 1.2)
    for g in zeros:
        C.dot(g, floor, 2.4, "fh")
    pn.append(C)

    # D — the clock: −θ(t)/π with the one-handed winding; Gram = integers
    D = Panel(1160, 340, 480, 230, (0, tmax), (-18, 3), "D  the lost clock: −θ(t)/π  (Gram points at integers)")
    D.frame()
    D.line(t, p["phi_signed"] / np.pi, t, tmax, 1.2)
    for n in range(-16, 3, 2):
        a, b = D.px(0, n), D.px(tmax, n)
        D.out.append(f'<line x1="{a[0]:.0f}" y1="{a[1]:.0f}" x2="{b[0]:.0f}" y2="{b[1]:.0f}" class="sr" stroke-width=".4"/>')
    D.text(2, 2.0, "turning point θ′ = 0 at t = 2π; rate −θ′ = −½ ln(t/2π)", "small")
    pn.append(D)

    # E, F, G — the spring in three canonical poses (scad-spatial: front, iso, top)
    pts = spring_points(p, floor)
    assert_model_framed([tuple(map(float, q)) for q in pts[::40]])
    poses = [("E  spring, front [90,0,0]: height = t", (90, 0, 0)),
             ("F  spring, iso [55,0,25]", (55, 0, 25)),
             ("G  spring, top [0,0,0] = panel B", (0, 0, 0))]
    for k, (name, pose) in enumerate(poses):
        X, Y, _ = project(pts, pose)
        s = float(max(np.ptp(X), np.ptp(Y) * 500 / 480)) / 2 * 1.1
        cx, cy = float(X.min() + X.max()) / 2, float(Y.min() + Y.max()) / 2
        P = Panel(40 + k * 560, 640, 500, 480, (cx - s, cx + s), (cy - s * 480 / 500, cy + s * 480 / 500), name, equal=True)
        P.frame()
        P.line(X, Y, t, tmax, 1.0)
        # gnomon: z axis (t) arrow from origin
        o = project(np.array([[0.0, 0.0, 0.0], [0, 0, 0.06 * tmax]]), pose)
        P.out.append(f'<polyline points="{P.px(o[0][0], o[1][0])[0]:.1f},{P.px(o[0][0], o[1][0])[1]:.1f} '
                     f'{P.px(o[0][1], o[1][1])[0]:.1f},{P.px(o[0][1], o[1][1])[1]:.1f}" class="ax" stroke-dasharray="4 3"/>')
        pn.append(P)

    legend = ('<text x="40" y="1160" class="small">colour = t (violet 0 → yellow %g).  Orthographic.  '
              'Spring: ρ = log₁₀|ζ| + 3 (clipped at 0), φ = arg ζ, z = 0.06·t.  Zeros sit on the axis.</text>' % tmax)
    body = "".join("".join(q.out) for q in pn)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">'
           f'<style>{CSS}</style><rect width="{W}" height="{H}" class="bg"/>{body}{legend}</svg>')
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)
    return path
