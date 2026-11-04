# Vendored verbatim from ~/.claude/skills/scad-spatial/assets/camera.py (author's own skill assets), 2026-11-03.
# Provenance: OURS (methodology skill); unchanged apart from this header.
# SPDX-License-Identifier: GPL-3.0-only
"""Pure-Python camera / frustum — the 3D analogue of nes-viewport's Viewport.

No numpy. Column-vector convention; matrices are row-major nested lists, applied
as M @ v. The one gate that matters is `assert_model_framed(points, aspect)` —
for every canonical pose, every model vertex must land inside NDC [-1,1]^3 and
the model must fill a healthy fraction of the frame (not a speck, not clipped).

    pts = cube_points(10)
    assert_model_framed(pts)          # raises on the first bad pose
    cam = Camera.look_at(eye=(30,-40,25), target=(0,0,0))
    ndc = [cam.project(p, aspect=4/3) for p in pts]
"""
from __future__ import annotations

import math
from dataclasses import dataclass

Vec = tuple[float, float, float]


# ── vec / mat ──────────────────────────────────────────────────────────────

def sub(a: Vec, b: Vec) -> Vec: return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def add(a: Vec, b: Vec) -> Vec: return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def scl(a: Vec, s: float) -> Vec: return (a[0]*s, a[1]*s, a[2]*s)
def dot(a: Vec, b: Vec) -> float: return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def cross(a: Vec, b: Vec) -> Vec:
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a: Vec) -> Vec:
    n = math.sqrt(dot(a, a)) or 1.0
    return (a[0]/n, a[1]/n, a[2]/n)


def matmul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def apply(M, v):  # v is (x,y,z); returns (x,y,z,w)
    x, y, z = v
    r = [M[i][0]*x + M[i][1]*y + M[i][2]*z + M[i][3] for i in range(4)]
    return tuple(r)


# ── camera ─────────────────────────────────────────────────────────────────

@dataclass
class Camera:
    view: list                      # world -> view (4x4)
    fovy_deg: float = 35.0
    near: float = 0.1
    far: float = 10_000.0

    @staticmethod
    def look_at(eye: Vec, target: Vec, up: Vec = (0.0, 0.0, 1.0),
                fovy_deg: float = 35.0) -> "Camera":
        f = norm(sub(target, eye))          # forward
        s = norm(cross(f, up))              # right
        u = cross(s, f)                     # true up
        # view: rows are [s; u; -f], translated by -R·eye
        V = [
            [s[0], s[1], s[2], -dot(s, eye)],
            [u[0], u[1], u[2], -dot(u, eye)],
            [-f[0], -f[1], -f[2], dot(f, eye)],
            [0, 0, 0, 1],
        ]
        return Camera(view=V, fovy_deg=fovy_deg)

    def proj(self, aspect: float):
        t = math.tan(math.radians(self.fovy_deg) / 2.0)
        n, fr = self.near, self.far
        return [
            [1.0/(aspect*t), 0, 0, 0],
            [0, 1.0/t, 0, 0],
            [0, 0, -(fr+n)/(fr-n), -2*fr*n/(fr-n)],
            [0, 0, -1, 0],
        ]

    def project(self, p: Vec, aspect: float = 4/3):
        """world point -> NDC (x,y,z) in [-1,1]^3 when inside the frustum."""
        clip = apply(matmul(self.proj(aspect), self.view), p)
        w = clip[3] if abs(clip[3]) > 1e-12 else 1e-12
        return (clip[0]/w, clip[1]/w, clip[2]/w)

    def in_frustum(self, p: Vec, aspect: float = 4/3, pad: float = 1e-6) -> bool:
        x, y, z = self.project(p, aspect)
        return (-1-pad <= x <= 1+pad and -1-pad <= y <= 1+pad and -1-pad <= z <= 1+pad)


# ── model helpers ──────────────────────────────────────────────────────────

def bbox(points):
    xs, ys, zs = zip(*points)
    return ((min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs)))


def centre_radius(points):
    lo, hi = bbox(points)
    c = tuple((lo[i]+hi[i])/2 for i in range(3))
    r = max(math.dist(p, c) for p in points) or 1.0
    return c, r


def cube_points(side: float = 10.0):
    h = side/2
    return [(sx*h, sy*h, sz*h) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]


# canonical poses: (name, euler_deg about world axes)  — matches OpenSCAD $vpr
CANON = {
    "front":  (90, 0, 0),
    "back":   (90, 0, 180),
    "top":    (0, 0, 0),
    "bottom": (180, 0, 0),
    "right":  (90, 0, 90),
    "left":   (90, 0, -90),
    "iso":    (55, 0, 25),
}


def _fit_distance(radius: float, fovy_deg: float, aspect: float, margin: float = 1.25) -> float:
    """Distance at which a bounding sphere of `radius` subtends < the fov on both
    axes, with headroom."""
    half_v = math.radians(fovy_deg) / 2.0
    half_h = math.atan(math.tan(half_v) * aspect)
    return radius / math.sin(min(half_v, half_h)) * margin


def _eye_for(pose_deg, centre, radius, dist_factor=3.0):
    rx, ry, rz = (math.radians(a) for a in pose_deg)
    # OpenSCAD $vpr: rotate the *scene* by (rx,ry,rz); camera sits on +view axis.
    # Equivalent camera direction = R_z R_y R_x applied to (0,-1,0) then scaled.
    def rotx(v): x, y, z = v; return (x, y*math.cos(rx)-z*math.sin(rx), y*math.sin(rx)+z*math.cos(rx))
    def roty(v): x, y, z = v; return (x*math.cos(ry)+z*math.sin(ry), y, -x*math.sin(ry)+z*math.cos(ry))
    def rotz(v): x, y, z = v; return (x*math.cos(rz)-y*math.sin(rz), x*math.sin(rz)+y*math.cos(rz), z)
    d = rotz(roty(rotx((0.0, -1.0, 0.0))))
    return add(centre, scl(norm(d), dist_factor * radius))


def assert_model_framed(points, aspect: float = 4/3, fovy_deg: float = 35.0,
                        min_fill: float = 0.10, max_fill: float = 0.98) -> None:
    """For every canonical pose: all points inside NDC, and the model's NDC
    bounding box fills between min_fill and max_fill of the frame on each axis.
    Raises AssertionError naming the pose and the failure."""
    c, r = centre_radius(points)
    dist = _fit_distance(r, fovy_deg, aspect)
    for name, pose in CANON.items():
        rx, ry, rz = (math.radians(a) for a in pose)

        def rotx(v): x, y, z = v; return (x, y*math.cos(rx)-z*math.sin(rx), y*math.sin(rx)+z*math.cos(rx))
        def roty(v): x, y, z = v; return (x*math.cos(ry)+z*math.sin(ry), y, -x*math.sin(ry)+z*math.cos(ry))
        def rotz(v): x, y, z = v; return (x*math.cos(rz)-y*math.sin(rz), x*math.sin(rz)+y*math.cos(rz), z)
        d = rotz(roty(rotx((0.0, -1.0, 0.0))))
        eye = add(c, scl(norm(d), dist))
        up = (0.0, 0.0, 1.0)
        if abs(dot(norm(sub(c, eye)), up)) > 0.999:      # looking straight down/up
            up = (0.0, 1.0, 0.0)
        cam = Camera.look_at(eye, c, up, fovy_deg=fovy_deg)
        ndc = [cam.project(p, aspect) for p in points]
        for (x, y, z), p in zip(ndc, points):
            assert -1.001 <= x <= 1.001 and -1.001 <= y <= 1.001 and -1.001 <= z <= 1.001, (
                f"pose '{name}': point {tuple(round(v, 2) for v in p)} projects "
                f"outside NDC at ({x:.2f},{y:.2f},{z:.2f}) — clipped or camera too close")
        xs = [p[0] for p in ndc]; ys = [p[1] for p in ndc]
        fx = (max(xs) - min(xs)) / 2.0
        fy = (max(ys) - min(ys)) / 2.0
        assert min_fill <= max(fx, fy) <= max_fill, (
            f"pose '{name}': model fills {max(fx, fy):.0%} of the frame "
            f"(want {min_fill:.0%}–{max_fill:.0%}) — dolly in/out or change fov")


if __name__ == "__main__":
    for side in (1, 10, 250):
        assert_model_framed(cube_points(side))
        print(f"ok  cube side {side}  — framed from every canonical pose")
    # a deliberately off-centre model must still frame (centre_radius handles it)
    pts = [(x+500, y-300, z+40) for x, y, z in cube_points(30)]
    assert_model_framed(pts)
    print("ok  off-centre model framed")
    # negative: a camera hard-coded too close should be caught by in_frustum
    cam = Camera.look_at(eye=(0, -2, 0), target=(0, 0, 0), fovy_deg=35)
    far_pt = (0, 0, 0)
    assert cam.in_frustum(far_pt)
    assert not cam.in_frustum((20, 0, 0)), "point far off-axis should be outside"
    print("ok  frustum test discriminates")
    print("all camera self-tests passed")
