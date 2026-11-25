# SPDX-License-Identifier: GPL-3.0-only
"""Tests for rhp.strange_path: the RK4 integrator's order, the port's layout, the attractor stays bounded, and the zeta path sits on the sphere with zeros at the south pole."""
import numpy as np
from rhp import strange_path as P

def test_rk4_is_fourth_order():
    f = lambda s: np.array([-s[0], 0.0, 0.0]); x0 = np.array([1.0, 0, 0]); T = 1.0
    def err(h):
        s = x0.copy()
        for _ in range(int(round(T / h))): s = P.rk4_step(f, s, h)
        return abs(s[0] - np.exp(-T))
    ratio = err(0.1) / err(0.05); assert 14 < ratio < 18                         # halving h divides the error by 2^4 = 16

def test_layout_is_x_z_y_and_line_distance_runs_zero_to_one():
    rhs = P.aizawa(); pos, line = P.generate_geometry(rhs, points=500, burn_in=300, scale=2.0)
    s = np.array([0.1, 0, 0.0]); 
    for _ in range(300): s = P.rk4_step(rhs, s, 0.01)
    s = P.rk4_step(rhs, s, 0.01)
    assert np.allclose(pos[:3], [2 * s[0], 2 * s[2], 2 * s[1]], atol=1e-5)       # slot 2 is z, slot 3 is y
    assert line[0] == 0 and abs(line[-1] - 499 / 500) < 1e-6 and np.all(np.diff(line) > 0)

def test_aizawa_orbit_is_bounded_and_fills_a_shell():
    pos, _ = P.generate_geometry(P.aizawa(), points=20000, burn_in=2000); X = pos.reshape(-1, 3)
    assert np.all(np.isfinite(X)) and np.abs(X).max() < 4.0
    r = np.linalg.norm(X, axis=1); assert r.max() / max(r.min(), 1e-9) > 2          # a hollow-looking shell with an axial tube, not a point

def test_sphere_point_poles_and_norm():
    assert np.allclose(P.sphere_point(0j), [0, 0, -1]) and abs(np.linalg.norm(P.sphere_point(3 - 4j)) - 1) < 1e-14
    assert P.sphere_point(1e9 + 0j)[2] > 1 - 1e-8                                   # w -> infinity is the north pole

def test_zeta_path_is_on_the_sphere_and_the_first_zero_is_the_south_pole():
    pos, line, ts, w = P.zeta_sphere_path(points=300, t0=14.0, t1=14.3)
    X = pos.reshape(-1, 3); assert np.max(np.abs(np.linalg.norm(X, axis=1) - 1)) < 1e-6
    i = int(np.argmin(np.abs(w))); assert abs(ts[i] - 14.134725) < 0.002             # nearest approach to w = 0 is the first Riemann zero
    assert np.linalg.norm(P.sphere_point(w[i]) - np.array([0, 0, -1])) < 0.01
