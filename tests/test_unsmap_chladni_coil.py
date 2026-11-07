# SPDX-License-Identifier: GPL-3.0-only
"""Tests: rhp.unsmap (number-free), rhp.chladni (parities, flux, lobe), rhp.coil (spiral laws)."""
import math
import numpy as np
import pytest
from scipy.special import fresnel

from rhp import unsmap as U, chladni as C, coil as K


# ── unsmap: shapes with known answers ──
def test_critical_line_is_a_great_circle_other_lines_are_small_circles():
    t = np.linspace(-80, 80, 4001)
    n, d, res = U.plane_fit(U.lift(0.5 + 1j * t)); assert abs(d) < 1e-9 and res < 1e-12
    n, d, res = U.plane_fit(U.lift(0.8 + 1j * t)); assert abs(d) > 0.3 and res < 1e-12       # still a circle (planar), not great
    th = np.linspace(0, 2 * np.pi, 2000)
    assert U.plane_fit(U.lift(2 + 0.7 * np.exp(1j * th)))[2] < 1e-12                          # a circle in the plane is a circle on the sphere


def test_genuine_cusp_survives_the_lift():
    c = np.linspace(-1, 1, 4001)
    z = 0.6 + 0.5j + 0.4 * (c ** 2 + 1j * c ** 3)
    f, g = U.cusp_survives_lift(z, c, 0.0)
    assert f < 1e-5 and g < 1e-5


def test_edge_on_rotation_vs_not():
    t = np.linspace(0, 60, 6000)
    r = U.edge_on_rotation_test(t, np.cos(2 * np.pi * 0.5 * t))
    assert r["amp_cv"] < 1e-2 and r["phase_linear_r2"] > 0.999999
    bad = U.edge_on_rotation_test(t, np.sin(t) * np.exp(-t / 10))
    assert bad["amp_cv"] > 0.5


def test_runaway_is_smooth_over_the_pole():
    t = np.linspace(-200, 200, 20001)
    w = 1 + 1j * t                                   # the line Re w = 1: flies to infinity both ways
    wmax, wmin, kg = U.pole_passage(w, t)
    assert wmax > 100 and kg < 1.5                   # sphere geodesic curvature finite (a small circle)


def test_spiral_laws_identify_the_spiral():
    ph = np.linspace(0, 12 * np.pi, 3000)
    a = U.spiral_laws((1 + .3 * ph) * np.exp(1j * ph)); assert a["archimedean"] > 0.9999 and a["logarithmic"] < 0.95
    l = U.spiral_laws(np.exp(.12 * ph) * np.exp(1j * ph)); assert l["logarithmic"] > 0.9999 and l["archimedean"] < 0.9
    x = np.linspace(0.05, 6, 3000); S, Cc = fresnel(x)
    e = U.spiral_laws(Cc + 1j * S, Cc[-1] + 1j * S[-1]); assert e["euler"] > 0.9999


def test_sphere_curvature_formula_on_a_small_circle():
    th = np.linspace(0, 2 * np.pi, 4001)
    P = np.column_stack([0.6 * np.cos(th), 0.6 * np.sin(th), 0.8 * np.ones_like(th)])      # small circle of angular radius a = arcsin 0.6: geodesic curvature cot a = 0.8/0.6
    kg = U.sphere_geodesic_curvature(P, th)
    assert np.allclose(np.abs(kg[10:-10]), 0.8 / 0.6, atol=1e-3)


# ── chladni ──
def test_mirror_is_conjugation_and_flux_vanishes_on_the_equator():
    wre, wim = C.parity_defect([0.05, 0.3], [14.0, 30.0]); assert wre < 1e-12 and wim < 1e-12
    a, b = C.flux_two_ways(0.7, 30.0); assert abs(a - b) < 1e-9 * (abs(a) + 1e-30)           # J two ways (Cauchy–Riemann)
    assert abs(C.flux_two_ways(0.5, 30.0)[0]) < 1e-20


def test_north_flux_positive_and_offline_pocket_detected():
    sig = np.linspace(0.3, 0.8, 11); t = np.linspace(13, 27, 36)
    J = C.rel_flux_grid(sig, t); assert C.north_violations(J, sig) == 0
    f = C.with_offline_quadruple(0.08, 25.0)
    s2 = np.linspace(0.5, 0.72, 23)[1:]; t2 = np.linspace(24, 26, 41)
    assert int((C.rel_flux_grid(s2, t2, f) <= 0).sum()) >= 1


def test_lobe_angle_mean_is_pi():
    z = np.loadtxt("data/zeros_1_to_3600.csv")
    la = C.lobe_angles(z[:600])
    assert abs(la.mean() - math.pi) < 0.01 and 0.9 < la.std() < 1.5                         # GUE spacing × π


# ── coil ──
def test_pancake_geometry_and_taps():
    x, y, phi, r = K.pancake()
    assert abs((r[-1] - r[0]) - 0.8 * 17) < 1e-9
    z = np.loadtxt("data/zeros_1_to_3600.csv")
    assert [v[0] for v in K.tap_zeros(z).values()] == [2, 6, 10, 34]


def test_index_radius_makes_hardy_angle_archimedean_height_does_not():
    z = np.loadtxt("data/zeros_1_to_3600.csv")
    rl = K.radius_law(z[:300])
    assert rl["index_vs_angle"] > 0.99998 and rl["height_vs_angle"] < rl["index_vs_angle"]
    p = K.hardy_pitch_t([20.0, 200.0, 2000.0]); assert p[0] > p[1] > p[2]                    # height-wound pitch tightens


@pytest.mark.parametrize("sigma", [0.5, 0.8, 1.0, 1.5])
def test_partial_sum_spiral_is_logarithmic_with_predicted_slope(sigma):
    meas, pred = K.tail_log_slope(100, N0=800, N1=8000, sigma=sigma)
    assert abs(meas - pred) < 2e-5


def test_axis_hole_widens_off_the_line():
    d = K.axis_distance(np.linspace(10, 40, 2500), (0.5, 0.6, 0.75, 1.0))
    assert d[0.5] < d[0.6] < d[0.75] < d[1.0]


def test_explicit_formula_trace_is_two_re_q():
    from rhp import diverge as D
    z = np.loadtxt("data/zeros_1_to_3600.csv")[:300]
    psi = D.psi_table(2000)
    u = np.linspace(3.0, 6.0, 2000)
    E = D.residual(u, psi)
    Q = K.explicit_path(u, z)
    assert np.corrcoef(E, 2 * Q.real)[0, 1] > 0.9
