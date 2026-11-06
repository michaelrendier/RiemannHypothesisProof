# SPDX-License-Identifier: GPL-3.0-only
"""Tests for rhp.sphere (the equatorial construction) and rhp.wrap (primes on the springs)."""
import math
import numpy as np
import pytest
from rhp import sphere as S, wrap as W, uns

ZEROS = np.loadtxt("data/zeros_1_to_3600.csv")


def test_critical_line_is_the_unit_circle():
    assert S.equator_residual(ZEROS) < 1e-12


@pytest.mark.parametrize("s", [0.3 + 10j, 0.7 + 25j, 0.2 + 100j, 0.9 + 7j])
def test_functional_equation_is_the_mirror(s):
    assert S.mirror_residual(s) < 1e-15


def test_mirror_is_inversion_in_the_unit_circle():
    s = 0.8 + 17.0j
    w = S.cayley(s)
    assert abs(S.cayley(1 - s.conjugate()) - 1 / w.conjugate()) < 1e-14      # w ↦ 1/w̄
    assert abs(S.cayley(1 - s) - 1 / w) < 1e-14                              # w ↦ 1/w


def test_hemispheres_and_poles():
    assert abs(S.cayley(2.0 + 3j)) < 1 and abs(S.cayley(-1.0 + 3j)) > 1
    assert abs(S.cayley(1.0 + 0j)) == 0                                       # s = 1: north pole
    r0, r1 = S.pole_residues()
    assert abs(r0 + 1) < 1e-6 and abs(r1 - 1) < 1e-6                         # Λ has residues −1 at 0, +1 at 1


def test_the_line_wraps_once():
    t = np.linspace(-1e4, 1e4, 20001)
    w = S.cayley(0.5 + 1j * t)
    assert np.allclose(np.angle(w), S.longitude_cayley(t), atol=1e-9) or np.allclose(
        np.exp(1j * np.angle(w)), np.exp(1j * S.longitude_cayley(t)), atol=1e-9)


def test_pi_factored_out_of_theta():
    g, p, tot = S.theta_split([10.0, 100.0, 1000.0])
    assert np.max(np.abs(g + p - tot)) < 1e-9
    assert abs(p[1] / p[0] - 10) < 1e-12                                      # the π part is linear in t


def test_zeta_is_pi_power_times_remainder():
    a, b = S.zeta_pi_factored(0.5 + 10j)
    assert abs(complex(a * b) - uns.zeta_line(10.0)) < 1e-12


def test_amplitude_ratio_cross():
    # equator: ratio 1 for every term; off the equator it is 1 only near n = √(t/2π)
    for n in (1, 2, 7, 30):
        assert abs(S.amplitude_ratio(n, 0.5, 250.0) - 1) < 1e-9
    t = 100.0
    N = math.sqrt(t / (2 * math.pi))
    assert abs(S.amplitude_ratio(round(N), 0.3, t) - 1) < 0.01
    assert abs(S.amplitude_ratio(1, 0.3, t) - 1) > 0.5
    assert S.polarisation(1.0) == 0.0


def test_ellipse_is_a_line_when_amplitudes_match():
    A, B = S.ellipse_axes(0.4, 0.4)
    assert B == 0 and abs(A - 0.8) < 1e-15


def test_every_twin_pair_straddles_the_zero_spot_of_6():
    n, ok = W.straddle_check(10 ** 6)
    assert n > 8000 and n == ok


def test_allowed_centres_count_matches_construction():
    for P in (5, 7, 11):
        m, ok = W.allowed_centres(P)
        assert len(ok) == W.allowed_centre_count(P) and 0 in ok


def test_wheel_lag_is_the_gap_in_the_beat():
    assert abs(W.wheel_phase_lag(101, 103, 10.0) - 10 * math.log(103 / 101)) < 1e-15
    assert abs(W.wheel_phase_lag(101, 103, W.antiphase_height(101, 103)) - math.pi) < 1e-12


def test_spring_tightens_and_primes_crowd():
    p = W.hardy_pitch([20.0, 200.0, 2000.0])
    assert p[0] > p[1] > p[2]
    d = W.prime_density_in_t([5.0, 10.0, 15.0])
    assert d[0] < d[1] < d[2]


def test_wheels_help_far_north_and_stop_helping_at_the_equator():
    t = np.linspace(10, 60, 41)
    Ks = [10, 100, 10000]
    north = W.euler_grip(1.5, t, Ks)
    eq = W.euler_grip(0.5, t, Ks)
    assert north[0] > north[1] > north[2]                                     # converges
    assert eq[2] > eq[1]                                                      # more wheels is worse than 100 at σ = ½
