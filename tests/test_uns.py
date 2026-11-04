# SPDX-License-Identifier: GPL-3.0-only
"""Regression tests for the flattening map (rhp.uns, rhp.flatten, rhp.render)."""
import math
import numpy as np
import pytest
from mpmath import mp, mpf, mpc, zeta, siegelz, siegeltheta

from rhp import uns, flatten, render

T, N = 50.0, 4001


@pytest.fixture(scope="module")
def p():
    return uns.pathway(T, N)


@pytest.fixture(scope="module")
def zeros():
    return uns.zeros_below(T)


def test_ten_zeros_below_50(zeros):
    assert len(zeros) == 10 and abs(zeros[0] - 14.134725) < 1e-5


def test_hardy_form(p):
    # ζ = e^{-iθ} Z with Z real: Im(e^{iθ}ζ) = 0, and Z agrees with mpmath's siegelz
    w = np.exp(1j * p["theta"]) * p["zeta"]
    assert np.max(np.abs(w.imag)) < 1e-10
    for i in (100, 1500, 3000):
        assert abs(p["Z"][i] - float(siegelz(mpf(p["t"][i])))) < 1e-9


def test_areal_velocity_identity():
    # Im(ζ̄ ζ′) = −θ′ Z² with the exact derivative of ζ
    for t in (3.0, 10.0, 14.0, 20.0, 27.7, 41.3):
        s = mpc(0.5, t)
        z = complex(zeta(s))
        dz = 1j * complex(zeta(s, derivative=1))
        L = (z.conjugate() * dz).imag
        Z = float(siegelz(mpf(t)))
        assert abs(L - (-uns.theta_prime(t) * Z * Z)) < 1e-9


def test_hub_direction_is_minus_theta(p, zeros):
    rows = flatten.hub_passes(p, zeros)
    assert len(rows) == 10 and all(r["agree"] for r in rows)
    # the petal sign flips at every simple zero
    assert all(r["sign_before"] == -r["sign_after"] for r in rows)


def test_gram_points_cross_real_axis():
    for g in uns.gram_points(40):
        z = uns.zeta_line(g)
        assert abs(z.imag) < 1e-9 and abs(uns.theta(g) / math.pi - round(uns.theta(g) / math.pi)) < 1e-9


def test_crossings_regression(p):
    c = flatten.false_crossings(p)
    hub = [r for r in c if r["hub"]]
    # ten passes through one point = C(10,2) pairwise crossings; the rest are false crossings
    assert len(hub) == math.comb(10, 2) == 45
    assert len(c) - len(hub) == 21
    assert all(r["dt"] > 0 for r in c)


def test_lost_clock(p, zeros):
    r = flatten.lost_clock(p, zeros)
    assert r["max_abs_dev"] < 1e-3
    assert abs(r["sign_change_t"] - 2 * math.pi) < 0.05


def test_mirror():
    assert flatten.mirror_check() < 1e-12


def test_gram_law_first_failure():
    g = flatten.gram_law(300.0)
    assert g["failures"][0]["n"] == 125 and g["failures"][0]["zeros"] == 0
    assert all(f["n"] >= 125 for f in g["failures"])


def test_viewport_sweep_covers_everything(p, zeros):
    c = flatten.false_crossings(p)
    rows = flatten.window_map(p, zeros, c, vh=500)
    assert rows[0]["t_from"] == 0.0 and abs(rows[-1]["t_to"] - T) < 1e-9
    assert sum(r["zeros"] for r in rows[:-1]) >= 8


def test_svg_and_csv(p, zeros, tmp_path):
    out = render.flattening_svg(p, zeros, uns.gram_points(T), tmp_path / "f.svg")
    s = open(out, encoding="utf-8").read()
    assert s.startswith("<svg") and "var(--" not in s
    render.write_csv(p, tmp_path / "p.csv")
    assert (tmp_path / "p.csv").read_text().splitlines()[0].startswith("t,re_zeta")
