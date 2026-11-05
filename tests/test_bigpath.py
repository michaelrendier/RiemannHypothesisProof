# SPDX-License-Identifier: GPL-3.0-only
"""Tests for the long pathway, the flow, the argument-principle accounting and the machinery (notebook 02)."""
import numpy as np
import pytest
from mpmath import mpf, siegelz

from rhp import bigpath, flow, epicycles as ep

ZEROS = np.loadtxt(__import__("pathlib").Path(__file__).resolve().parents[1] / "data" / "zeros_1_to_3600.csv")


def test_zero_cache():
    assert len(ZEROS) == 3066 and abs(ZEROS[0] - 14.134725141734693) < 1e-9 and ZEROS[-1] < 3600
    assert np.all(np.diff(ZEROS) > 0)


def test_rs_matches_mpmath():
    rng = np.random.default_rng(1)
    for x in rng.uniform(100, 3600, 6):
        assert abs(bigpath.Z_rs(np.array([x]))[0] - float(siegelz(mpf(x)))) < 2e-3


def test_grid_zeros_match_cache():
    p = bigpath.pathway_big(3600.0, 0.02)
    z = bigpath.zeros_from_Z(p["t"], p["Z"])
    assert len(z) == 3065 and np.max(np.abs(z - ZEROS[ZEROS > 20])) < 1e-2


def test_flow_symmetry():
    du, dv = flow.symmetry_residuals([0.3 + 5j, 0.1 + 14.5j, 0.45 + 30j, 0.35 + 47j])
    assert du < 1e-12 and dv < 1e-12


def test_flow_sign_grid():
    sig = np.array([0.1, 0.3, 0.45, 0.55, 0.7, 0.9])
    U, V = flow.flow_grid(sig, np.linspace(0.5, 40, 25))
    assert np.all(U * (sig - 0.5)[None, :] > 0)


def test_argument_principle_accounting():
    mids = (ZEROS[:-1] + ZEROS[1:]) / 2
    rows = flow.argument_accounting([mids[3], mids[100], mids[1000], mids[2500]], ZEROS)
    assert all(abs(r["gap"]) < 1e-6 and r["step_max"] < 0.5 for r in rows)


def test_afe_modulus_balance():
    t = np.linspace(20, 80, 6)
    A, cB, chi, zt = ep.afe_parts(0.5, t)
    assert np.max(np.abs(np.abs(A) - np.abs(cB))) < 1e-10 and np.allclose(np.abs(chi), 1, atol=1e-10)
    A4, cB4, _, _ = ep.afe_parts(0.4, t)
    assert np.all(np.abs(cB4) > np.abs(A4) * 0.8) and np.median(np.abs(cB4) / np.abs(A4)) > 1.0


def test_wheels_converge_for_sigma_gt_1():
    t = np.linspace(0, 20, 201)
    tab = dict(ep.wheel_table(1.5, t, [1, 10, 100]))
    assert tab[1] > tab[10] > tab[100] and tab[100] < 2e-2


def test_explicit_formula_improves_with_zeros():
    x = np.linspace(2.0, 60.0, 500)
    psi = ep.psi_exact(x)
    e = [np.sqrt(np.mean((ep.psi_zeros(x, ZEROS[:K]) - psi) ** 2)) for K in (10, 100, 3066)]
    assert e[0] > e[1] > e[2] and e[2] < 0.15
