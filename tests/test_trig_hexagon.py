# SPDX-License-Identifier: GPL-3.0-only
"""Tests for rhp.trig_hexagon: every identity on the hexagon, the A2 structure, the rotation invariances, and trig from the rotation group."""
import numpy as np
import pytest
from rhp import trig_hexagon as T

ANG = np.random.default_rng(1).uniform(0.15, 1.4, 200)

def test_reciprocals_products_pythagoras_hold_at_random_angles():
    for t in ANG:
        v = T.values(t); assert T.reciprocal_ok(v) and T.product_ok(v) and T.pythag_ok(v)

def test_quotients_follow_from_the_product_rule():
    for t in ANG:
        v = T.values(t); assert abs(v['tan'] - v['sin'] / v['cos']) < 1e-9 and abs(v['cot'] - v['cos'] / v['sin']) < 1e-9

def test_log_space_is_A2():
    M, ok_map, ok_order, ok_sum = T.rotation_matrix_action()
    assert ok_map and ok_order and ok_sum                                      # R^6 = I, 60° maps each root to the next, neighbours add to the vertex between
    A = T.a2_coordinates(0.5)
    for t in ANG[:20]:
        v = T.values(t)
        for k in T.RING: assert abs(np.log(abs(v[k])) - (A[k][0] * np.log(np.sin(t)) + A[k][1] * np.log(np.cos(t)))) < 1e-9

def test_which_symmetries_preserve_which_identities():
    inv = T.invariance()
    holds = lambda fam: {s for s, ok in inv[fam] if ok}
    assert len(holds('product')) == 12                                         # the product rule is invariant under the whole dihedral group of the hexagon
    assert len(holds('reciprocal')) == 12                                      # so are the three opposite pairs
    # MEASURED (the 3-fold guess was wrong): the three Pythagorean forms have different hypotenuse roles (1 for sin/cos, sec for tan, csc for cot),
    # so no rotation of the labels preserves them; only the identity and the co-function reflection theta -> pi/2 - theta (sin<->cos, tan<->cot, sec<->csc) do.
    assert holds('pythagorean') == {'rotate 0°', 'reflect+240°'}

def test_all_six_from_the_rotation_matrix():
    for t in ANG:
        v = T.six_from_rotation(t); w = T.values(t)
        assert all(abs(v[k] - w[k]) < 1e-9 for k in T.RING)

def test_sum_double_angle_and_de_moivre_are_the_group_law():
    rng = np.random.default_rng(2)
    for a, b in rng.uniform(0.1, 0.7, (100, 2)):
        s = T.sum_formulas_from_rotation(a, b)
        assert abs(s['sin'] - (np.sin(a) * np.cos(b) + np.cos(a) * np.sin(b))) < 1e-12
        assert abs(s['cos'] - (np.cos(a) * np.cos(b) - np.sin(a) * np.sin(b))) < 1e-12
        assert abs(s['tan'] - (np.tan(a) + np.tan(b)) / (1 - np.tan(a) * np.tan(b))) < 1e-9
    for a in rng.uniform(0.1, 0.7, 50):
        d = T.double_angle_from_rotation(a); assert abs(d['sin'] - 2 * np.sin(a) * np.cos(a)) < 1e-12 and abs(d['cos'] - (np.cos(a) ** 2 - np.sin(a) ** 2)) < 1e-12
        assert np.allclose(T.de_moivre(a, 7), T.R(7 * a), atol=1e-12)
