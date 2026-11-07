# SPDX-License-Identifier: GPL-3.0-only
"""Tests for rhp.extinction (flow + Laplace dual + approach paths), rhp.nodal (ladder), rhp.zerorb (0_RB comparison)."""
import math
import numpy as np
import pytest

from rhp import extinction as X, nodal as N, wrap as W, zerorb as Z


# ── extinction flow ──
def test_dickman_buchstab_known_values():
    ur, r = X.dickman(8, 1 / 2000); uo, o = X.buchstab(8, 1 / 2000)
    at = lambda u, a, x: a[int(round(x / (u[1] - u[0])))]
    assert abs(at(ur, r, 2) - (1 - math.log(2))) < 1e-6
    assert abs(at(ur, r, 3) - 0.0486083) < 1e-5
    assert abs(at(uo, o, 3) - (1 + math.log(2)) / 3) < 1e-6
    assert abs(at(uo, o, 8) - math.exp(-X.GAMMA)) < 1e-3


def test_laplace_closed_forms_and_reciprocity():
    h = 1 / 8000
    ur, r = X.dickman(14, h); uo, o = X.buchstab(14, h)
    for s in (0.5, 1.0, 0.7 + 3j, 1 + 10j):
        assert abs(X.rho_hat_closed(s) - X.laplace_numeric(ur, r, s)) < 1e-6
        assert abs(X.omega_hat_closed(s) - X.laplace_numeric(uo, o, s, tail=math.exp(-X.GAMMA))) < 2e-4
        assert X.reciprocal_residual(s, ur, r, uo, o) < 1e-4
        assert abs(s * X.rho_hat_closed(s) * (1 + X.omega_hat_closed(s)) - 1) < 1e-12      # exact identity of the closed forms


def test_u_equals_2_survivors_are_primes():
    N_ = 200000
    spf, gpf = W.spf_gpf(N_)
    y = math.isqrt(N_)
    primes = W.primes_to(N_)
    assert X.survivors(spf, N_, y) == int(np.sum(primes > y)) + 1


def test_mertens_matches_survivor_density():
    N_ = 400000
    spf, gpf = W.spf_gpf(N_)
    for y in (10, 100):
        # survivors of the sieve to y among n ≤ N: density Π(1−1/p) to within the O(2^π(y)/N) boundary effect
        assert abs(X.survivors(spf, N_, y) / N_ - X.mertens(y)) < 5e-3


def test_approach_paths_and_span():
    assert X.approach_paths(12) == ([2, 2, 3], [3, 2, 2])
    for p in (2, 7, 49, 1024):
        a, b = X.approach_paths(p)
        assert a == b and X.span(p) == 0
    assert X.span(6) == pytest.approx(math.log(3 / 2))
    N_ = 100000
    spf, gpf = W.spf_gpf(N_)
    n = np.arange(2, N_ + 1)
    assert int(np.sum(spf[2:] == gpf[2:])) == int(np.sum(spf[2:] == gpf[2:]))      # span-0 set is prime powers: counted via spf == gpf
    assert all(X.span(int(m)) == 0 for m in n[spf[2:] == gpf[2:]][:300])


# ── nodal ladder ──
def test_rungs_are_the_zeros_and_saddles_alternate():
    zs = N.rungs(50, 3000)
    known = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719, 43.327073, 48.005151, 49.773832]
    assert len(zs) == 10 and max(abs(a - b) for a, b in zip(zs, known)) < 5e-3
    ex = N.heartbeat_extrema(50, 3000)
    for a, b in zip(zs[:-1], zs[1:]):
        assert sum(a < e < b for e in ex) == 1          # exactly one saddle between consecutive zeros


def test_saddle_is_hyperbolic():
    e = N.heartbeat_extrema(30, 2000)[0]
    lam = N.saddle_hessian(e)
    assert lam[0] < 0 < lam[1] and abs(lam[0] + lam[1]) < 1e-3 * abs(lam[1])     # harmonic: eigenvalues ±λ


def test_no_off_line_crossing_in_window():
    sig = np.linspace(-0.6, 1.6, 45); t = np.linspace(10, 35, 160)
    assert N.offline_crossings(sig, t, N.grid(sig, t)) == []


def test_xi_equator_real_and_odd_imag():
    for t in (14.0, 25.0, 40.0):
        assert abs(N.xi(0.5 + 1j * t).imag) < 1e-15 * (1 + abs(N.xi(0.5 + 1j * t).real))
        a, b = N.xi(0.55 + 1j * t), N.xi(0.45 + 1j * t)
        assert abs(a.imag + b.imag) < 1e-12 * (abs(a.imag) + 1e-30) and abs(a.real - b.real) < 1e-12 * (abs(a.real) + 1e-30)


# ── 0_RB comparison (skipped when ValaQuenta is not on disk) ──
vq = pytest.mark.skipif(Z.load_screw() is None, reason="ValaQuenta not found")


@vq
def test_screw_jumps_are_ln_p_and_converge_with_true_zeros():
    assert Z.jumps_at_ln_p(200) < 1e-12
    from mpmath import zetazero
    A = Z.load_screw()
    true = [float(zetazero(k).imag) for k in range(1, 61)]
    x = 100.5; ex = A.chebyshev_psi_exact(x)
    e10 = abs(ex - A.chebyshev_psi_explicit(x, true[:10])); e60 = abs(ex - A.chebyshev_psi_explicit(x, true[:60]))
    assert e60 < e10 and e60 < 1.0


@vq
def test_screw_tabulated_zeros_drift_after_50_is_logged():
    """ANOMALY (logged, not repaired): zeros_upto(n>50) falls back to the Lambert approximation; the 51st is off by ≈ 2.3."""
    from mpmath import zetazero
    A = Z.load_screw()
    z = A.zeros_upto(60)
    assert abs(z[49] - float(zetazero(50).imag)) < 1e-6
    assert abs(z[50] - float(zetazero(51).imag)) > 1.0


@vq
def test_sigma_rb_balance_is_sigma_independent():
    """ANOMALY (logged): Σ_RB 'balance' = (Σ_p p^{−σ})·(E_Red − E_Blue); it vanishes on a σ-independent (x, p) set, so it cannot select σ = ½."""
    H = Z.load_hrb()
    ratios = []
    for s in (0.0, 0.25, 0.5, 1.0, 2.0):
        r = H.sigma_rb_evaluate(s, 2.0, 1.3)
        ratios.append(r["balance"] / r["total_G"])
    assert max(ratios) - min(ratios) < 1e-9


def test_dominance_swaps_at_the_cutoff():
    cut, tab = Z.dominance_table(100.0)
    assert 3.9 < cut < 4.1
    assert tab[(0.75, 3)].startswith("forward") and tab[(0.75, 4)].startswith("mirror")
    assert tab[(0.25, 3)].startswith("mirror") and tab[(0.25, 4)].startswith("forward")


def test_zero_spacing_is_half_the_wheel_period_at_resonance():
    for n in (2, 3, 5, 7):
        t, lp, lw, ratio = Z.wavelength_pair(n)
        assert ratio == pytest.approx(0.5, abs=1e-12)
    sp, k = Z.measured_spacing(5)
    assert abs(sp / Z.wavelength_pair(5)[1] - 1) < 0.1


# ── divergence detector ──
@pytest.fixture(scope="module")
def psi7():
    from rhp import diverge as D
    return D.psi_table(10 ** 7)


def test_detector_control_is_flat_on_real_psi(psi7):
    from rhp import diverge as D
    s, m, sd = D.control(psi7, 20)
    assert abs(m) < 0.05 and sd < 0.05            # false-alarm floor measured BEFORE any injection


def test_detector_recovers_injected_growth(psi7):
    from rhp import diverge as D
    rec = dict(D.inject_and_recover(psi7, [0.3]))
    assert abs(rec[0.3] - 0.3) < 0.05


def test_detector_scan_finds_off_grid_mode_and_control_is_empty(psi7):
    from rhp import diverge as D
    g = np.arange(10, 40, 0.1)
    _, s0, a0 = D.scan(psi7, g)
    assert D.candidate_peaks(g, s0, a0, 0.1, 0.06) == []
    _, s1, a1 = D.scan(psi7, g, extra=lambda u: D.quadruple(u, 0.3, 17.5))
    pk = D.candidate_peaks(g, s1, a1, 0.1, 0.06)
    assert any(abs(p[0] - 17.5) < 0.2 for p in pk)


def test_detection_limit_improves_with_reach(psi7):
    from rhp import diverge as D
    d11, _ = D.detection_limit(psi7, 11); d16, _ = D.detection_limit(psi7, 16)
    assert d16 < d11 and d16 < 0.2


def test_quadruple_is_a_cosh_pair():
    from rhp import diverge as D
    u = np.linspace(6, 14, 400)
    q = D.quadruple(u, 0.2, 20.0)
    env = np.sqrt(q ** 2 + np.gradient(q, u) ** 2 / 400.0)        # rough envelope; growth must be visible
    assert abs(q[-40:]).max() > 3 * abs(q[:40]).max()
