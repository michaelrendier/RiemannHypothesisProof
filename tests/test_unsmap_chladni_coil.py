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


# ── the pasted standing-wave summation, claim by claim (rhp.standing) ──
def test_J_is_order_four_with_no_fixed_point_and_circle_invariant():
    from rhp import standing as S
    c = S.J_cycle(0.6 + 0.3j)
    assert abs(c[4] - c[0]) < 1e-12 and abs(c[2] + c[0]) < 1e-12          # J⁴ = id, J² = −z
    free, cmin, cmax = S.J_fixed_point_defect()
    assert free > 1.0 and abs(cmin - math.sqrt(2)) < 1e-9 and abs(cmax - math.sqrt(2)) < 1e-9   # r = 1 invariant, not fixed


def test_J_is_not_a_symmetry_of_xi_but_the_mirror_is():
    from rhp import standing as S
    ss = [complex(0.5 + dx, t) for dx, t in [(0.1, 20.0), (-0.2, 31.0), (0.3, 44.0)]]
    assert np.all(np.abs(np.log(S.xi_ratio_under(S.J, ss))) > 5)
    assert np.allclose(S.xi_ratio_under(S.mirror, ss), 1, atol=1e-9)


def test_energy_integral_is_not_pi():
    from rhp import standing as S
    for k in (1.0, 14.134725):
        closed, area = S.energy_gemini(k)
        assert abs(closed - S.energy_numeric(k)) < 1e-6 and abs(closed - area) > 0.3


def test_flux_identity_exact_and_mean_equals_theta_prime_for_large_sigma():
    from rhp import standing as S
    for sg, t in [(0.7, 30.0), (1.3, 50.0)]:
        assert S.flux_identity_defect(sg, t) < 1e-8
    a, b, c = S.remainder_terms(1.3, 50.0)
    assert abs(a) < 1e-2 and abs(b) < 1e-3                                  # only Re ζ′/ζ is not O(1/t)


def test_engine_J_matches_and_only_radius_closes():
    from rhp import standing as S
    r = S.engine_J_agrees()
    if r is None:
        pytest.skip("ValaQuenta not found")
    assert r["max_diff"] < 1e-12 and r["radius_closes_after_two"] and abs(r["angle_after_two"] - math.pi) < 1e-12


# ── the Heartbeat: what the standing wave contains ──
def test_zero_circle_area_budget_matches_hadamard_constant():
    from rhp import heartbeat as H
    z = np.loadtxt("data/zeros_1_to_3600.csv")
    s, st, ex = H.budget(z)
    assert abs(st - ex) < 1e-7 and abs(ex - (1 + 0.5772156649015329 / 2 - 0.5 * math.log(4 * math.pi))) < 1e-15
    assert abs(H.SUM_INV_RHO + H.GROUND_SHARE - 1) < 1e-15
    cum = np.cumsum(H.tone_radius(z) ** 2) / H.SUM_INV_RHO
    assert abs(cum[0] - 0.2164) < 5e-4 and abs(cum[9] - 0.586) < 2e-3 and abs(cum[33] - 0.7598) < 2e-3


def test_mean_power_is_parseval_value_on_exact_psi():
    from rhp import heartbeat as H, diverge as D
    psi = D.psi_table(10 ** 7)
    m, _ = H.mean_power_measured(psi, 6.0, 16.1)
    assert abs(m / H.MEAN_POWER - 1) < 0.03
    assert abs(H.MEAN_POWER - (2 + 0.5772156649015329 - math.log(4 * math.pi))) < 1e-15


def test_add_is_spiral_similarity_and_pair_is_real():
    from rhp import heartbeat as H
    g = 14.134725141735
    sc, rot = H.spiral_similarity(g, 1.0)
    assert abs(sc - math.exp(0.5)) < 1e-12 and abs((rot - g) % (2 * math.pi)) < 1e-9 or abs((rot - g) % (2 * math.pi) - 2 * math.pi) < 1e-9
    u = np.linspace(0, 5, 500)
    f, m, s = H.counter_rotating_pair(u, g, 0.0)
    assert np.max(np.abs(s.imag)) < 1e-14 and np.allclose(np.abs(f), np.abs(m))
    f2, _, _ = H.counter_rotating_pair(u, g, 0.15)
    assert abs(f2[-1]) / abs(f2[0]) > 2          # off the line the common radius grows


# ── 65-wheel coil ──
def test_wheel65_exact_primality_below_317_squared():
    from rhp import wheel65 as w
    from rhp.wrap import sieve
    sv = sieve(w.LIMIT - 1)
    assert len(w.PRIMES) == 65 and w.PRIMES[-1] == 313
    assert all(w.is_prime(n) == bool(sv[n]) for n in range(2, w.LIMIT))


def test_wheel65_313_catches_only_two_composites_below_1e5():
    from rhp import wheel65 as w
    hit = [n for n in range(314, 100000) if w.ringing(n) and w.ringing(n)[0] == 313]
    assert hit == [97969, 99221]
    assert 317 * 317 > 99999


def test_wheel65_response_is_ramanujan_sum():
    import numpy as np
    from rhp import wheel65 as w
    n = np.array([2, 30, 313 * 317, 100003])
    assert list(w.response(n)) == [2 - 65, 2 + 3 + 5 - 65, 313 - 65 + 0, -65][:2] + [313 - 65, -65]


# ── the Heartbeat waveform from the two sets ──
def test_wave_partition_matches_engine_and_is_exact():
    from rhp import wave as w
    e = w.engine_two_trees(100_000)
    assert (e["telperion_prime"], e["laurelin_composite"], e["mingling_0_and_1"]) == w.own_partition(100_000) == (9592, 90407, 2)
    assert e["exact"]


def test_wave_trees_are_exact_mirrors_up_to_floor_term():
    import numpy as np
    from rhp import wave as w
    _, psi = w.lam_table(2_000_000)
    u = np.linspace(6, 14.5, 4000)
    d, b = w.antisymmetry_defect(u, psi)
    assert d <= b + 1e-12


def test_wave_zero_tones_rebuild_the_sieve_wave_and_power():
    import math
    import numpy as np
    from rhp import wave as w, diverge as dv
    _, psi = w.lam_table(2_000_000)
    Z = dv.load_zeros()
    u = np.linspace(6, 14.5, 6000)
    c100, _ = w.reconstruction_error(u, psi, Z[:100])
    c1000, r1000 = w.reconstruction_error(u, psi, Z[:1000])
    assert c100 > 0.9 and c1000 > 0.97 and r1000 < 0.06
    wt, _ = w.trees(u, psi)
    assert abs(w.wind_stats(wt)["rms"] - math.sqrt(w.MEAN_POWER)) < 0.02


def test_wave_factor_step_peaks_at_e_squared():
    import math
    from rhp import wave as w
    best = max(range(2, 200), key=lambda p: w.factor_step_ratio(p, 1.0)[0])
    assert best == 7                                    # nearest integer to e² = 7.389
    assert abs(w.factor_step_ratio(math.e ** 2, 1.0)[0] - 2 / math.e) < 1e-12


def test_wave_tone_is_phase_shifted_sine_of_amplitude_two_over_rho():
    import numpy as np
    from rhp import wave as w
    g = 14.134725
    u = np.linspace(0, 20, 5000)
    t = w.tone(u, g)
    assert abs(t.max() - 2 / abs(complex(0.5, g))) < 1e-4
    ref = -2 * (0.5 * np.cos(g * u) + g * np.sin(g * u)) / (0.25 + g * g)      # −2[½cos + γ sin]/|ρ|²
    assert np.allclose(t, ref)


# ── the pathway of the first K primes ──
def test_primepath_survivors_after_65_primes_are_not_only_primes():
    import numpy as np
    from rhp import primepath as pp
    from rhp.wrap import sieve, primes_to
    assert abs(pp.survivor_density(65) / pp.mertens_estimate(313) - 1) < 0.02
    N = np.arange(1_000_000, 1_200_000); ok = np.ones(len(N), bool)
    for p in primes_to(313):
        ok &= (N % int(p) != 0)
    prime = sieve(1_200_000)[1_000_000:1_200_000]
    assert ok.mean() > 0.09 and ok.mean() < 0.10
    assert (ok & ~prime).sum() / ok.sum() > 0.2          # a quarter of the survivors are composites with all factors > 313


def test_primepath_prime_sum_follows_pnt_smooth_part_and_keeps_growing():
    import numpy as np
    from rhp import primepath as pp
    P = pp.first_primes(100_000)
    z = pp.prime_sum_path(3.0, P)
    sm = pp.smooth_prime_sum(3.0, P)
    assert abs(abs(z[-1]) / abs(sm[-1]) - 1) < 0.02
    assert abs(z[-1]) > 5 * abs(z[999])                  # larger at 10⁵ primes than at 10³


def test_primepath_wheel_sphere_radius_grows_like_sqrt_lnln():
    import math
    from rhp import primepath as pp
    r65, r1e4 = pp.wheel_sphere_radius(pp.first_primes(65)), pp.wheel_sphere_radius(pp.first_primes(10_000))
    assert abs(r65 - 1.42215) < 1e-4 and abs(r1e4 - 1.64598) < 1e-4
    assert abs(r1e4 ** 2 - (math.log(math.log(104729)) + 0.2615)) < 0.05      # Mertens: Σ1/p ≈ ln ln p + M


def test_primepath_zero_heights_need_about_two_decimals():
    import math
    import numpy as np
    from rhp import wave as w, diverge as dv
    Z = dv.load_zeros(1000)
    _, psi = w.lam_table(1_000_000)
    u = np.linspace(6, 13.5, 3000)
    wt, _ = w.trees(u, psi)
    corr = lambda Zr: float(np.corrcoef(wt, w.heartbeat(u, Zr) - math.log(2 * math.pi) * np.exp(-u / 2))[0, 1])
    c8, c2, c0 = corr(Z), corr(np.round(Z, 2)), corr(np.round(Z, 0))
    assert c8 - c2 < 0.01 and c0 < 0.5


def test_primepath_primes_are_spectral_lines_of_log_zeta():
    import math
    from rhp import primepath as pp
    pk = pp.spectral_lines(sigma=1.2, T=150.0, dt=0.05, umax=2.3, thresh=0.1)
    us = [u for u, a in pk]
    for n in (2, 3, 4, 5, 7):                                  # primes 2,3,5,7 and the prime power 4
        assert min(abs(u - math.log(n)) for u in us) < 0.04, n
    a2 = max(a for u, a in pk if abs(u - math.log(2)) < 0.04); a3 = max(a for u, a in pk if abs(u - math.log(3)) < 0.04)
    assert abs(a2 / a3 - (2 ** -1.2) / (3 ** -1.2)) < 0.25    # heights follow p^{−σ} (finite-window leakage)


# ── Two Trees inside ζ (twotrees_zeta) and the Gemini sims (geminisim) ──
def test_twotrees_dirichlet_identities():
    from rhp import twotrees_zeta as tz
    r = tz.identity_residuals(complex(2, 3), N=100000)
    assert r["split"] < 1e-12 and r["convolution"] < 1e-12      # T + L = ζ − 1 ;  −ζ′ = ζ·T  (1∗Λ = ln)
    assert r["series_L"] < 1e-6                                   # Laurelin sum = closed form


def test_twotrees_line_heights_split_the_integers():
    from rhp import twotrees_zeta as tz
    lam = tz.vonmangoldt_table(40)
    h = tz.line_heights(sigma=1.2, T=150.0, dt=0.05, ns=(2, 3, 5, 6, 7, 10))
    unit = lambda n: n ** -1.2
    for n in (2, 3, 5, 7):                                        # primes: T, ln ζ and L all show the line
        assert abs(h["T"][n] / (lam[n] * unit(n)) - 1) < .1
    for n in (6, 10):                                             # composites: ζ and L show it, T and ln ζ do not
        assert h["zeta"][n] > 0 and h["L"][n] / unit(n) > .8 and h["T"][n] < .06 * unit(n)


def test_geminisim_curves():
    from rhp import geminisim as G
    assert G.fib_ratio_error(12) < 1e-4
    (to, ro), (ti, ri) = G.transition()
    assert abs(ro[-1] * ri[0] - (1 + G.GAP) * ro[-1] * (1 / ro[-1])) < 0.02       # crossing at r = 1 is continuous up to the gap
    th, r = G.revised2("Heartbeat", 1.0)
    assert len(th) == 600 and abs(G.turns_per_index(1.0) - 1 / np.pi ** 2) < 1e-12
