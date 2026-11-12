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


def test_original_sim_worm_is_phase_modulation_rosette_is_amplitude_modulation():
    from scipy.special import jv
    from rhp import geminisim as G
    for t in (0.6, 1.0):                                  # whole cycles in the 10π window: lines are exactly Bessel J_n(5)
        assert max(abs(m - j) for _, j, m in G.worm_sidebands(t)) < 1e-9
    c, up, dn = G.rosette_sidebands(0.5)
    assert abs(c - 1) < 1e-3 and abs(up - .1) < 1e-3 and abs(dn - .1) < 1e-3     # carrier 1, sidebands 1 ± k at depth/2
    th, r = G.original_phase(2, 1.0)
    assert abs(th[-1] - 100 * np.pi) < 1e-9 and r.max() <= 1.2 * np.sqrt(100 * np.pi)


def test_prime_wheels_never_close_a_rosette():
    from rhp import geminisim as G
    cv = G.convergents_ratio(3, 2, 6)
    assert [(h, k) for h, k, e in cv][-3:] == [(3, 2), (8, 5), (19, 12)] or (19, 12) in [(h, k) for h, k, e in cv]
    assert all(e > 0 for h, k, e in cv)                    # ln3/ln2 is irrational: no convergent is exact


# ── Boundary Condition Engine (rhp/boundary.py, rhp/boundary_engine.py) ──
def test_boundary_engine_menu_drives_the_map_layer():
    from rhp import boundary_engine as E
    log = E.selftest(verbose=False)
    assert log and all(s == "PASS" for s, _ in log) and len(log) > 40


def test_worm_turnaround_law_is_exact():
    from rhp import boundary as B
    for tau in (.3, .6, 1.0):
        tb = np.linspace(0, 2 * np.pi / tau, 200001)                       # one full wobble cycle
        rate = 1 + 5 * tau * np.cos(tau * tb)                              # dθ/dθ_b of the worm, ω = τ, β = 5
        assert abs(np.mean(rate < 0) - B.retro_fraction(5, tau)) < 2e-3
    assert B.retro_fraction(5, .15) == 0 and B.worm_retro_threshold() == .2


def test_wheel_transition_turnarounds_stop_at_beta_omega_one():
    from rhp import boundary as B
    al, c = B.transition_curve("O2", "O3", "wheel", m=21, n=3000)
    thr = 1 - np.sqrt(1 / (5 * B.PATTERNS["O2"].tau0))                       # β(α)·ω(α) = 5τ(1−α)² = 1
    assert c["angular_reversals"][al < thr - .05].min() > 0 and c["angular_reversals"][al > thr + .05].max() == 0


def test_zeta_turns_one_way_and_wobble_decays():
    from rhp import boundary as B
    zt = B.zeta_turning(B.PP.GAMMA1, 5000)
    assert zt["last_sense_flip_prime"] < 100                                  # beyond the first few primes every step turns the same way
    meas, pred = B.zeta_wobble()
    assert abs(meas / pred - 1) < .15                                         # β_ζ(p) = ½ p^{-½}


def test_uclock_zeta_loops_uniform_in_u_and_log_spiral_radius():
    from rhp import boundary as B
    L = B.zeta_loops(B.PP.GAMMA1, 60000)
    assert abs(L["spacing_u"].mean() / L["expected_spacing"] - 1) < .01 and L["spacing_u"].std() / L["spacing_u"].mean() < .06   # uniform in u
    assert L["spacing_index"][-1] > 100 * L["spacing_index"][0]                                                                  # wildly non-uniform in index
    L5 = B.zeta_loops(5.0, 60000)
    assert abs(np.median(L5["growth"]) / np.median(L5["smooth_growth"]) - 1) < .06                                              # away from zeros: the smooth spiral
    assert np.allclose(B.worm_loop_growth(), [1.25, 1.33, 1.49, 1.95, 19.0], rtol=.03)


def test_zeros25k_recreation_in_house():
    import json, os
    p = "data/zeros_25000_hp.tsv"
    if not os.path.exists(p): return
    g = np.array([float(l.split("\t")[1]) for l in open(p).read().split("\n")[1:] if l.strip()])
    assert len(g) == 25000 and np.all(np.diff(g) > 0.03)                       # strictly increasing, min gap 0.0353
    assert abs(g[0] - 14.134725141734693790) < 1e-15 and abs(g[24999] - 21942.5924301325946845) < 1e-12
    r = json.load(open("experiments/zeros25k/results.json"))
    assert r["A"]["failures"] == 0 and r["B"]["count_wrong"] == 0              # every root a sign change; Backlund count = k at 500 heights
    assert r["A"]["table_err_max"] < 1e-13                                      # the SemPriCon table is good to 5e-14


def test_sedenion_slot_machinery_and_norm_identity():
    import os, json
    from rhp import sedenion_slots as S
    for N, nd in ((16, 42), (32, 294), (64, 1518)):
        A = S.Alg(N); I = np.eye(N)
        assert all(np.allclose(A.Lmat(k) @ A.Lmat(k), -I) for k in range(1, N))        # every slot is a circle
        D = S.deform_matrix(N); assert int(D.sum() // 2) == nd and (N - 1) * (N - 2) // 2 - nd == {16: 63, 32: 171, 64: 435}[N]
    P = S.primes_upto(2000)
    assert (P[63], P[64], P[127], P[128]) == (311, 313, 719, 727)                      # 64th, 65th, 128th, 129th prime
    pk = P[:120]; t = 14.134725
    p1 = S.euler_path(S.rule_slots("same", pk, 16), pk, t, 16); z = S.complex_path(pk, t)
    assert np.abs(p1[:, 0] + 1j * p1[:, 1] - z).max() < 1e-12                          # control: one slot = complex zeta exactly
    rng = np.random.default_rng(3)
    for N in (16, 64):                                                                   # norm identity: any imaginary-slot assignment has |v_K| = |zeta_K|
        pr = S.euler_path(rng.integers(1, N, len(pk)), pk, t, N)
        assert np.allclose(np.linalg.norm(pr, axis=1), np.abs(z), rtol=1e-12)
    sib = os.path.join(os.path.dirname(__file__), "..", "..", "SedenionSpectralRelativity")
    if os.path.exists(os.path.join(sib, "emerger_spectrum.py")):
        import sys; sys.path.insert(0, sib); import emerger_spectrum as E
        S16 = S.sign_table(16)
        assert all(int(E.cd_mul(E.e(i), E.e(j))[i ^ j]) == S16[i, j] for i in range(16) for j in range(16))


def test_sedenion_slot_rule_results_recorded():
    import os, json
    p = "experiments/sedenion_slot_rule.json"
    if not os.path.exists(p): return
    r = json.load(open(p))
    for N in ("16", "32", "64", "128"):
        for rule in ("seq0", "seq1"):                                                    # the sequential rules are not distinguished from random assignments
            assert all(abs(v["z"]) < 3.5 for v in r["E1"][N][rule].values())


def test_ping_to_200k_recorded():
    import os, json
    p = "experiments/zeros25k/ping_accuracy_n200000.json"
    if not os.path.exists(p): return
    r = json.load(open(p)); a = lambda K, b: r[K]["blocks"][b]["balanced_acc"]
    assert a("3066", "300-999") == 1.0 and a("3066", "1000-1999") > 0.99 and a("3066", "2000-2999") < 0.96      # perfect to ~10³ at 3066 zero-tones
    assert a("25000", "6000-7999") == 1.0 and a("25000", "8000-9999") > 0.995 and a("25000", "20000-29999") < 0.9   # perfect to ~10⁴ at 25,000
    assert all(v["cohen_d"] < 0 for K in r for v in r[K]["stragglers"].values())                                   # stragglers read less prime-like, not more


def test_slot0_empty_keeps_the_euler_norm_exact():
    import os, json
    p = "experiments/sedenion_slot_rule.json"
    if not os.path.exists(p): return
    r = json.load(open(p))
    for N in ("16", "32", "64", "128"):
        for rule, n0 in r["slot0_primes"][N].items():
            rl = r["E2"][N][rule]["ratio_log"]["value"]
            assert (n0 == 0) == (rl < 1e-12), (N, rule, n0, rl)     # |v_K| = |zeta_K| exactly iff no prime sits on the identity slot


def test_nodal_extractor_core_and_results():
    import os, json, math
    from rhp import nodal_extract as NX
    x = np.linspace(-2, 2, 401); X, Y = np.meshgrid(x, x)
    L, c = NX.zero_set_2d(x, x, X ** 2 + Y ** 2 - 1)
    assert len(L) == 1 and c == [True] and abs(NX.polyline_length(L[0]) - 2 * math.pi) < 3e-4                  # one closed curve, the unit circle
    xs = np.linspace(-2, 2, 321); XX, YY = np.meshgrid(xs, xs); Zc = (XX + 1j * YY) ** 2 - 1
    Rl, _ = NX.zero_set_2d(xs, xs, Zc.real); Il, _ = NX.zero_set_2d(xs, xs, Zc.imag)
    pts = NX.seg_intersections(Rl, Il); assert len(pts) == 2 and np.allclose(sorted(pts[:, 0]), [-1, 1], atol=2e-3)      # zeros of z² − 1
    psi, lam = NX.psi_table(400); ex, dr, kd = NX.exact_prime_nodes(100, 300, psi)
    xv = []; 
    for n in range(100, 301): xv += [n - 1e-9, n, n + 0.25, n + 0.5, n + 0.75, n + 1 - 1e-9]
    brute = 0; prev = None
    for xx in xv:
        v = psi[int(math.floor(xx))] - xx; s = 1 if v > 0 else -1
        if prev is not None and s != prev: brute += 1
        prev = s
    assert abs(len(ex) - brute) <= 2                                                                              # exact nodes agree with brute sampling
    p1 = "experiments/nodal_extract/r1.json"; p2 = "experiments/nodal_extract/r2.json"
    if os.path.exists(p1):
        r = json.load(open(p1))
        assert r["1000-1500"]["exact_nodes"] == 39 and r["1000-1500"]["per_K"]["25000"]["exact_matched"] > .99 and r["3000-3500"]["per_K"]["25000"]["exact_matched"] > .99
        assert r["1000-1500"]["per_K"]["25000"]["nodes"] < r["1000-1500"]["per_K"]["25000"]["rice_predicted"] / 2   # Rice overpredicts at large K
    if os.path.exists(p2):
        r = json.load(open(p2))
        assert r["rungs_found"] == r["zeros_known_below_tmax"] == 21 and r["unmatched_known"] == 0 and r["offline_crossings"] == 0
        assert r["re_domains_pos"] + r["re_domains_neg"] == r["rungs_found"] + 1 and r["im_lines_offline"] == 20


def test_fermat_side_probe_recorded():
    import os, json, math
    p = "experiments/fermat_side_probe.json"
    if not os.path.exists(p): return
    r = json.load(open(p))
    assert abs(r["L_chi_m4_zeros_below_100"][0] - 6.0209489) < 1e-6 and len(r["L_chi_m4_zeros_below_100"]) == 50
    assert r["r2_over_4_n1_32"][:10] == [1, 1, 0, 1, 2, 0, 0, 1, 1, 2]                              # sums of two squares: n = 5 and 10 have 2 Gaussian representations
    assert abs(r["LambdaK"]["5"] - 2 * math.log(5)) < 1e-12 and r["LambdaK"]["3"] == 0 and abs(r["LambdaK"]["9"] - 2 * math.log(3)) < 1e-12
    assert all(v["max_ratio"] <= 1.0 for v in r["weil_bound"].values())                              # Weil: Riemann hypothesis for the Fermat curves over F_p
    assert r["two_squares_density"]["1000000"]["density"] < r["two_squares_density"]["1000"]["density"]


def test_nshape_shape():
    import os, json, math
    from rhp import nshape as NS
    om = NS.big_omega(1000); assert om[1:21].tolist() == [0, 1, 1, 2, 1, 2, 1, 3, 2, 2, 1, 3, 1, 2, 2, 4, 1, 3, 1, 3]
    assert [NS.e8_shell(m) for m in (2, 4)] == [240, 2160] and [NS.count_hurwitz_norm(n) for n in (1, 2, 3, 5)] == [24, 24, 96, 144]
    assert [NS.euler_factor(N) for N in (65, 221, 1073)] == [5, 17, 29]
    P = NS.primes_upto(2000); isp = np.zeros(2001, bool); isp[P] = True; odd = np.arange(2001) % 2 == 1
    h = NS.hyperbola_points(2000); assert ((h == (odd & ~isp & (np.arange(2001) > 1)))[odd & (np.arange(2001) > 1)]).all()
    p = "experiments/nshape_shape.json"
    if os.path.exists(p):
        r = json.load(open(p))
        assert r["hyperbola"]["lattice_reach_equals_odd_composites"] and r["hyperbola"]["direct_search_ok"]
        assert r["circle"]["theory_matches_all_odd_N"] and r["circle"]["euler_factor_failures"] == 0 and r["circle"]["all_2plus_are_composite"]
        assert abs(r["tower"]["hole_fraction_among_primes_below_1e6"]["Zi"] - 0.5) < 0.01 and r["onion"]["10000000"]["mode"] == 3


def test_bitflip():
    import os, json
    from sympy import jacobi_symbol, factorint
    from rhp import bitflip as BF
    rng = np.random.default_rng(5); n = rng.integers(1, 10 ** 6, 4000) * 2 + 1; a = rng.integers(-10 ** 5, 10 ** 5, 4000)
    assert (BF.jacobi(a, n) == np.array([jacobi_symbol(int(x) % int(m), int(m)) for x, m in zip(a, n)])).all()
    assert (BF.modpow(a, (n - 1) // 2, n) == np.array([pow(int(x), (int(m) - 1) // 2, int(m)) for x, m in zip(a, n)])).all()
    P = np.array([p for p in range(3, 3000) if all(p % q for q in range(2, int(p ** .5) + 1))], np.int64)
    assert all(BF.two_ways(P, d)[2].all() for d in BF.DISC)                       # primes agree on every base (Euler's criterion)
    for c in (1857241, 2433601, 6189121):                                        # the three Carmichael survivors of all seven bits below 10⁷
        f = factorint(c); assert all(BF.two_ways(np.array([c]), d)[2][0] for d in BF.DISC) and all((c - 1) % (q - 1) == 0 for q in f)
    p = "experiments/bitflip_test.json"
    if os.path.exists(p):
        r = json.load(open(p)); assert r["primes_agree_all_7"] == r["primes_total"] and r["composites_agreeing_all_7"] == 3
        assert r["per_base_composites_agreeing"]["-4"] == 750 and r["minimal_subset_zero_survivors_below_100000"][0] == 2


def test_zdtree():
    import os, json
    from rhp import zdtree as T, sedenion_slots as SS
    prev = None
    for k, zd in ((4, 42), (5, 294), (6, 1518)):
        n = 2 ** k; Z, per, tot = T.zd_incidence(k); carry, clean = T.classes(per, tot)
        assert int(Z[np.triu_indices(n, 1)].sum()) == zd and len(carry) == 2 ** k - k - 5
        assert clean == list(range(1, 8)) + [2 ** j for j in range(3, k)]                                  # octonion differences + pure doublings
        assert (Z[1:, 1:] == SS.deform_matrix(n)[1:, 1:]).all()                                              # ZD incidence = anticommutation graph
        if prev is not None: assert (Z[: n // 2, : n // 2][1:, 1:] == prev[1:, 1:]).all()                    # old pairs never change status
        prev = Z
    p = "experiments/zdtree_probe.json"
    if os.path.exists(p):
        r = json.load(open(p)); assert r["8"]["zd_pairs"] == 29886 and r["8"]["carrying_delta"] == 243 and all(r[k].get("survivor_recursion_ok", True) for k in r)


def test_normcheck():
    import os, json
    from rhp import normcheck as nc
    for N in (1, 2, 4, 8): assert nc.level_report(N, 400)["all"]["max_abs"] < 1e-12                       # Hurwitz: composition exact
    r16 = nc.level_report(16, 800)
    assert r16["all"]["max_abs"] > 0.2 and all(r16[k]["max_abs"] < 1e-12 for k in ("LL", "LU", "UL", "UU"))   # fails only for factors that mix the halves
    assert nc.level_report(32, 400)["LL"]["max_abs"] > 0.2                                                  # T_32's lower half is T_16: the failure is inherited
    v = nc.basis_pair_sums(16); assert set(v) == {0, 4, 8} and v[0] == v[8] == 336                          # norm 0 and norm-doubling in equal number
    assert nc.zero_product_halves(16) == {("LU", "LU"): 336}


def test_j2():
    import numpy as np
    from rhp import j2, sedenion_slots as SS
    for N in (8, 16, 32):
        Jl, Jp = j2.J_left(N), j2.J_perm(N)
        assert j2.order(Jl) == 4 and j2.projective_order(Jl) == 2 and j2.swaps_colour(Jl, N) and np.allclose(Jl @ Jl, -np.eye(N))     # swap by multiplication: a 4-cycle, an involution on rays
        assert j2.order(Jp) == 2 and j2.swaps_colour(Jp, N) and int(np.nonzero(Jp[:, 0])[0][0]) == 4                                   # the strict involution moves 1 (e₀) to e₄: not an automorphism
        assert min(j2.automorphism_defects(Jl, N, 40)) > 1
    assert j2.orbit_of_identity(16, 4, 4) == [(0, 1), (4, 1), (0, -1), (4, -1), (0, 1)]
    assert j2.klein_orbit_size(0.5 + 14.13j) == 2 and j2.klein_orbit_size(0.7 + 14.13j) == 4
    P = SS.primes_upto(5000)[1:]
    for m in (1, 2, 3): assert j2.decomposition_check(P, SS._DISC[:m])[1] == 0
    sh = j2.rule_class_shares(16, "charx", SS.primes_upto(20000)[:500]); assert sh["zero"] == 0 and sh["deform"] == 0 and sh["anti"] > 0.8
    assert j2.rule_class_shares(16, "char", SS.primes_upto(20000)[:500])["zero"] > 0.03


def test_beurling():
    import os, json
    from rhp import beurling as B
    from rhp.sedenion_slots import primes_upto
    P = primes_upto(2_000_000); xs = np.exp(np.linspace(np.log(1e5), np.log(2e6), 60))
    one = np.ones(len(P), np.int64)
    assert np.abs(B.psi_beurling(P, one, xs) / xs - 1).max() < 0.005                                      # true primes: ψ/x → 1
    m = B.gated_multiplicity(P, 1.0); r = B.psi_beurling(P, m, xs) / xs
    assert np.abs(r - B.gated_prediction(xs, 1.0)).max() < 0.02 and r.max() - r.min() > 1.0             # gated: PNT fails, oscillation as predicted from the harmonics
    a = B.generalized_integers(P[P <= 2000], m[P <= 2000], 2000); assert (a[1], a[2], a[3], a[4], a[5]) == (1, 2, 2, 3, 0)   # m(2) = m(3) = 2, m(5) = 0: (1−x)^−2 gives 3 products equal to 4
    assert B.z5_irreducibles_brute(300) == B.z5_irreducibles(300)[2] == 83                              # irreducibles of Z[√−5]
    n_all, n_pr, n_irr = B.z5_irreducibles(10 ** 6); assert n_irr > 1.5 * n_all
    p = "experiments/beurling_probe.json"
    if os.path.exists(p):
        r = json.load(open(p)); assert 0.19 < r["true"]["rms_windows"][0] < 0.25 and abs(r["true"]["zero_sum_prediction"] - 0.215) < 0.01


def test_mertens_wave():
    import os, json
    from rhp import mertens_wave as W
    g, zp = W.load_zp(K=3000); c = W.residues(g, zp); mu = W.mobius_table(5000); M = np.cumsum(mu.astype(np.int64))
    assert mu[1] == 1 and mu[4] == 0 and mu[6] == 1 and mu[30] == -1
    x = np.array([100.5, 1000.5]); assert np.abs(W.m_reconstruct(x, g, c) - M[[100, 1000]]).max() < 0.7           # the explicit formula
    n = np.arange(2, 600); assert np.mean(W.classify(W.mu_ping(n, g, c, 3000)) == mu[n]) == 1.0                    # 3000 tones: μ perfect to ≈ 600
    assert abs(W.convolution_residual(200, g, c, 3000)) < 0.01                                                      # Σ M(x/n) = 1
    p = "experiments/mertens_probe.json"
    if os.path.exists(p):
        r = json.load(open(p)); d = r["exact"]["decades"]; assert all(v["corr_M_psi"] < -0.7 for v in d.values())   # the two waves are in antiphase
        assert abs(r["slope"]["partial_variance"]["25000"] ** 0.5 - 0.17) < 0.01 and r["exact"]["max_M_over_sqrt"] < 1.0


def test_zerophase():
    import math, os, json
    from rhp import zerophase as Z
    u = np.linspace(0.5, 6, 1500); g1 = 14.134725
    b = Z.bearing(Z.lift_stereo(Z.zero_spiral(u, g1)))
    assert b.std() < 1e-3 and abs(math.degrees(b.mean()) - math.degrees(math.atan(2 * g1))) < 1e-3                  # loxodrome: constant bearing arctan 2γ
    assert np.degrees(Z.bearing(Z.lift_linear(u, g1 * u, 6.0)).std()) > 0.5                                         # the linear lift is not conformal
    assert abs(Z.pitch_deg(g1) - 2.026) < 0.01 and abs(Z.polygon_h(4) - 1) < 1e-12 and abs(Z.polygon_h(6) - 2 ** -0.5) < 1e-12
    w, lx = Z.weights(20000); ts = np.array([14.134725, 17.5, 21.022040]); A = Z.scan_zero(ts, w, lx)
    assert A[0] > 0.7 and A[2] > 0.7 and A[1] < 0.4                                                                 # lines at the zeros, quiet between
    p = "experiments/zerophase_probe.json"
    if os.path.exists(p):
        r = json.load(open(p)); assert r["zero_scan"]["zeros_in_range"] == 29 and r["zero_scan"]["peak_height_min"] > 0.8 and r["zero_scan"]["max_offset"] < 0.1


def test_omega_dstar():
    import mpmath as mp
    from rhp import omega_dstar as O
    r = O.omega_routes(40); v = [r[k] for k in ("lambertw", "iteration", "newton", "bisection_digit_curves")]
    assert max(abs(a - b) for a in v for b in v) < mp.mpf(10) ** -35                                         # four routes, one number
    d = O.d_star_faces(30); assert abs(d["d_star"] - d["d_star_equals_minus_log10_omega"]) < mp.mpf(10) ** -25 and d["closure_residual"] < mp.mpf(10) ** -25
    assert abs(float(d["d_star"]) - 0.246307201473) < 1e-11 and abs(float(d["gap_in_omega_units"]) - 7.0736e-4) < 1e-7
    m = O.map_pair(); assert abs(float(m["product"]) - 1) < 1e-6 and m["forward_signs"][:4] == [1, -1, 1, -1]  # multipliers −Ω and −1/Ω, alternating
    n = O.n_shape_crossing(); assert float(n["closes"]) < 1e-25                                                # the N-shape closes at E = 2Ω
    w, wu = O.lambert_branches(); assert abs(float(w["W0_fold"]) + 1) < 1e-13 and abs(float(w["Wm1_fold"]) + 1) < 1e-13   # square-root singularity at −1/e: ~√ε accuracy, not ε
    assert O.boltzmann_orientation()["root_exists"] is False


def test_exp_cardioid():
    import mpmath as mp
    from rhp import omega_dstar as O
    om = mp.lambertw(1).real
    a = O.exp_cardioid(-om); assert abs(a["c"] + 1) < mp.mpf(10) ** -25 and abs(a["multiplier"] + om) < mp.mpf(10) ** -25      # c = −1: the entropy map, multiplier −Ω
    b = O.exp_cardioid(1); assert abs(b["c"] - 1 / mp.e) < mp.mpf(10) ** -25 and abs(b["multiplier"] - 1) < mp.mpf(10) ** -25   # the cusp: fold
    d = O.exp_cardioid(-1); assert abs(d["c"] + mp.e) < mp.mpf(10) ** -25 and abs(d["multiplier"] + 1) < mp.mpf(10) ** -25      # the flip
    z = O.exp_cardioid(mp.e ** (1j * 0.7)); assert abs(abs(z["multiplier"]) - 1) < mp.mpf(10) ** -25                           # the whole boundary has |multiplier| = 1


def test_lambda_family():
    import mpmath as mp
    from rhp import omega_dstar as O
    f = O.lambda_family(); om = mp.lambertw(1).real
    assert abs(f["paper"]["x_star"] - om) < mp.mpf(10) ** -25 and abs(f["paper"]["D_star"] + om / mp.log(10)) < mp.mpf(10) ** -25          # λ = 1: Ω, D* = −d*
    assert abs(f["flip"]["x_star"] - 1 / mp.e) < mp.mpf(10) ** -25 and abs(f["flip"]["mult"] + 1) < mp.mpf(10) ** -25                       # λ = e: x* = 1/e, multiplier −1
    assert abs(f["fold"]["x_star"] - mp.e) < mp.mpf(10) ** -12 and abs(f["fold"]["mult"] - 1) < mp.mpf(10) ** -12                          # λ = −1/e: tangency at x = e (sqrt-eps accuracy)
    assert abs(f["vacuum (limit)"]["x_star"] - 1) < mp.mpf(10) ** -15 and abs(f["vacuum (limit)"]["mult"]) < mp.mpf(10) ** -15
    assert abs(float(O.d_star_faces(30)["d_star"]) / float(mp.log10(mp.e)) - float(om)) < 1e-12                                              # d*/log₁₀e = Ω: the fraction of the way from vacuum to flip
    k, c = O.bulb_tangent_point(1, 4); assert abs(k - 1j) < 1e-25 and abs(c - 1j * mp.e ** (-1j)) < 1e-25
