# SPDX-License-Identifier: GPL-3.0-only
"""Ω = W(1) and d* reproduced from the two curves of the First Age experiment (2026-04-09), blind, then ignored.   :CALCULATED (mpmath, 50 digits) unless marked.

Sources read (2026-11-11): Ainulindale/AgeFirst/original_2026/ResearchPaper/SMNNIP_Omega_Constant_Paper.txt, SMNNIP_Maths_Notes_Log.txt;
Ainulindale/wiki/17 (the four d*), 50 (causality ceiling), 73 §4 (N-shape: J_f = e^{−σE}, J_b = e^{−(1−σ)E}, J_f·J_b = e^{−E}), 96 (the Horizon); ValaQuenta modules/constants `derive_omega_zs`, `derive_d_star`.

The experiment as written: x = T/T_Planck.
  forward entropy ceiling  S_f(x)/k_B = e^{−x}     (falling)    in decimal digits  D_f(x) = log₁₀ e^{−x} = −x/ln 10
  inverse entropy floor    S_i(x)/k_B = x          (rising)                        D_i(x) = log₁₀ x
  they cross once, where x = e^{−x}  ⇒  x* = Ω,  T* = Ω T_Planck,  and the digit coordinate of the crossing is  D* = log₁₀ Ω = −Ω/ln 10 = −d*.
The step S_f = e^{−x} is ASSERTED in the paper ('substituting the Bekenstein entropy'), not derived.  Everything below follows from it.
"""
from __future__ import annotations
import mpmath as mp

T_PLANCK = 1.416784e32          # K, as quoted in the First Age paper
D_STAR_SPEC = 0.24600           # the measured (SPARC / Berry–Keating) value quoted in wiki 17; not reproducible here


def omega_routes(dps=50):
    """Ω by four routes that share no step: closed form W(1); the contraction x → e^{−x}; Newton on x e^x − 1; bisection on the digit curves D_f = D_i (no W, no e^x solved)."""
    mp.mp.dps = dps; out = {}
    out["lambertw"] = mp.lambertw(1).real
    x = mp.mpf("0.5"); n = 0
    while True:
        y = mp.e ** (-x); n += 1
        if abs(y - x) < mp.mpf(10) ** (-(dps - 2)): break
        x = y
    out["iteration"] = y; out["iteration_steps"] = n
    x = mp.mpf("0.5")
    for _ in range(40): x = x - (x * mp.e ** x - 1) / ((x + 1) * mp.e ** x)
    out["newton"] = x
    Df = lambda t: -t / mp.log(10); Di = lambda t: mp.log10(t)
    a, b = mp.mpf("0.01"), mp.mpf(1)
    for _ in range(int(dps * 3.4) + 10):
        m = (a + b) / 2
        if (Df(m) - Di(m)) > 0: a = m               # at small x the floor log₁₀x is far below the ceiling; the ceiling is above until the crossing
        else: b = m
    out["bisection_digit_curves"] = (a + b) / 2
    return out


def d_star_faces(dps=50):
    """The four faces of wiki 17 from Ω: d*, d*_taut, d*_ln10 (the Translator) and the closure d*·ln10 = Ω."""
    mp.mp.dps = dps; om = mp.lambertw(1).real; ln10 = mp.log(10)
    dtaut = om / ln10
    return dict(omega=om, d_star=dtaut, d_star_equals_minus_log10_omega=-mp.log10(om), d_star_ln10=dtaut * ln10, closure_residual=abs(dtaut * ln10 - om),
                d_star_spec=D_STAR_SPEC, gap_in_omega_units=abs(D_STAR_SPEC * ln10 - om), gap_in_d_units=abs(D_STAR_SPEC - dtaut))


def crossing_temperature(dps=30):
    """T* = Ω T_Planck and the digit coordinates quoted in the paper (log₁₀ T* = 31.905…); the 'match' |D_max − D_min| in float64 and at 30 digits."""
    mp.mp.dps = dps; om = mp.lambertw(1).real; T = om * mp.mpf(T_PLANCK)
    dmin = mp.log10(T); dmax = mp.log10(mp.e ** (-om) * mp.mpf(T_PLANCK))
    import math
    o = float(om); f_match = abs(math.log10(o * T_PLANCK) - math.log10(math.exp(-o) * T_PLANCK))
    return dict(T_star=T, D_min=dmin, D_max=dmax, match_mp=abs(dmax - dmin), match_float64=f_match)


def map_pair(n=30):
    """The forward map f(x) = e^{−x} (entropy: falls toward the ceiling) and its inverse f⁻¹(x) = −ln x (inertia: runs away from it), both fixed at Ω.
    Multipliers f'(Ω) = −Ω and (f⁻¹)'(Ω) = −1/Ω: product 1; both alternate sides every step (period-2 approach / departure).  Returns the iterates and ratios."""
    mp.mp.dps = 40; om = mp.lambertw(1).real
    x = mp.mpf("0.6"); fw = [x]
    for _ in range(n): x = mp.e ** (-x); fw.append(x)
    y = om + mp.mpf("1e-12"); bw = [y]
    for _ in range(12): y = -mp.log(y); bw.append(y)
    r_f = [(fw[i + 1] - om) / (fw[i] - om) for i in range(n - 1)]; r_b = [(bw[i + 1] - om) / (bw[i] - om) for i in range(len(bw) - 1)]
    return dict(omega=om, forward_ratio_last=r_f[-1], backward_ratio_first=r_b[0], product=r_f[-1] * r_b[0], forward_signs=[int(mp.sign(v - om)) for v in fw[:10]],
                forward_ratio_theory=-om, backward_ratio_theory=-1 / om)


def n_shape_crossing(dps=30):
    """wiki 73 §4: J_f(σ) = e^{−σE}, J_b(σ) = e^{−(1−σ)E}, product e^{−E} for every σ; they cross at σ = ½ at height e^{−E/2}.  The crossing height equals the half-energy E/2
    (J = E/2: the N-shape closes on itself) exactly at E = 2Ω.  READING: the same fixed-point equation x = e^{−x} with x = E/2."""
    mp.mp.dps = dps; om = mp.lambertw(1).real; E = 2 * om
    return dict(E=E, J_cross=mp.e ** (-E / 2), half_E=E / 2, product=mp.e ** (-E), omega_squared=om ** 2, closes=abs(mp.e ** (-E / 2) - E / 2))


def lambert_branches(dps=30):
    """W₀(0) = 0, W₀(1) = Ω, W₀(−1/e) = W₋₁(−1/e) = −1 (the fold), and the same divided by ln 10 (the UNS reading)."""
    mp.mp.dps = dps; ln10 = mp.log(10)
    v = dict(W0_0=mp.lambertw(0).real, W0_1=mp.lambertw(1).real, W0_fold=mp.lambertw(-1 / mp.e).real, Wm1_fold=mp.lambertw(-1 / mp.e, -1).real)
    return v, {k: x / ln10 for k, x in v.items()}


def boltzmann_orientation(dps=30):
    """Check on the asserted step.  For a Planck-energy excitation the Boltzmann weight at temperature T is e^{−T_P/T} = e^{−1/x}, not e^{−x}.  Does x = e^{−1/x} have a real root?
    (ln x = −1/x ⇒ x ln x = −1, but min over x>0 of x ln x is −1/e > −1.)  Returns (min of x ln x, whether a root exists)."""
    mp.mp.dps = dps; xmin = 1 / mp.e; return dict(min_x_ln_x=xmin * mp.log(xmin), root_exists=bool(xmin * mp.log(xmin) <= -1))
