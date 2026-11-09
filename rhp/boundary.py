# SPDX-License-Identifier: GPL-3.0-only
"""Boundary Condition Engineering, made measurable: every pattern in the five Gemini event-horizon files, the transitions between them,
and the turn-around diagnostics used to compare them with the ζ pathway.

Sources (read-only; Ainulindale/AgeFirst/original_2026/ResearchPaper/):
  [O]  event_horizon_sim.py                       the ORIGINAL (Gemini, 143 lines): Orientation / Metric swap (worm) / Fermat frequency (rosette) / primes
  [C]  FirstDraft/FirstAge/Gemini/event_horizon_sim.py   Claude Sonnet, April 2026: 5 phases, golden-angle prime spiral, horizon invariance
  [G]  Gemini/Gemini_event_horizon_sim.py         z ↦ 1 + 1/z drawn in polar
  [R]  Gemini/Gemini_Revised_2_EHS.py             Inertial / Metric Swap / Heartbeat / Star with a LERP TransitionManager
  [X]  Gemini/Gemini_Transition_EHS.py            inertial path crossing r = 1 and inverting
Every pattern is a function (params, τ) → (θ, r) with τ ∈ (0, 1] the file's own animation clock.  [Z] patterns are the ζ side (this repo).
Labels: transcriptions are ESTABLISHED as to what the files do; the diagnostics and comparisons are OURS-CALCULATED; nothing here is evidence about RH.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from . import primepath as PP, twotrees_zeta as TZ, geminisim as G

PHI = G.PHI; INV_PHI = PHI - 1; PI_INV = 2 / np.pi; E_INV = 2 / np.e; GOLD_ANG = 2 * np.pi / PHI ** 2
_PR = G._primes(9000)                       # first 9000 primes, shared
_PRIMES_8000 = _PR[_PR <= 8000]
_PRIMES_10000 = _PR[_PR <= 10000]


def _unwrap_polar(z):
    return np.unwrap(np.angle(z)), np.abs(z)


# ── [O] the original ─────────────────────────────────────────────────────────────────────────────
def o_orientation(p, tau):
    return G.original_phase(0, tau, pi_f=p["pi_f"])[0], G.original_phase(0, tau, pi_f=p["pi_f"])[1]


def o_worm(p, tau):
    return G.original_phase(1, tau, e_f=p["e_f"])


def o_rosette(p, tau):
    return G.original_phase(2, tau, pi_f=p["pi_f"])


def o_primes(p, tau):
    pr = _PRIMES_10000[_PRIMES_10000 < 10000 * max(tau, 0.02)]
    return pr * p["pi_f"], np.sqrt(pr) * p["e_f"]


# ── [C] Claude's five phases ─────────────────────────────────────────────────────────────────────
def c_inertial(p, tau):
    th = np.linspace(0, 8 * np.pi, 3000); r = np.linspace(12, .1, 3000) * (1 - .3 * tau)
    return th, r * (1 + .08 * np.sin(th * PI_INV * 4))


def c_metric(p, tau):
    th = np.linspace(0, 12 * np.pi, 2000); ind = INV_PHI * (.5 + .5 * tau)
    return th, (.5 + .5 * np.exp(-ind * th / (4 * np.pi))) * 12


def c_heartbeat(p, tau):
    th = np.linspace(0, 20 * np.pi, 2000); ph = th + .95 * th * (1 - .3 * tau)
    return th, np.sqrt(th + .1) * (1 + .3 * np.sin(ph))


def c_golden(p, tau):
    n = min(len(_PRIMES_8000), max(10, int(tau * len(_PRIMES_8000) * 1.2))); idx = np.arange(1, n + 1)
    return idx * GOLD_ANG, np.sqrt(_PRIMES_8000[:n])


def c_horizon(p, tau):
    n = min(len(_PRIMES_8000), max(20, int(tau * 500 + 50))); idx = np.arange(1, n + 1); pr = _PRIMES_8000[:n]
    pr_r = np.sqrt(pr); th = idx * GOLD_ANG
    z = (.5 + pr_r / pr_r.max() * .48) + 1j * np.sin(th) * (pr_r / pr_r.max()) * .8
    return np.unwrap(np.angle(z - .5)), np.abs(z - .5)          # drawn about the horizon line x = ½


# ── [G] [R] [X] ──────────────────────────────────────────────────────────────────────────────────
def g_phiflow(p, tau):
    z = G.phi_flow(max(2, int(2 + 12 * tau))); return np.linspace(0, np.pi, len(z)), z


def r_phase(name):
    def f(p, tau):
        return G.revised2(name, .1 + 9.9 * tau)
    return f


def x_cross(p, tau):
    (to, ro), (ti, ri) = G.transition(h=.5 + 1.5 * tau)
    return np.concatenate([to, ti]), np.concatenate([ro, ri])


# ── the wheel family: one state vector that contains the worm AND the rosette (and everything between) ───────
WHEEL_KEYS = ("Th", "beta", "om", "d", "k", "m", "a")


def wheel_draw(s, n=5000):
    """θ_b = Th·u, angle = θ_b + β sin(ω θ_b), radius = (1−m)·a/(10 − 9.9u + 0.01) + m·√θ_b·(1 + d sin(k θ_b)).
    β = 5, ω = τ, d = 0, m = 0 is the WORM; β = 0, d = 0.2, k = (2/π)τ, m = 1, Th = 100π is the ROSETTE."""
    u = np.linspace(0, 1, n); tb = s["Th"] * u
    ang = tb + s["beta"] * np.sin(s["om"] * tb)
    r = (1 - s["m"]) * s["a"] / (10 - 9.9 * u + .01) + s["m"] * np.sqrt(tb) * (1 + s["d"] * np.sin(s["k"] * tb))
    return ang, r


def worm_state(p, tau):
    return dict(Th=10 * np.pi, beta=5.0, om=tau, d=0.0, k=0.0, m=0.0, a=p["e_f"] * tau)


def rosette_state(p, tau):
    return dict(Th=100 * np.pi, beta=0.0, om=0.0, d=0.2, k=p["pi_f"] * tau, m=1.0, a=0.0)


def w_family(p, tau):
    return wheel_draw({k: p[k] for k in WHEEL_KEYS})


# ── [Z] the ζ side ────────────────────────────────────────────────────────────────────────────────
ZT0 = (PP.GAMMA1 - 3.0) / 27.0                      # τ at which the ζ tile sits on the first zero


def z_height(p, tau):
    return 3.0 + 27.0 * tau


def z_euler(p, tau):
    t = z_height(p, tau); L = PP.log_product_path(t, PP.first_primes(int(p["K"])), .5)
    return _unwrap_polar(L - L[0] * 0)


def _dirichlet_path(t, N, weight):
    n = np.arange(2, N + 1, dtype=float); lam = TZ.vonmangoldt_table(N)[2:]
    return np.cumsum(weight(lam) * n ** (-(.5 + 1j * t)))


def z_telperion(p, tau):
    return _unwrap_polar(_dirichlet_path(z_height(p, tau), int(p["N"]), lambda lam: lam))


def z_laurelin(p, tau):
    return _unwrap_polar(_dirichlet_path(z_height(p, tau), int(p["N"]), lambda lam: 1 - lam))


# ── registry ──────────────────────────────────────────────────────────────────────────────────────
@dataclass
class Pattern:
    id: str; label: str; source: str; fn: object; params: dict = field(default_factory=dict)
    tau0: float = .6; kind: str = "polar"; note: str = ""
    state: object = None          # wheel-family patterns: tau,p → state dict


PI = dict(pi_f=PI_INV, e_f=E_INV)
PATTERNS = {q.id: q for q in [
    Pattern("O1", "Orientation", "[O] original ph.1", o_orientation, dict(PI), .4, note="inward spiral, θ scaled by the 2/π slider"),
    Pattern("O2", "Worm (metric swap)", "[O] original ph.2", o_worm, dict(PI), .6, "wheel", "phase modulation β = 5", worm_state),
    Pattern("O3", "Rosette (Fermat frequency)", "[O] original ph.3", o_rosette, dict(PI), np.pi / 4, "wheel", "amplitude modulation 0.2", rosette_state),
    Pattern("O4", "Prime spiral", "[O] original ph.4", o_primes, dict(PI), 1.0, note="θ = p·(2/π), r = √p·(2/e)"),
    Pattern("C1", "Inertial", "[C] Claude ph.1", c_inertial, {}, .5), Pattern("C2", "Metric swap (inductance)", "[C] Claude ph.2", c_metric, {}, .5),
    Pattern("C3", "Heartbeat (co-directed)", "[C] Claude ph.3", c_heartbeat, {}, .5), Pattern("C4", "Golden prime spiral", "[C] Claude ph.4", c_golden, {}, 1.0),
    Pattern("C5", "Horizon invariance (Fermat lattice)", "[C] Claude ph.5", c_horizon, {}, 1.0),
    Pattern("G1", "φ-flow  z ↦ 1 + 1/z", "[G] oldest", g_phiflow, {}, 1.0, note="Fibonacci ratios → φ"),
    Pattern("R1", "Inertial", "[R] Revised_2", r_phase("Inertial"), {}, .1), Pattern("R2", "Metric Swap", "[R] Revised_2", r_phase("Metric Swap"), {}, .1),
    Pattern("R3", "Heartbeat", "[R] Revised_2", r_phase("Heartbeat"), {}, .1), Pattern("R4", "Star (phyllotaxis)", "[R] Revised_2", r_phase("Star"), {}, .1),
    Pattern("X1", "Crossing r = 1 and inversion", "[X] Transition", x_cross, {}, .4, note="inverted branch r → (1+0.0007)/r"),
    Pattern("W1", "Wheel family (state lerp)", "this repo", w_family, dict(Th=100 * np.pi, beta=2.5, om=.5, d=.1, k=.3, m=.5, a=.3), 1.0, "wheel",
            "worm ⇄ rosette continuation"),
    Pattern("Z1", "ζ Euler path  log P_K", "[Z] this repo", z_euler, dict(K=2000), ZT0, note="s = ½ + it, t = 3 + 27τ"),
    Pattern("Z2", "Telperion  Σ Λ(n) n⁻ˢ", "[Z] this repo", z_telperion, dict(N=20000), ZT0, note="T = −ζ′/ζ partial sums"),
    Pattern("Z3", "Laurelin  Σ (1−Λ) n⁻ˢ", "[Z] this repo", z_laurelin, dict(N=20000), ZT0, note="L = ζ − 1 + ζ′/ζ partial sums"),
]}

# the grid the viewer scrolls over: one row per source file
GRID = [["O1", "O2", "O3", "O4", None], ["C1", "C2", "C3", "C4", "C5"], ["R1", "R2", "R3", "R4", None],
        ["G1", "X1", "W1", None, None], ["Z1", "Z2", "Z3", None, None]]
ROW_NAMES = ["[O] original_2026/event_horizon_sim.py", "[C] Gemini/event_horizon_sim.py (Claude)", "[R] Gemini_Revised_2_EHS.py",
             "[G][X] Gemini_event_horizon_sim.py · Gemini_Transition_EHS.py · wheel family", "[Z] the ζ pathway — whole, Telperion, Laurelin"]


def points(pid, tau=None, params=None):
    """Complex points of pattern pid at clock τ (default: the tile's own τ₀)."""
    q = PATTERNS[pid]; p = {**q.params, **(params or {})}
    th, r = q.fn(p, q.tau0 if tau is None else tau)
    return np.asarray(r) * np.exp(1j * np.asarray(th))


def polar(pid, tau=None, params=None):
    q = PATTERNS[pid]; p = {**q.params, **(params or {})}
    th, r = q.fn(p, q.tau0 if tau is None else tau)
    return np.asarray(th, float), np.asarray(r, float)


# ── transition mappings ───────────────────────────────────────────────────────────────────────────
TRANSITION_MODES = ("polar", "cartesian", "wheel")


def _resample(a, n):
    return np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(a)), a)


def transition(pa, pb, alpha, mode="polar", n=None, tau_a=None, tau_b=None):
    """Points of the A → B transition at α ∈ [0, 1].
    polar      Gemini's TransitionManager: index-matched LERP of θ and of r.
    cartesian  LERP of the complex points.
    wheel      LERP of the WHEEL-FAMILY state (Th, β, ω, d, k, m, a) — only when both patterns carry a state (O2, O3, W1); else falls back to polar.
    Returns (z, mode_used)."""
    A, B = PATTERNS[pa], PATTERNS[pb]
    n = n or (5000 if mode == "wheel" else 3000)            # 5000 = the originals' own point count, so wheel endpoints equal them exactly
    ta = A.tau0 if tau_a is None else tau_a; tb = B.tau0 if tau_b is None else tau_b
    if mode == "wheel" and A.state and B.state:
        sa, sb = A.state({**A.params}, ta) if A.id != "W1" else {k: A.params[k] for k in WHEEL_KEYS}, \
                 B.state({**B.params}, tb) if B.id != "W1" else {k: B.params[k] for k in WHEEL_KEYS}
        s = {k: (1 - alpha) * sa[k] + alpha * sb[k] for k in WHEEL_KEYS}
        ang, r = wheel_draw(s, n); return r * np.exp(1j * ang), "wheel"
    tha, ra = polar(pa, ta); thb, rb = polar(pb, tb)
    if mode == "cartesian":
        za = ra * np.exp(1j * tha); zb = rb * np.exp(1j * thb)
        return (1 - alpha) * _resample(za.real, n) + alpha * _resample(zb.real, n) + 1j * ((1 - alpha) * _resample(za.imag, n) + alpha * _resample(zb.imag, n)), "cartesian"
    th = (1 - alpha) * _resample(tha, n) + alpha * _resample(thb, n); r = (1 - alpha) * _resample(ra, n) + alpha * _resample(rb, n)
    return r * np.exp(1j * th), "polar" if mode != "wheel" else "polar(fallback)"


# ── turn-around diagnostics ───────────────────────────────────────────────────────────────────────
def turn_stats(z, tol=1e-12):
    """What 'turning around' means, measured two ways on a polyline z:
    angular_reversals  sign changes of the angular momentum about the origin Im(conj(z)·dz)  — the curve runs BACKWARDS in angle (retrograde)
    curvature_reversals sign changes of the turning of the tangent Im(conj(dz_k)·dz_{k+1})  — the tangent starts rotating the other way
    cusps              consecutive steps that fold back by more than 150°
    winding            net turns of the tangent;  length  arclength."""
    z = np.asarray(z, complex); dz = np.diff(z)
    L = np.imag(np.conj(z[:-1]) * dz); K = np.imag(np.conj(dz[:-1]) * dz[1:])
    def flips(x):
        x = x[np.abs(x) > tol * (np.abs(x).max() + 1e-300)]; s = np.sign(x); return int(np.sum(s[1:] != s[:-1])) if len(s) > 1 else 0
    ph = np.unwrap(np.angle(dz)); d = np.diff(ph)
    return dict(angular_reversals=flips(L), curvature_reversals=flips(K), cusps=int(np.sum(np.abs(np.angle(dz[1:] * np.conj(dz[:-1]))) > np.deg2rad(150))),
                winding=float((ph[-1] - ph[0]) / (2 * np.pi)), length=float(np.sum(np.abs(dz))), steps=len(dz))


def transition_curve(pa, pb, mode="polar", m=41, n=3000):
    """turn_stats along α ∈ [0, 1] — where the transition makes or loses turn-arounds."""
    al = np.linspace(0, 1, m); rows = [turn_stats(transition(pa, pb, a, mode, n)[0]) for a in al]
    return al, {k: np.array([r[k] for r in rows]) for k in rows[0]}


def retro_fraction(beta, omega):
    """Fraction of θ_b the worm runs backwards in angle: dθ/dθ_b = 1 + βω cos(ωθ_b) < 0.  Zero if βω ≤ 1, else arccos(1/(βω))/π.  EXACT."""
    a = beta * omega
    return 0.0 if a <= 1 else float(np.arccos(1 / a) / np.pi)


def worm_retro_threshold(beta=5.0):
    """The cycle time below which the worm never turns back: t < 1/β (= 0.2 for the original's β = 5)."""
    return 1 / beta


# ── how the ζ pathway turns ───────────────────────────────────────────────────────────────────────
def zeta_turning(t=PP.GAMMA1, K=20000):
    """The ζ Euler path z_K = log P_K(½ + it).  Its k-th step is −ln(1 − p^{-s}) = p^{-s}(1 + p^{-s}/2 + …), direction −t ln p − (small).
    Consecutive steps turn by −t·ln(p_{k+1}/p_k) (mod 2π) — ONE sense for every k once t·ln(p_{k+1}/p_k) < π, i.e. beyond the first few primes.
    Returns the diagnostics and the exact carrier law."""
    P = PP.first_primes(K); z = PP.log_product_path(t, P, .5); dz = np.diff(np.concatenate([[0], z]))
    phi = np.unwrap(np.angle(dz)); carrier = -t * np.log(P)
    resid = phi - carrier; resid = resid - np.round(resid[len(resid) // 2] / (2 * np.pi)) * 2 * np.pi
    turn = np.diff(phi); sense = np.sign(np.sin(-t * np.diff(np.log(P))))
    st = turn_stats(z)
    first_flip = int(np.max(np.nonzero(sense[1:] != sense[:-1])[0]) + 2) if np.any(sense[1:] != sense[:-1]) else 0
    return dict(stats=st, P=P, last_sense_flip_prime_index=first_flip, last_sense_flip_prime=int(P[min(first_flip, len(P) - 1)]),
                tangent_resid_rms_by_decade={f"p<{10 ** d}": float(np.sqrt(np.mean(resid[(P >= 10 ** (d - 1)) & (P < 10 ** d)] ** 2))) for d in (2, 3, 4)
                                              if np.any((P >= 10 ** (d - 1)) & (P < 10 ** d))})


def tone_share(x, nfft=None):
    """Fraction of the variance of a signal in its single strongest non-DC Fourier line — 1 for a pure tone."""
    x = np.asarray(x, float); x = x - x.mean()
    f = np.abs(np.fft.rfft(x * np.hanning(len(x)), nfft)) ** 2
    return float(f[1:].max() / f[1:].sum())


def turning_rate_tone_share(z, resample=4096):
    """Tone share of the tangent-angle RATE (after removing the carrier trend) of a curve, taken against arclength.  Worm: one dominant tone."""
    z = np.asarray(z, complex); dz = np.diff(z); s = np.concatenate([[0], np.cumsum(np.abs(dz))])
    ph = np.unwrap(np.angle(dz)); sm = .5 * (s[1:] + s[:-1])
    sg = np.linspace(sm[0], sm[-1], resample); phg = np.interp(sg, sm, ph)
    rate = np.gradient(phg, sg); rate = rate - np.polyval(np.polyfit(sg, rate, 3), sg)
    return tone_share(rate)


def zeta_vs_worm(ts=(3.0, PP.GAMMA1, 30.0, 100.0), K=20000):
    """The comparison table: how the ζ Euler path and the worm turn.  Rows: ζ at several heights t; the worm at several τ."""
    rows = []
    for t in ts:
        z = PP.log_product_path(t, PP.first_primes(K), .5); s = turn_stats(z); zt = zeta_turning(t, K)
        rows.append(dict(kind="zeta", at=f"t={t:g}", angular=s["angular_reversals"], curvature=s["curvature_reversals"], cusps=s["cusps"], winding=s["winding"],
                         last_sense_flip_prime=zt["last_sense_flip_prime"], tone=turning_rate_tone_share(z), resid=zt["tangent_resid_rms_by_decade"]))
    for tau in (.15, .6, 1.0):
        z = points("O2", tau); s = turn_stats(z)
        rows.append(dict(kind="worm", at=f"τ={tau:g}", angular=s["angular_reversals"], curvature=s["curvature_reversals"], cusps=s["cusps"], winding=s["winding"],
                         last_sense_flip_prime=None, tone=turning_rate_tone_share(z), resid=None, retro=retro_fraction(5, tau)))
    return rows


def zeta_wobble(t=PP.GAMMA1, K=9592, lo=1e4, hi=1e5):
    """rms of the ζ step-angle residual (step angle + t ln p) for p ∈ [lo, hi) against the analytic ½·p^{-½}/√2 (the second term of −ln(1 − p^{-s}) = p^{-s}(1 + p^{-s}/2 + …)).
    The ζ 'phase-modulation index' is β_ζ(p) = ½ p^{-½} — it decays; the worm's is the constant 5."""
    P = PP.first_primes(K); z = PP.log_product_path(t, P, .5); dz = np.diff(np.concatenate([[0], z]))
    resid = np.unwrap(np.angle(dz)) + t * np.log(P); m = (P >= lo) & (P < hi)
    resid = resid - np.round(np.median(resid[m]) / (2 * np.pi)) * 2 * np.pi
    return float(np.sqrt(np.mean((resid[m] - resid[m].mean()) ** 2))), float(np.mean(.5 / np.sqrt(P[m])) / np.sqrt(2))
