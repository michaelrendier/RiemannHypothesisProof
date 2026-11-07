# SPDX-License-Identifier: GPL-3.0-only
"""The divergence detector: does any mode of the prime-counting residual grow (a DILATE) where only rotations are allowed?

Setting (ESTABLISHED, von Mangoldt):  ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²).   With u = ln x and E(u) = (x − ψ(x) − ln 2π − …)/√x,
  E(u) = Σ_ρ x^{σ−½} e^{iγu}/ρ     (each zero ρ = σ + iγ; the sum over ρ and ρ̄ is real).
On the critical line (σ = ½) every term is a pure rotation in u (gain 1): isometries only.  An off-line zero ρ₀ = ½+δ+iγ₀ brings its whole
quadruple {ρ₀, ρ̄₀, 1−ρ₀, 1−ρ̄₀}: at frequency γ₀ the amplitude is a·e^{δu} + b·e^{−δu} — a DILATE (length-changing) component, graded by δ.
Lineage (generational-lineage §5, first signature): 'a quantity that changes length where only isometries were in play'.

Detector: a Hann-windowed lock-in at frequency γ over consecutive windows of u gives amplitude A_k(γ); the statistic is the slope of ln A against u.
  CONTROL first (§8): on the real ψ the slope at every true zero is ≈ 0 — the measured scatter is the false-alarm floor.
  Then INJECT a known quadruple at δ and recover its slope: the smallest δ the detector can see at this reach u_max.
Nothing here is evidence about RH: a null result on the real ψ over u ≤ u_max says only that no growth ≥ the floor is present there.
"""
from __future__ import annotations

import math
import pathlib
import numpy as np

DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "zeros_1_to_3600.csv"


def load_zeros(k=None):
    z = np.loadtxt(DATA)
    return z if k is None else z[:k]


def psi_table(X):
    """ψ(n) for n = 0..X (Chebyshev, exact): cumulative Σ ln p over prime powers."""
    s = np.ones(X + 1, bool); s[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    lam = np.zeros(X + 1)
    for p in np.flatnonzero(s):
        q = int(p); lp = math.log(p)
        while q <= X:
            lam[q] = lp
            q *= int(p)
    return np.cumsum(lam)


def residual(u, psi):
    """E(u) = (x − ψ(x) − ln 2π − ½ln(1−x⁻²))/√x at x = e^u (ψ read at floor(x), so the sample sits between jumps)."""
    x = np.exp(u)
    ps = psi[np.minimum(np.floor(x).astype(np.int64), len(psi) - 1)]
    return (x - ps - math.log(2 * math.pi) - 0.5 * np.log1p(-x ** -2.0)) / np.sqrt(x)


def quadruple(u, delta, gamma0, phase=0.0, scale=1.0):
    """The contribution to E(u) of an off-line zero at ½+δ+iγ₀ with its three partners (conjugate, mirror 1−ρ, mirror-conjugate).
    δ = 0 reduces to the ordinary pair on the line (counted twice: ρ and 1−ρ̄ coincide), so use δ > 0."""
    out = np.zeros_like(u)
    for sg in (+delta, -delta):                     # σ − ½ = ±δ
        rho = complex(0.5 + sg, gamma0)
        out += 2 * scale * np.real(np.exp(sg * u) * np.exp(1j * (gamma0 * u + phase)) / rho)
    return out


def lockin(u, E, gamma, ua, ub):
    """Hann-windowed amplitude of E at frequency γ on [ua, ub]:  |2∫ E w e^{−iγu}du / ∫w|."""
    m = (u >= ua) & (u <= ub)
    uu, ee = u[m], E[m]
    w = np.hanning(len(uu))
    return abs(2 * np.sum(ee * w * np.exp(-1j * gamma * uu)) / np.sum(w))


def amplitude_track(u, E, gamma, u0, u1, width, step):
    """(centres, amplitudes) of the lock-in over sliding windows of the given width in u."""
    cs, As = [], []
    a = u0
    while a + width <= u1 + 1e-9:
        cs.append(a + width / 2); As.append(lockin(u, E, gamma, a, a + width)); a += step
    return np.array(cs), np.array(As)


def growth_slope(u, E, gamma, u0=6.0, u1=16.0, width=3.0, step=1.0):
    """d ln A / du at frequency γ (least squares over the windows): 0 for an on-line tone, ≈ δ for a lone off-line one."""
    c, A = amplitude_track(u, E, gamma, u0, u1, width, step)
    return float(np.polyfit(c, np.log(A), 1)[0])


def control(psi, K=20, u0=6.0, u1=16.0, n=60000):
    """Slopes at the first K true zeros on the real ψ: returns (slopes, mean, sd).  sd is the false-alarm floor."""
    u = np.linspace(u0, u1, n)
    E = residual(u, psi)
    g = load_zeros(K)
    s = np.array([growth_slope(u, E, x, u0, u1) for x in g])
    return s, float(s.mean()), float(s.std(ddof=1))


def inject_and_recover(psi, deltas, gamma0=14.134725, u0=6.0, u1=16.0, n=60000, scale=1.0):
    """Add a known quadruple at frequency γ₀ to the real E and measure the slope there: (δ, slope) pairs."""
    u = np.linspace(u0, u1, n)
    E = residual(u, psi)
    return [(d, growth_slope(u, E + quadruple(u, d, gamma0, scale=scale), gamma0, u0, u1)) for d in deltas]


# ── what the detector can see, and where ──────────────────────────────────────────────────────
def detection_limit(psi, u1, k_sigma=3.0, K=20, u0=6.0, n=60000, gamma0=14.134725):
    """Smallest δ whose injected quadruple (at γ₀) gives a slope ≥ k_sigma × (control sd on [u0,u1]).  Bisection on δ.
    Because the mirror partner cancels the first-order growth, the slope of ln(a·e^{δu}+b·e^{−δu}) is O(δ²·u) for δu ≪ 1,
    so the limit falls like 1/u_max, not like the noise."""
    _, _, sd = control(psi, K, u0, u1, n)
    u = np.linspace(u0, u1, n)
    E = residual(u, psi)
    lo, hi = 1e-4, 1.0
    for _ in range(18):
        mid = math.sqrt(lo * hi)
        s = growth_slope(u, E + quadruple(u, mid, gamma0), gamma0, u0, u1)
        lo, hi = (mid, hi) if s < k_sigma * sd else (lo, mid)
    return math.sqrt(lo * hi), sd


def scan(psi, gammas, u0=6.0, u1=16.0, n=60000, extra=None):
    """Slope and mean amplitude at each frequency γ of E (+ an optional injected term `extra(u)`): the candidate map.
    Returns (gammas, slopes, mean_amplitude)."""
    u = np.linspace(u0, u1, n)
    E = residual(u, psi)
    if extra is not None:
        E = E + extra(u)
    sl, am = [], []
    for g in gammas:
        c, A = amplitude_track(u, E, g, u0, u1, 3.0, 1.0)
        sl.append(np.polyfit(c, np.log(A), 1)[0]); am.append(A.mean())
    return np.asarray(gammas), np.array(sl), np.array(am)


def candidates(gammas, slopes, amps, slope_floor, amp_floor):
    """Frequencies whose slope exceeds the control floor AND whose amplitude is a real tone (above amp_floor)."""
    return [float(g) for g, s, a in zip(gammas, slopes, amps) if s > slope_floor and a > amp_floor]


def candidate_peaks(gammas, slopes, amps, slope_floor, amp_floor):
    """Group contiguous candidate frequencies and report each group's strongest frequency: [(γ_peak, amplitude, slope)]."""
    idx = [i for i in range(len(gammas)) if slopes[i] > slope_floor and amps[i] > amp_floor]
    groups, cur = [], []
    for i in idx:
        if cur and i != cur[-1] + 1:
            groups.append(cur); cur = []
        cur.append(i)
    if cur:
        groups.append(cur)
    out = []
    for gr in groups:
        j = max(gr, key=lambda i: amps[i])
        out.append((float(gammas[j]), float(amps[j]), float(slopes[j])))
    return out
