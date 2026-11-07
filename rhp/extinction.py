# SPDX-License-Identifier: GPL-3.0-only
"""rhp.extinction — the extinction flow: the two books (Telperion / Laurelin) as CONTINUOUS flows in τ = ln y, and their Laplace dual in s.

Working axis (the Archimedes Screw's own):  u = ln x.  The cutoff ('time' of the flow) is  τ = ln y.  The ratio  U = u/τ = ln x / ln y.

ESTABLISHED (textbook; de Bruijn 1951, Buchstab 1937, Tenenbaum ch. III.5-6):
  Telperion (sieve) survivors   Φ(x,y) = #{n ≤ x : spf(n) > y}      Φ(x,y) ~ x·ω(U)/ln y      (U·ω)′ = ω(U−1),  ω = 1/U on [1,2]
  Laurelin  (un-sieve) smooth   Ψ(x,y) = #{n ≤ x : gpf(n) ≤ y}      Ψ(x,y) ~ x·ρ(U)           U·ρ′ = −ρ(U−1),           ρ = 1 on [0,1]
  Both are delay-differential equations with delay exactly 1 in U; ω → e^{−γ}, ρ → 0 superexponentially.
OURS-CALCULATED (derived here from the two equations above by Laplace transform in U; checked numerically in tests; I have
  not found it stated in this form, so it is labelled OURS and may turn out to be known):
      ρ̂(s) = e^{−E₁(s)}/s ,     1 + ω̂(s) = e^{+E₁(s)}      ⇒      s·ρ̂(s) · (1 + ω̂(s)) = 1
  the two books are RECIPROCAL under the Laplace transform in U — the same relation as ζ(s)·(1/ζ(s)) = 1 in the Dirichlet domain
  (1/ζ = Σ μ(n)n^{−s} is the sieve; ζ is the un-sieve).  E₁ is the exponential integral; its log branch at s = 0 sits where
  the normalised variable s′ = s − 1 meets ζ's pole (the north pole of the equatorial sphere).
"""
from __future__ import annotations

import math
import numpy as np
from mpmath import mp, e1 as _e1, mpc, exp as _exp

GAMMA = 0.5772156649015329


# ── the two delay-differential equations, solved on a fine grid ───────────────────────────────
def buchstab(umax=14.0, h=1.0 / 2000):
    """ω(U) on [0, umax]:  0 for U<1,  1/U on [1,2],  U·ω(U) = 1 + ∫₁^{U−1} ω  beyond.  Returns (grid, ω)."""
    n = int(round(umax / h))
    u = np.arange(n + 1) * h
    w = np.zeros(n + 1)
    k1 = int(round(1 / h))
    w[k1:] = 1.0 / u[k1:]                           # first guess on [1, 2]; recomputed beyond by recurrence
    cum = np.zeros(n + 1)                           # cum[k] = ∫₀^{u_k} ω
    for k in range(1, n + 1):
        if u[k] > 2 + 1e-12:
            j = k - k1                              # index of U−1
            w[k] = (1.0 + cum[j] - cum[k1]) / u[k]  # ∫₁^{U−1} ω  = cum[j] − cum[k1]
        elif u[k] >= 1 - 1e-12:
            w[k] = 1.0 / u[k]
        cum[k] = cum[k - 1] + 0.5 * h * (w[k] + w[k - 1])
    return u, w


def dickman(umax=14.0, h=1.0 / 2000):
    """ρ(U) on [0, umax]:  1 on [0,1],  U·ρ′(U) = −ρ(U−1).  Trapezoid march.  Returns (grid, ρ)."""
    n = int(round(umax / h))
    u = np.arange(n + 1) * h
    r = np.ones(n + 1)
    k1 = int(round(1 / h))
    for k in range(k1, n):
        a = r[k - k1] / u[k] if u[k] > 0 else 0.0
        # implicit trapezoid in ρ(U−1): at k+1 the delayed value r[k+1−k1] is already known
        b = r[k + 1 - k1] / u[k + 1]
        r[k + 1] = r[k] - 0.5 * h * (a + b)
    return u, r


# ── the flow measured on the integers ─────────────────────────────────────────────────────────
def survivors(spf, N, y):
    """Telperion: #{2 ≤ n ≤ N : spf(n) > y}, plus 1 for n = 1 (nothing has touched it)."""
    return int(np.sum(spf[2:N + 1] > y)) + 1


def smooth(gpf, N, y):
    """Laurelin: #{2 ≤ n ≤ N : gpf(n) ≤ y}, plus 1 for n = 1 (the empty product)."""
    return int(np.sum(gpf[2:N + 1] <= y)) + 1


def mertens(y):
    """Π_{p ≤ y}(1 − 1/p) — the sieve's survivor density; ≈ e^{−γ}/ln y (Mertens 1874, established)."""
    n = int(y)
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    p = np.flatnonzero(s).astype(float)
    return float(np.prod(1 - 1 / p))


# ── the Laplace dual ──────────────────────────────────────────────────────────────────────────
def laplace_numeric(u, f, s, tail=None):
    """∫₀^umax f(u) e^{−su} du by trapezoid (+ an analytic tail  c·e^{−s·umax}/s  if tail=c is given)."""
    g = f * np.exp(-s * u)
    h = u[1] - u[0]
    val = h * (g.sum() - 0.5 * (g[0] + g[-1]))
    if tail is not None:
        val += tail * np.exp(-s * u[-1]) / s
    return val


def rho_hat_closed(s):
    """ρ̂(s) = e^{−E₁(s)}/s."""
    with mp.workdps(30):
        s = mpc(s)
        return complex(_exp(-_e1(s)) / s)


def omega_hat_closed(s):
    """ω̂(s) = e^{+E₁(s)} − 1."""
    with mp.workdps(30):
        return complex(_exp(_e1(mpc(s))) - 1)


def reciprocal_residual(s, u_rho, rho, u_om, om):
    """|s·ρ̂(s)·(1 + ω̂(s)) − 1| computed from the NUMERIC ρ and ω (not the closed forms): tests the reciprocity itself."""
    r_hat = laplace_numeric(u_rho, rho, s)
    w_hat = laplace_numeric(u_om, om, s, tail=math.exp(-GAMMA))
    return abs(s * r_hat * (1 + w_hat) - 1)


def log_gain(s):
    """Re E₁(s′): ln|1 + ω̂(s′)| — the gain of the sieve book at the complex frequency s′; −ln|s′ρ̂(s′)| is the same number."""
    with mp.workdps(20):
        return float(_e1(mpc(s)).real)


# ── order of approach: the two paths by which n is reached ────────────────────────────────────
def factor_sequence(n):
    """Prime factors of n with multiplicity, ascending (the order Telperion strips them: smallest first)."""
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d); n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def approach_paths(n):
    """(Telperion path, Laurelin path): the same multiset of primes taken in the two opposite orders — ascending = sieve strips the
    smallest first; descending = the wheels complete n by switching on the largest last.  They coincide iff n is a prime power."""
    a = factor_sequence(n)
    return a, a[::-1]


def span(n):
    """ln(gpf/spf): how many e-folds of scale n's arrival crosses between dying in the sieve (τ = ln spf) and being born in the
    wheels (τ = ln gpf).  Zero exactly for prime powers (both clocks fire at once)."""
    a = factor_sequence(n)
    return math.log(a[-1] / a[0])
