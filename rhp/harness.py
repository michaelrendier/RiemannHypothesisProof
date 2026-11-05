# SPDX-License-Identifier: GPL-3.0-only
"""Operator test harness: score ANY candidate spectrum against the zeros.

A candidate is a sorted array of 'energies' E_k.  Four scores, against the zeros γ_k (data/zeros_1_to_3600.csv):
  weyl      max_k |θ(E_k)/π + 1 − (k − ½)|        (0 for the carrier; = max|S| for the zeros)
  spacing   variance of unfolded nearest-neighbour spacings   (zeros / GUE ≈ 0.18–0.27, Poisson 1, a ladder 0)
  numvar    number variance Σ²(L) of the unfolded levels in windows of L spacings  (L = 5)
  locate    mean |E_k − γ_k| in mean spacings        (0 = the zeros; ~0.2 = the carrier; large = unrelated)
Controls: the carrier (θ = (k−3/2)π), a ladder (n+½)ħ (the ValaQuenta H_NN spectrum and the xp-in-a-periodic-u-box spectrum),
a Poisson spectrum, and a GUE spectrum unfolded onto the zeros' density.
"""
from __future__ import annotations
import numpy as np
from . import bigpath, clock


def unfold(E):
    """Smooth counting function at each level: u_k = θ(E_k)/π + 1."""
    return bigpath.theta_fast(np.asarray(E, float)) / np.pi + 1


def numvar(u, L=5.0, nwin=4000, seed=0):
    rng = np.random.default_rng(seed)
    starts = rng.uniform(u[0], u[-1] - L, nwin)
    counts = np.searchsorted(u, starts + L) - np.searchsorted(u, starts)
    return float(counts.var())


def score(E, zeros):
    E = np.sort(np.asarray(E, float))
    u = unfold(E)
    sp = np.diff(u)
    k = np.arange(1, len(E) + 1)
    weyl = float(np.max(np.abs(u - (k - 0.5))))               # u = θ/π + 1; the carrier has u = k − ½ exactly
    n = min(len(E), len(zeros))
    mean_sp = np.mean(np.diff(zeros[:n]))
    loc = float(np.mean(np.abs(E[:n] - zeros[:n])) / mean_sp) if n else np.nan
    return dict(weyl=weyl, spacing_var=float(sp.var()), numvar5=numvar(u), locate=loc)


def gue_spectrum(N, seed=0):
    """GUE eigenvalues (size N) unfolded by the semicircle law onto the zeros' smooth density, then mapped back to energies."""
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    H = (A + A.conj().T) / (2 * np.sqrt(2 * N))              # spectrum ≈ semicircle on [−√2·… ] -> unfold by CDF
    ev = np.linalg.eigvalsh(H)
    R = np.sqrt(2.0) * 1.0                                     # semicircle radius for this normalisation
    x = np.clip(ev / R, -1, 1)
    cdf = 0.5 + (x * np.sqrt(1 - x * x) + np.arcsin(x)) / np.pi
    u = cdf * N                                                # unfolded levels, mean spacing 1
    return u


def energies_from_u(u):
    """Invert u = θ(E)/π + 1 − 3/2 + offset: place unfolded levels on the zeros' smooth density (carrier map)."""
    from .clock import carrier_zeros
    base = carrier_zeros(int(np.ceil(u.max())) + 2)
    return np.interp(u, np.arange(1, len(base) + 1) - 1.0 + 0.5, base)


def poisson_u(N, seed=1):
    rng = np.random.default_rng(seed)
    return np.sort(rng.uniform(0, N, N))


def candidates(zeros, N=None):
    N = N or len(zeros)
    out = {"zeros (reference)": zeros[:N]}
    out["carrier  θ(E)=(k−3/2)π"] = clock.carrier_zeros(N)
    hbar = (zeros[N - 1] - zeros[0]) / (N - 1)
    out["ladder (n+½)ħ  [H_NN, xp in a periodic u-box]"] = zeros[0] + hbar * (np.arange(N) + 0.5)
    out["Poisson"] = energies_from_u(poisson_u(N))
    out["GUE (unfolded)"] = energies_from_u(gue_spectrum(N)[:N])
    return out
