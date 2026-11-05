# SPDX-License-Identifier: GPL-3.0-only
"""Bifurcation margins: how close the local extrema of Z(t) between zeros come to touching zero (Lehmer's phenomenon).

A double zero of the real function Z splits under a small real perturbation into two real zeros (stays on the line)
or a symmetric pair (leaves it).  The margin of the gap between γ_n and γ_{n+1} is the height of Z's extremum there,
relative to the typical extremum height around it.
"""
from __future__ import annotations
import numpy as np
from mpmath import mp, mpf, siegelz
from . import bigpath


def gap_margins(p, zeros):
    t, Z = p["t"], np.abs(p["Z"])
    idx = np.searchsorted(t, zeros)
    peaks = np.array([Z[idx[i]:idx[i + 1]].max() if idx[i + 1] > idx[i] else np.nan for i in range(len(zeros) - 1)])
    at = np.array([t[idx[i] + np.argmax(Z[idx[i]:idx[i + 1]])] if idx[i + 1] > idx[i] else np.nan for i in range(len(zeros) - 1)])
    loc = np.array([np.nanmedian(peaks[max(0, i - 50):i + 51]) for i in range(len(peaks))])
    rel = peaks / loc
    gaps = np.diff(zeros)
    sp = np.pi / np.maximum(bigpath.theta_fast(zeros[:-1] + gaps / 2) * 0 + (0.5 * np.log(zeros[:-1] / (2 * np.pi))), 1e-9)
    return dict(gamma=zeros[:-1], gap=gaps, gap_norm=gaps / sp, peak=peaks, peak_t=at, rel=rel)


def refine_peak(a, b):
    """Maximise |Z| on [a, b] (the gap between two consecutive zeros) with mpmath (golden section); returns (t*, |Z(t*)|)."""
    mp.dps = 15
    gr = (np.sqrt(5) - 1) / 2
    c, d = b - gr * (b - a), a + gr * (b - a)
    f = lambda x: -abs(float(siegelz(mpf(x))))
    fc, fd = f(c), f(d)
    for _ in range(60):
        if fc < fd: b, d, fd = d, c, fc; c = b - gr * (b - a); fc = f(c)
        else: a, c, fc = c, d, fd; d = a + gr * (b - a); fd = f(d)
    x = (a + b) / 2
    return x, -f(x)
