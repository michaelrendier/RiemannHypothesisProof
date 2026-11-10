# SPDX-License-Identifier: GPL-3.0-only
"""The nodal-line extractor: the zero set of a standing-wave field, its crossings, its domains.   (cymatics: where the sand collects)

A Chladni figure is the zero set of a field built from standing waves.  Here the fields are
  1-D  the prime wave  W(u) = (ψ(x) − x)/√x,  u = ln x, built from the zero-tones  (heartbeat)   → nodes = the 'prime crossings' ψ(x) = x
  2-D  the real and imaginary parts of ξ(σ + it) on the strip                                    → nodal lines Re ξ = 0, Im ξ = 0; their intersections are the zeros
Everything is extracted by the code below (no plotting library): marching squares with linear interpolation and the saddle rule, segments stitched
into polylines (closed / open, length), exact segment–segment intersections, sign-domain labelling (Courant's nodal domains), 1-D nodes with direction,
and Rice's law for the node rate of a band-limited sum of standing tones.  Labels: ESTABLISHED = textbook; :CALCULATED = floating point from this file.
"""
from __future__ import annotations
import math
import numpy as np


# ------------------------------------------------------------------ 1-D
def sign_nodes(x, f):
    """Nodes of a sampled function: linear interpolation of every sign change.  Returns (positions, direction) with direction +1 (upward) / −1 (downward)."""
    x = np.asarray(x, float); f = np.asarray(f, float)
    s = np.sign(f); idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    pos = x[idx] - f[idx] * (x[idx + 1] - x[idx]) / (f[idx + 1] - f[idx])
    return pos, np.sign(f[idx + 1] - f[idx]).astype(int)


def refine_nodes(fun, x, f, tol=1e-12):
    """Bisection (regula falsi) of each bracket of a sampled sign change against the callable `fun`, to tol."""
    s = np.sign(f); out = []
    for i in np.nonzero(s[:-1] * s[1:] < 0)[0]:
        a, b = float(x[i]), float(x[i + 1]); fa, fb = float(f[i]), float(f[i + 1])
        for _ in range(100):
            c = b - fb * (b - a) / (fb - fa); fc = fun(c)
            if fc == 0 or abs(b - a) < tol: break
            if fa * fc < 0: b, fb = c, fc
            else: a, fa = c, fc
        out.append(c)
    return np.array(out)


def psi_table(X):
    """ψ(n) for n = 0 … X (Chebyshev, exact cumulative sum of Λ)."""
    lam = np.zeros(X + 1); sv = np.ones(X + 1, bool); sv[:2] = False
    for p in range(2, int(X ** .5) + 1):
        if sv[p]: sv[p * p::p] = False
    for p in np.nonzero(sv)[0]:
        q = int(p)
        while q <= X: lam[q] = math.log(p); q *= int(p)
    return np.cumsum(lam), lam


def exact_prime_nodes(x0, x1, psi=None):
    """Exact nodes of ψ(x) − x on [x0, x1] (integers x0 ≤ n ≤ x1).  ψ is constant between integers and jumps up by Λ(n) AT n, so
    – on (n, n+1): ψ(n) − x is decreasing; it has a node at x = ψ(n) iff n < ψ(n) < n + 1                    (direction −1)
    – at x = n: the jump crosses zero iff ψ(n−1) − n and ψ(n) − n have opposite signs (direction ±1)           (the node sits AT the integer)
    Returns (positions, direction, kind) with kind 'jump' / 'slide'."""
    X = int(x1) + 2
    if psi is None: psi, _ = psi_table(X)
    pos, dr, kd = [], [], []
    for n in range(max(2, int(math.ceil(x0))), int(x1) + 1):
        a, b = psi[n - 1] - n, psi[n] - n
        if a * b < 0: pos.append(float(n)); dr.append(1 if b > a else -1); kd.append("jump")
        if n < psi[n] < n + 1: pos.append(float(psi[n])); dr.append(-1); kd.append("slide")
    return np.array(pos), np.array(dr), kd


def zero_tone_wave(u, zeros, chunk=100):
    """W_K(u) = −2 Σ_k cos(γ_k u − arg ρ_k)/|ρ_k|  − ln 2π e^{−u/2}  (the first K zero-tones, rhp.wave.heartbeat plus the constant term); chunked."""
    u = np.asarray(u, float); g = np.asarray(zeros, float); rho = 0.5 + 1j * g; out = np.empty(len(u))
    for a in range(0, len(u), chunk):
        out[a:a + chunk] = -2 * np.real(np.exp(1j * np.outer(u[a:a + chunk], g)) / rho).sum(axis=1)
    return out - math.log(2 * math.pi) * np.exp(-u / 2)


def rice_rate(zeros):
    """Rice's law (ESTABLISHED for a stationary Gaussian process; here a sum of incommensurate tones, so an approximation): expected zero-crossings per unit u
    of W_K = Σ a_k cos(γ_k u − φ_k), a_k = 2/|ρ_k|:  (1/π) √( Σ a_k² γ_k² / Σ a_k² )."""
    g = np.asarray(zeros, float); a2 = 4 / (0.25 + g * g)
    return float(math.sqrt((a2 * g * g).sum() / a2.sum()) / math.pi)


def match_nodes(exact, found, tol):
    """Fraction of `exact` nodes with a node of `found` within ±tol (and the reverse): how well two node sets agree."""
    exact = np.sort(exact); found = np.sort(found)
    if len(exact) == 0 or len(found) == 0: return 0.0, 0.0
    def near(a, b):
        i = np.searchsorted(b, a); lo = np.abs(a - b[np.clip(i - 1, 0, len(b) - 1)]); hi = np.abs(a - b[np.clip(i, 0, len(b) - 1)])
        return float(np.mean(np.minimum(lo, hi) <= tol))
    return near(exact, found), near(found, exact)


# ------------------------------------------------------------------ 2-D marching squares
def zero_set_2d(x, y, F):
    """Zero set of F on the grid x (columns) × y (rows).  Returns a list of polylines (arrays of (x, y)) and a parallel list of `closed` flags.
    Marching squares; the crossing on each cell edge is the linear interpolation; a saddle cell (opposite corners agree) is resolved by the sign of the cell mean.
    Segments are stitched through shared edge crossings (an edge is identified by its two grid nodes)."""
    F = np.asarray(F, float); ny, nx = F.shape
    segs = []                                                   # (point A, point B, key A, key B)
    def cross(i, j, k):                                         # k: 0 bottom, 1 right, 2 top, 3 left edge of cell (i, j)
        a, b = {0: ((i, j), (i, j + 1)), 1: ((i, j + 1), (i + 1, j + 1)), 2: ((i + 1, j), (i + 1, j + 1)), 3: ((i, j), (i + 1, j))}[k]
        t = F[a] / (F[a] - F[b])
        return (x[a[1]] + t * (x[b[1]] - x[a[1]]), y[a[0]] + t * (y[b[0]] - y[a[0]])), (min(a, b), max(a, b))
    for i in range(ny - 1):
        for j in range(nx - 1):
            c = (F[i, j], F[i, j + 1], F[i + 1, j + 1], F[i + 1, j])          # corners 0..3 going round
            code = sum((1 << q) for q in range(4) if c[q] > 0)
            if code in (0, 15): continue
            edges = [k for k, (p, q) in enumerate(((0, 1), (1, 2), (3, 2), (0, 3))) if (c[p] > 0) != (c[q] > 0)]
            if len(edges) == 2: pairs = [tuple(edges)]
            else:
                mean_pos = (sum(c) / 4) > 0
                pairs = [(0, 3), (1, 2)] if (code == 5) == mean_pos else [(0, 1), (2, 3)]
            for p, q in pairs:
                (A, ka), (B, kb) = cross(i, j, p), cross(i, j, q)
                segs.append((A, B, ka, kb))
    ends = {}
    for n, (A, B, ka, kb) in enumerate(segs):
        ends.setdefault(ka, []).append(n); ends.setdefault(kb, []).append(n)
    used = [False] * len(segs)
    def chain(n, key):                                          # leave segment n through `key`, follow the line, return the points beyond and the last key
        pts = []
        while True:
            nxt = [m for m in ends[key] if m != n and not used[m]]
            if not nxt: return pts, key
            m = nxt[0]; used[m] = True; A, B, ka, kb = segs[m]
            if ka == key: pts.append(B); key = kb
            else: pts.append(A); key = ka
            n = m
    out, cl = [], []
    for n in range(len(segs)):
        if used[n]: continue
        used[n] = True; A, B, ka, kb = segs[n]
        fwd, kf = chain(n, kb); back, kback = chain(n, ka)
        arr = np.array(list(reversed(back)) + [A, B] + fwd)
        closed = (kf == ka) or (kback == kb)
        if closed and not np.allclose(arr[0], arr[-1]): arr = np.vstack([arr, arr[0]])
        out.append(arr); cl.append(bool(closed))
    return out, cl


def polyline_length(p):
    return float(np.sum(np.hypot(*np.diff(np.asarray(p), axis=0).T)))


def seg_intersections(polysA, polysB):
    """Exact intersection points of two sets of polylines (segment against segment)."""
    pts = []
    for A in polysA:
        a0, a1 = A[:-1], A[1:]
        for B in polysB:
            b0, b1 = B[:-1], B[1:]
            da = (a1 - a0)[:, None, :]; db = (b1 - b0)[None, :, :]; w = a0[:, None, :] - b0[None, :, :]
            den = da[..., 0] * db[..., 1] - da[..., 1] * db[..., 0]
            with np.errstate(divide="ignore", invalid="ignore"):
                s = (db[..., 0] * w[..., 1] - db[..., 1] * w[..., 0]) / den          # parameter along A
                t = (da[..., 0] * w[..., 1] - da[..., 1] * w[..., 0]) / den          # parameter along B
            ok = (np.abs(den) > 1e-14) & (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
            for i, j in zip(*np.nonzero(ok)): pts.append(a0[i] + s[i, j] * (a1[i] - a0[i]))
    return dedupe(np.array(pts)) if pts else np.zeros((0, 2))


def dedupe(pts, tol=1e-9):
    """Merge points closer than tol (a crossing at a polyline joint is reported once per adjoining segment pair)."""
    out = []
    for p in np.asarray(pts, float):
        if not any(np.hypot(*(p - q)) < tol for q in out): out.append(p)
    return np.array(out)


def sign_domains(F):
    """Nodal domains of the sampled field: connected components (4-neighbour) of the sets F > 0 and F < 0.  Returns (n_positive, n_negative, label array)."""
    from scipy import ndimage
    F = np.asarray(F)
    lp, npos = ndimage.label(F > 0); ln_, nneg = ndimage.label(F < 0)
    lab = np.where(F > 0, lp, np.where(F < 0, -ln_, 0))
    return int(npos), int(nneg), lab


def sand(F, eps=None):
    """Cymatic 'sand': density ∝ exp(−F²/(2 eps²)), heaviest on the nodal set (eps defaults to 5 % of the field's rms)."""
    F = np.asarray(F, float); eps = eps or 0.05 * float(np.sqrt(np.mean(F ** 2)))
    return np.exp(-F ** 2 / (2 * eps ** 2))
