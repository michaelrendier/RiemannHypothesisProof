# SPDX-License-Identifier: GPL-3.0-only
"""The Heartbeat waveform, derived from the two sets (the Two Trees partition), and its 'wind' statistics.

THE TWO SETS (engine: GenerationalLineage `two_trees`, AbrikosovTree/engine/two_trees.py):  every integer is prime (Telperion), composite (Laurelin) or
one of {0,1} (Mingling); B(n) + R(n) + M(n) = 1 exactly.  A partition forces the two deviations from their smooth parts to be EXACTLY opposite.

STEP 1 (weights).  Give each integer n one unit of weight and split it between the trees:  Telperion takes Λ(n) (ln p on a prime power p^k, else 0),
  Laurelin takes the remainder 1 − Λ(n).  Cumulative:  T(x) = Σ_{n≤x} Λ(n) = ψ(x),   L(x) = Σ_{2≤n≤x}(1 − Λ(n)) = ⌊x⌋ − 1 − ψ(x).
  So  T(x) + L(x) = ⌊x⌋ − 1  exactly (the engine's B+R+M=1 in weighted form).  Smooth parts: T̄ = x, L̄ = −{x}-sized (mean 0).
STEP 2 (native variable).  u = ln x, and divide by e^{u/2} = √x (the ½ of the line is this flare):  W_T(u) = (ψ(x) − x)/√x,  W_L(u) = (L(x) − 0)/√x ≈ −W_T.
  W_T + W_L = (⌊x⌋ − 1 − x)/√x = −({x} + 1)/√x → 0.   Exact antisymmetry up to a term ≤ 2/√x.
STEP 3 (the zero expansion; ESTABLISHED, von Mangoldt):  ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²), so with ρ = ½ + iγ
  W_T(u) = − Σ_ρ e^{iγu}/ρ  − ln 2π · e^{−u/2} − …  = − 2 Σ_{γ>0} Re[e^{iγu}/ρ] = −2 Σ_{γ>0} [½ cos γu + γ sin γu]/|ρ|².
  Each zero pair is TWO equal-radius circles of radius 1/|ρ|, turning opposite ways (e^{+iγu}/ρ and its conjugate): a standing wave; their sum is the real tone
  −2cos(γu − arg ρ)/|ρ|, a sine of amplitude 2/|ρ| with a phase set by the ½ (arg ρ = π/2 − arctan(1/2γ)).
THE FOUNDATIONAL FUNCTION:  H(u) = −2 Σ_{γ>0} cos(γu − arg ρ_γ)/|ρ_γ|;   Telperion carries +H, Laurelin carries −H.
Statistics (ESTABLISHED under the assumption the γ are linearly independent): mean power ⟨H²⟩ = Σ_ρ 1/|ρ|² = 2 + γ_E − ln 4π = 0.046191.
Labels: steps 1–2 are identities (ESTABLISHED, trivial); step 3 is ESTABLISHED; the numbers below are OURS-CALCULATED.  'Wind' / 'gust' are OURS readings.
"""
from __future__ import annotations

import math
import pathlib
import sys
import numpy as np

from .heartbeat import MEAN_POWER
from .unsmap import analytic_signal

ENGINE = pathlib.Path(__file__).resolve().parents[2] / "GenerationalLineage"


def engine_two_trees(N=100_000):
    """The GenerationalLineage engine's own `two_trees(N)` (the reference the partition is checked against)."""
    if str(ENGINE) not in sys.path:
        sys.path.insert(0, str(ENGINE))
    from engine.lineage import two_trees
    return two_trees(N)


def own_partition(N):
    """(primes, composites, mingling) counted independently of the engine."""
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    p = int(s.sum())
    return p, N + 1 - p - 2, 2


def lam_table(X):
    """Λ(n) for n = 0..X (von Mangoldt) and the cumulative ψ."""
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
    return lam, np.cumsum(lam)


def trees(u, psi):
    """Telperion's and Laurelin's waves from the SETS (no zeros): W_T = (ψ(x) − x)/√x,  W_L = (⌊x⌋ − 1 − ψ(x))/√x."""
    u = np.asarray(u, float)
    x = np.exp(u)
    fl = np.minimum(np.floor(x).astype(np.int64), len(psi) - 1)
    ps = psi[fl]
    return (ps - x) / np.sqrt(x), (fl - 1 - ps) / np.sqrt(x)


def antisymmetry_defect(u, psi):
    """max |W_T + W_L| and the bound max 2/√x: the two trees' deviations cancel up to the floor term."""
    wt, wl = trees(u, psi)
    return float(np.max(np.abs(wt + wl))), float(np.max(2 * np.exp(-np.asarray(u) / 2)))


def tone(u, gamma):
    """One zero pair's real standing tone: −2 cos(γu − arg ρ)/|ρ|."""
    u = np.asarray(u, float)
    rho = complex(0.5, gamma)
    return -2 * np.cos(gamma * u - np.angle(rho)) / abs(rho)


def heartbeat(u, zeros):
    """H(u) = Σ_{γ in zeros} tone(u, γ)  (the first len(zeros) zero pairs); excludes the e^{−u/2} ln 2π remainder."""
    u = np.asarray(u, float)
    g = np.asarray(zeros, float)
    rho = 0.5 + 1j * g
    out = np.empty(len(u))
    for a in range(0, len(u), 2000):                      # chunked: 2000 × K complex at a time
        out[a:a + 2000] = -2 * np.real(np.exp(1j * np.outer(u[a:a + 2000], g)) / rho).sum(axis=1)
    return out


def reconstruction_error(u, psi, zeros):
    """(correlation, rms of W_T − (H − ln2π e^{−u/2})) between the sieve wave and the zero-sum wave at the sample points (u between jumps)."""
    wt, _ = trees(u, psi)
    h = heartbeat(u, zeros) - math.log(2 * math.pi) * np.exp(-np.asarray(u) / 2)
    c = float(np.corrcoef(wt, h)[0, 1])
    return c, float(np.sqrt(np.mean((wt - h) ** 2)))


def wind_stats(w):
    """Moments of a wave: rms, skew, excess kurtosis, peak |w|, and gust factor peak/rms."""
    w = np.asarray(w, float)
    m = w.mean(); sd = w.std()
    z = (w - m) / sd
    return dict(mean=float(m), rms=float(math.sqrt(np.mean(w ** 2))), sd=float(sd), skew=float(np.mean(z ** 3)),
                ex_kurt=float(np.mean(z ** 4) - 3), peak=float(np.max(np.abs(w))), gust=float(np.max(np.abs(w)) / math.sqrt(np.mean(w ** 2))))


def envelope(w):
    """Instantaneous amplitude of the wave (Hilbert analytic signal): the 'wind speed' with the carrier removed."""
    return analytic_signal(w)[0]


def factor_step_ratio(p, rms):
    """A prime p contributes a jump ln p to ψ at x = p; in W units that is ln p/√p.  Returns (step in W units, step/rms)."""
    s = math.log(p) / math.sqrt(p)
    return s, s / rms


def exceedance(w, thr):
    """Fraction of the sample with |w| ≥ thr: how often a gust is at least a factor's own step."""
    return float(np.mean(np.abs(np.asarray(w)) >= thr))


def jump_fraction(p, zeros, K):
    """Fraction of the true jump ln p of ψ at x = p recovered by the first K zero pairs: [S_K(p+½) − S_K(p−½)]/ln p … read across a unit window,
    the smooth part (+1 from x) removed:  jump_K = 1 − Σ_{k≤K} Δ(x^ρ/ρ)·2Re."""
    g = np.asarray(zeros[:K], float)
    rho = 0.5 + 1j * g
    d = 2 * np.real(((p + 0.5) ** rho - (p - 0.5) ** rho) / rho).sum()
    return float((1.0 - d) / math.log(p))
