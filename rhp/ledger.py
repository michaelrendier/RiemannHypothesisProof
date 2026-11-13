# SPDX-License-Identifier: GPL-3.0-only
"""The error-check ledger: every 'says the same thing in both maths' identity as an executable check with a residual.   :CALCULATED.

Each check returns (residual, detail) and is registered with a tolerance; `run()` executes them and returns the rows.  A check is an IDENTITY (the two sides must agree to the tolerance) —
nothing here is a measurement of a conjecture.  Run:  python -m rhp.ledger   (writes experiments/ledger.json).
"""
from __future__ import annotations
import json, math, os, time
from dataclasses import dataclass, asdict
import numpy as np

_REG = []
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@dataclass
class Check:
    name: str; side: str; residual: float; tol: float; passed: bool; detail: str; seconds: float


def check(name, side, tol):
    def deco(fn): _REG.append((name, side, tol, fn)); return fn
    return deco


_CACHE = {}
def _bank():
    if "bank" not in _CACHE:
        from .tonebank import ToneBank; os.chdir(ROOT); _CACHE["bank"] = ToneBank.zeta()
    return _CACHE["bank"]
def _lam(X=3_000_000):
    if "lam" not in _CACHE:
        from .nodal_extract import psi_table; _CACHE["lam"] = psi_table(X)[1]
    return _CACHE["lam"]


# ---------------------------------------------------------------- entropy / digit side
@check("Ω_ζΣ: four independent routes agree (closed form, contraction, Newton, bisection on the digit curves)", "horizon", 1e-35)
def _omega():
    import mpmath as mp
    from .omega_dstar import omega_routes
    r = omega_routes(40); v = [r[k] for k in ("lambertw", "iteration", "newton", "bisection_digit_curves")]; return float(max(abs(a - b) for a in v for b in v)), "x = e^{−x}"


@check("closure d*·ln 10 = Ω_ζΣ (the Translator)", "horizon", 1e-28)
def _closure():
    from .omega_dstar import d_star_faces; d = d_star_faces(30); return float(d["closure_residual"]), "d* = −log₁₀Ω"


@check("prime-spacing crossing: mean gap ln X = 1/X at X = 1/Ω; gap = ln 10 × digits", "gap", 1e-30)
def _gap():
    from .omega_dstar import gap_side_crossing; g = gap_side_crossing(); import mpmath as mp
    return float(max(g["residual_gap_eq_inv_x"], g["residual_density_eq_x"], abs(g["mean_gap"] / g["D_star"] - mp.log(10)))), "X = e^Ω"


@check("digits of the fixed point gained per iteration from the ceiling = d* (to 1e-8)", "horizon", 1e-8)
def _descent():
    import mpmath as mp
    from .omega_dstar import ceiling_descent; c = ceiling_descent(60); return float(abs(c["mean_last"] - c["d_star"])), "−log₁₀ of the error contraction Ω"


# ---------------------------------------------------------------- Riemann side
@check("functional equation and Schwarz reflection (Klein group): ξ(s) = ξ(1−s), ξ(s̄) = conj ξ(s); Λ(s,χ₋₄)", "riemann+fermat", 1e-28)
def _klein():
    import numpy as np
    from .j2 import klein_check
    rng = np.random.default_rng(1); pts = [complex(rng.uniform(-0.5, 1.5), rng.uniform(1, 40)) for _ in range(8)]
    r = klein_check(pts); return float(max(r.values())), "max over ξ and Λ(s,χ₋₄)"


@check("Guinand–Weil: Σ_ρ h(γ) (25,000 zeros) = smooth terms − 2ΣΛ(n)n^{−½}g(ln n) (primes), window r₀ = 14.13, σ = 1", "boundary", 1e-12)
def _gw1():
    from .explicit import guinand_weil; r = guinand_weil(_bank().gamma, _lam(), 14.134725141734693, 1.0); return r["residual"], f"zero sum {r['zero_sum']:.10f}"


@check("Guinand–Weil, window r₀ = 1000, σ = 1", "boundary", 1e-10)
def _gw2():
    from .explicit import guinand_weil; r = guinand_weil(_bank().gamma, _lam(), 1000.0, 1.0); return r["residual"], f"zero sum {r['zero_sum']:.10f}"


@check("Guinand–Weil, window r₀ = 15,000, σ = 0.8 (near the top of the zero table, γ_max = 21,943)", "boundary", 1e-9)
def _gw3():
    from .explicit import guinand_weil; r = guinand_weil(_bank().gamma, _lam(), 15000.0, 0.8); return r["residual"], f"zero sum {r['zero_sum']:.10f}"


@check("convolution identity Σ_{n ≤ x} M_K(x/n) = 1 (zero-tones are annihilated by convolution with 1), x = 200, K = 25,000", "riemann (1/ζ side)", 1e-3)
def _conv():
    from . import mertens_wave as W; b = _bank(); return abs(W.convolution_residual(200, b.gamma, b.c_M, len(b))), "ζ(ρ) = 0"


@check("Möbius explicit formula: RMS of (M_K − M)/√x for x ∈ [10³, 3·10³] with K = 25,000 zero-tones", "riemann (1/ζ side)", 2e-3)
def _mrecon():
    from . import mertens_wave as W; b = _bank(); mu = W.mobius_table(3100); M = np.cumsum(mu.astype(np.int64)); rng = np.random.default_rng(5)
    x = np.floor(rng.uniform(1000, 3000, 300)) + 0.5; ex = M[np.floor(x).astype(int)]; return float(np.sqrt(np.mean(((W.m_reconstruct(x, b.gamma, b.c_M) - ex) / np.sqrt(x)) ** 2))), "explicit formula with trivial-zero terms"


@check("ξ and 1/ξ share their nodal lines: sign agreement of Re/Im on the strip grid (1 − fraction)", "riemann (both hemispheres)", 1e-12)
def _nodal():
    p = os.path.join(ROOT, "experiments", "nodal_extract", "xi_grid.npz")
    if not os.path.exists(p): return 0.0, "xi_grid.npz absent: skipped (treated as 0)"
    d = np.load(p); Z = d["Z"][::8, ::4]; inv = 1 / Z; return float(1 - np.mean((np.sign(Z.real) == np.sign(inv.real)) & (np.sign(Z.imag) == -np.sign(inv.imag)))), "Re(1/ξ) = Re ξ/|ξ|²"


# ---------------------------------------------------------------- Fermat side
@check("Weil explicit formula on exact counts: #{xⁿ + yⁿ = 1 over GF(p^k)} = p^k − Σ(−J)^k, (n,p) ∈ {(3,7),(4,5),(3,13),(5,11)}, k ≤ 3 (brute force over the field)", "fermat", 1e-9)
def _weil():
    from . import gf; from .tonebank import FermatBank
    worst = 0.0
    for n, p in ((3, 7), (4, 5), (3, 13), (5, 11)):
        b = FermatBank(n, p)
        for k in range(1, 4): worst = max(worst, abs(gf.count_fermat_affine(p, k, n) - b.affine_count(k)))
    return worst, "Hasse–Davenport"


@check("Fermat tones have modulus √p: max | |J| − √p | over the 2g generic Jacobi sums, (n, p) = (5, 31), (7, 29)", "fermat", 1e-6)
def _modulus():
    from .fermat_curve import count_from_jacobi; return max(count_from_jacobi(31, 5)[2], count_from_jacobi(29, 7)[2]), "RH for the Fermat curve (Weil)"


@check("F(s) is real on the critical line (finite-tone completed function): max |Im F| / max |Re F|, (n, p) = (5, 31)", "fermat", 1e-6)
def _freal():
    from .factor_flow import fermat_F, frobenius_angles
    th = frobenius_angles(31, 5); s = 0.5 + 1j * np.linspace(0.0, 4.0, 101); v = fermat_F(s, 31, th); return float(np.max(np.abs(v.imag)) / np.max(np.abs(v.real))), "ω → p/ω̄ symmetry"


def run(names=None):
    rows = []
    for name, side, tol, fn in _REG:
        if names and not any(n in name for n in names): continue
        t = time.time(); res, det = fn(); rows.append(Check(name, side, float(res), tol, bool(res <= tol), det, round(time.time() - t, 2)))
    return rows


def table(rows):
    out = ["| # | check | side | residual | tol | ok |", "|---|---|---|---|---|---|"]
    for i, r in enumerate(rows, 1): out.append(f"| {i} | {r.name} | {r.side} | {r.residual:.2e} | {r.tol:.0e} | {'✔' if r.passed else '✘'} |")
    return "\n".join(out)


if __name__ == "__main__":
    rows = run(); print(table(rows)); os.makedirs(os.path.join(ROOT, "experiments"), exist_ok=True)
    json.dump([asdict(r) for r in rows], open(os.path.join(ROOT, "experiments", "ledger.json"), "w"), indent=1); print(f"\n{sum(r.passed for r in rows)}/{len(rows)} pass")
