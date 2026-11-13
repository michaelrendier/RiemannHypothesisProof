# SPDX-License-Identifier: GPL-3.0-only
"""One object for the tones of a zeta-type function, on either side of the boundary.   :CALCULATED.

ToneBank.zeta(K)         the Riemann side: zero ordinates γ_k, ζ'(ρ_k), the ψ-wave weights 1/ρ and the M-wave weights 1/(ρ ζ'(ρ)); conventions as in rhp.wave and rhp.mertens_wave (u = ln x).
FermatBank(n, p)         the Fermat side: the (n−1)(n−2) Frobenius angles θ_j (multiplicities kept), ω_j = −J_j = √p e^{iθ_j}; wave W_F(u) = −Σ cos(θ_j u/ln p), u = k ln p.
Both expose  .gamma (ordinates; for the finite-field bank θ_j/ln p), .mult, .wave(u) and .zero_sum(h) so the explicit-formula checks (rhp.explicit) and the nodal tools (rhp.nodal_extract, rhp.two_sided) take either.
"""
from __future__ import annotations
import math
import numpy as np
from . import mertens_wave as MW, nodal_extract as NX, factor_flow as FL


class ToneBank:
    kind = "zeta"

    def __init__(self, gamma, zp):
        self.gamma = np.asarray(gamma, float); self.zp = np.asarray(zp, complex); self.mult = np.ones(len(self.gamma), int)
        self.rho = 0.5 + 1j * self.gamma; self.c_psi = 1 / self.rho; self.c_M = MW.residues(self.gamma, self.zp)

    @classmethod
    def zeta(cls, K=None, path="data/zetaprime_25000.tsv"):
        g, z = MW.load_zp(path, K); return cls(g, z)

    def __len__(self): return len(self.gamma)

    def psi_wave(self, x, K=None):
        """(ψ(x) − x)/√x from the first K zero-tones (without the −ln 2π term's e^{-u/2} ... kept: −ln 2π/√x)."""
        K = K or len(self); x = np.asarray(x, float)
        return NX.zero_tone_wave(np.log(x), self.gamma[:K])

    def mobius_wave(self, x, K=None):
        """M_K(x)/√x (explicit formula, with the s = 0 and trivial-zero terms)."""
        K = K or len(self); x = np.asarray(x, float); return MW.m_reconstruct(x, self.gamma, self.c_M, K) / np.sqrt(x)

    def zero_sum(self, h, K=None):
        """Σ over ALL zeros (γ and −γ) of h(γ), h even and real."""
        K = K or len(self); return float(2 * np.sum(h(self.gamma[:K])))

    def resolution(self, n):
        """number of zero-tones the ping needs to resolve integers near n (empirical law K ≈ n/3)."""
        return int(math.ceil(n / 3))


class FermatBank:
    kind = "fermat"

    def __init__(self, n, p):
        self.n, self.p = int(n), int(p); self.theta = FL.frobenius_angles(p, n); self.L = math.log(p); self.period = 2 * math.pi / self.L
        key = np.round(np.mod(self.theta, 2 * math.pi), 7); uniq, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
        self.gamma = uniq / self.L; self.mult = cnt

    def __len__(self): return len(self.theta)

    def wave(self, u):
        return FL.fermat_wave(u, self.p, self.theta)

    def omegas(self):
        return math.sqrt(self.p) * np.exp(1j * self.theta)

    def affine_count(self, k):
        """predicted #{x^n + y^n = 1 over F_{p^k}} = p^k − Σ_{a,b ≥ 1} (−J_p(a,b))^k  (Hasse–Davenport), with the trivial a + b ≡ 0 terms included."""
        from .fermat_curve import jacobi_sums
        J = jacobi_sums(self.p, self.n); return float(self.p ** k - sum((-v) ** k for v in J.values()).real)
