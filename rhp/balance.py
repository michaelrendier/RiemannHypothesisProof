# SPDX-License-Identifier: GPL-3.0-only
"""The modulus balance of the two halves of the approximate functional equation, off the line.

ζ(s) ≈ A + χB.  D(σ,t) = ln|χ| + ln|B| − ln|A|.  On σ = ½: D ≡ 0 (B = Ā, |χ| = 1).  Off the line D/(½−σ) = ln(t/2π)+O(1) + c(σ,t),
with c driven by the logarithmic derivatives of the partial sums: it can vanish only where A or B nearly vanishes.
NOTE (chased): for the EXACT functions ζ(s) = χ(s)ζ(1−s) makes |ζ(s)| = |χ||ζ(1−s)| an identity; the balance below concerns the TRUNCATED halves.
"""
from __future__ import annotations
import numpy as np


def ln_abs_chi(sig, t):
    """ln|χ(σ+it)| for t ≥ 20 via Stirling for lnΓ(1−s):  χ = 2^s π^{s−1} sin(πs/2) Γ(1−s)."""
    s = sig + 1j * t
    z = 1 - s
    lg = (z - 0.5) * np.log(z) - z + 0.5 * np.log(2 * np.pi) + 1 / (12 * z) - 1 / (360 * z**3) + 1 / (1260 * z**5)
    ln_sin = (np.pi * t / 2) - np.log(2.0)                    # |sin(πs/2)| = ½ e^{πt/2}(1 + O(e^{−πt}))
    return np.real(s * np.log(2) + (s - 1) * np.log(np.pi) + lg) + ln_sin


def halves(sig, t):
    """A, B (complex) for arrays sig, t (broadcast). N = ⌊√(t/2π)⌋ terms."""
    sig, t = np.broadcast_arrays(np.asarray(sig, float), np.asarray(t, float))
    s = sig + 1j * t
    N = np.maximum(1, np.floor(np.sqrt(t / (2 * np.pi))).astype(int))
    A = np.zeros(s.shape, complex); B = np.zeros(s.shape, complex)
    for n in range(1, int(N.max()) + 1):
        m = N >= n
        A[m] += np.exp(-s[m] * np.log(n)); B[m] += np.exp((s[m] - 1) * np.log(n))
    return A, B


def D_field(sig, t):
    A, B = halves(sig, t)
    return ln_abs_chi(sig, t) + np.log(np.abs(B)) - np.log(np.abs(A)), A, B


def margin_field(sig, t):
    """M = D/(½−σ): positive everywhere in σ<½ ⇒ no off-line balance crossing."""
    D, A, B = D_field(sig, t)
    return D / (0.5 - np.asarray(sig, float)), A, B


def offline_crossings(ts, sigmas):
    """Sign changes of M(σ, t) along σ ∈ sigmas (σ < ½) for each t. Returns list of (t, σ*, |A+χB|, phase mismatch)."""
    out = []
    for t in ts:
        M, A, B = margin_field(sigmas, np.full(len(sigmas), t))
        sc = np.nonzero(np.sign(M[:-1]) != np.sign(M[1:]))[0]
        for i in sc:
            # linear interpolation of the root
            w = M[i] / (M[i] - M[i + 1]); s_star = sigmas[i] + w * (sigmas[i + 1] - sigmas[i])
            Mi, Ai, Bi = margin_field(np.array([s_star]), np.array([t]))
            chi = np.exp(ln_abs_chi(s_star, t))                # modulus only
            resid = abs(Ai[0]) - chi * abs(Bi[0])
            out.append((float(t), float(s_star), float(abs(Ai[0])), float(abs(Bi[0]) * chi)))
    return out
