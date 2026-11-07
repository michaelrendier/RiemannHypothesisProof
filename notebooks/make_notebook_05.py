# SPDX-License-Identifier: GPL-3.0-only
"""Builds notebooks/05_extinction_flow_nodal_ladder_zerorb.ipynb."""
import nbformat as nbf
nb = nbf.v4.new_notebook()
M = lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))
C = lambda s: nb.cells.append(nbf.v4.new_code_cell(s.strip("\n")))

M("""
# 05 — The extinction flow, its Laplace dual, the order of approach, the nodal ladder, and what 0_RB says about ζ

**Assistant note.** Written by an AI assistant (Claude) with the author, from `rhp/extinction.py`, `rhp/nodal.py`, `rhp/zerorb.py` and the figures in `rhp/atlas/`
(`flow.svg`, `net.svg`, `nodal.svg`, `flow_player.html`). Labels: ESTABLISHED, OURS, THEORETICAL, CONJECTURE, OPEN; `:CALCULATED` = this code reproduces it.
**This sets up a space and machinery; it proves nothing about RH.**

**Vantage point.** Working axis u = ln x (the Archimedes Screw's own); time of the flow τ = ln y; dual variable s′ = s − 1 (Laplace conjugate of u).
Aperture: integers ≤ 2·10⁶ for the flow; ξ on σ ∈ [−1.5, 2.5], t ≤ 50 for the ladder.
""")
C("""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd().parent))
import numpy as np, math
from IPython.display import SVG, display
from rhp import extinction as X, nodal as N, wrap as W, zerorb as Z
""")
M("## 1  The two books as continuous flows (ESTABLISHED: Buchstab ω, Dickman ρ, Mertens)")
C("""
ur, r = X.dickman(); uo, o = X.buchstab()
at = lambda u, a, x: a[int(round(x / (u[1] - u[0])))]
print("rho(2) =", at(ur, r, 2), " exact 1-ln2 =", 1 - math.log(2))
print("omega(3) =", at(uo, o, 3), " exact (1+ln2)/3 =", (1 + math.log(2)) / 3, " limit e^-gamma =", math.exp(-X.GAMMA))
N_ = 1_000_000
spf, gpf = W.spf_gpf(N_)
for U in (2, 3):
    y = N_ ** (1 / U)
    print(f"U={U}: survivors*ln y/N = {X.survivors(spf, N_, y) * math.log(y) / N_:.4f}   omega = {at(uo, o, U):.4f}   | smooth/N = {X.smooth(gpf, N_, y) / N_:.4f}  rho = {at(ur, r, U):.4f}")
display(SVG("../rhp/atlas/flow.svg"))
""")
M("""
Finite-N points sit below ω and above ρ by O(1/ln y); that is the known slow convergence, not a failure.
Open `../rhp/atlas/flow_player.html` in a browser to drag τ and watch each integer die (blue→grey) and be born (orange ring).
""")
M("""
## 2  The Laplace dual: the two books are reciprocal  (OURS-CALCULATED)
ρ̂(s) = e^{−E₁(s)}/s and 1 + ω̂(s) = e^{+E₁(s)}, so **s·ρ̂(s)·(1 + ω̂(s)) = 1** — the Laplace twin of ζ·(1/ζ) = 1 (1/ζ = Σμ(n)n⁻ˢ is the sieve).
Derived from (uω)′ = ω(u−1) and uρ′ = −ρ(u−1); checked below from the NUMERIC solutions, error first order in the grid step. Not found stated this way; may be known.
""")
C("""
for h in (1/500, 1/2000, 1/8000):
    a, b = X.dickman(14, h), X.buchstab(14, h)
    print(f"h=1/{round(1/h)}:", [f"{X.reciprocal_residual(s, a[0], a[1], b[0], b[1]):.1e}" for s in (0.5, 1.0, 0.7+3j, 1+10j)])
""")
M("## 3  Order of approach: the two paths by which a number is reached (OURS reading; the numbers are measured)")
C("""
for n in (12, 360, 49, 97, 1001):
    print(n, X.approach_paths(n), "span ln(gpf/spf) = %.3f" % X.span(n))
w = np.log(gpf[2:] / spf[2:])
print("integers <= 1e6 with span 0 (prime powers):", int(np.sum(w == 0)), "of", len(w), "; mean span of the rest %.2f e-folds" % w[w > 0].mean())
from scipy.stats import spearmanr
n_ = np.arange(2, 200002)
print("rank correlation of size with birth order (gpf): %.3f ; with death order (spf): %.3f" % (spearmanr(n_, gpf[2:200002])[0], spearmanr(n_, spf[2:200002])[0]))
""")
M("""
A number is the same multiset of primes read in two opposite orders: ascending (the sieve strips the smallest first) and descending (the wheels finish it by switching on the largest last).
They coincide only for prime powers. The size order of the integers is almost unrelated to either arrival order.
**READING, not a theorem:** the reversal of the order of approach is the candidate for what the north–south mirror s↦1−s does to a term (0_RB: R̂† = B̂).
""")
M("## 4  The six faces of (u, τ, σ′, t)")
C('display(SVG("../rhp/atlas/net.svg"))')
M("""
## 5  The nodal ladder, both hemispheres
ESTABLISHED: ξ(½+it) = F(t) is real; ∂ξ/∂σ = −iF′ there. OURS-CALCULATED here: for t ≤ 50 the curves Re ξ = 0 (through the zeros) and Im ξ = 0 (through F's extrema) alternate as a ladder of
near-horizontal rungs that cross only on the equator; the extrema are hyperbolic saddles (Hessian eigenvalues ±λ); no off-line crossing in the window tested.
""")
C("""
zs = N.rungs(50, 3000); ex = N.heartbeat_extrema(50, 3000)
print(len(zs), "zero rungs;", len(ex), "saddle rungs; one saddle between each consecutive pair:", all(sum(a < e < b for e in ex) == 1 for a, b in zip(zs[:-1], zs[1:])))
print("Hessian eigenvalues at first saddle:", N.saddle_hessian(ex[0]))
sig = np.linspace(-0.6, 1.6, 45); t = np.linspace(10, 35, 160)
print("off-line Re=Im=0 cells in sigma in [-0.6,1.6], t in [10,35]:", N.offline_crossings(sig, t, N.grid(sig, t)))
display(SVG("../rhp/atlas/nodal.svg"))
""")
M("""
## 6  What 0_RB says about ζ (located 2026-11-07) and the tests against it
Where ζ lives in 0_RB: **Σ_RB = Σ_p p^{−σ}[R̂_p⊗∂̂ + ∂̂†⊗B̂_p]** (`h_rb_hat/maths.py`) — ζ is the coupling p^{−σ}; R̂† = B̂ is the functional equation; Red = x·p, Blue = ½p² + ℘.
The **Archimedes Screw** (`archimedes_screw/maths.py`) uses u = ln x, pitch ln p, ψ jumps by ln p at u = ln p, zeros are tones. **Prime/Factor wavelength** and **Ring-theory spine** are in `.clauderc_canonical_maths`.
""")
C("""
print("screw: max |jump - ln p| over prime powers <= 200:", Z.jumps_at_ln_p(200))
for r in Z.screw_reconstruction((100.5,), (10, 50)): print("tabulated zeros  x=%.1f K=%d  exact %.3f  explicit %.3f  |err| %.3f" % r)
print("Sigma_RB balance / sum(G) at x=2, p=1.3, by sigma:", [round(v[1] / Z.load_hrb().sigma_rb_evaluate(v[0], 2.0, 1.3)['total_G'], 6) for v in Z.sigma_rb_balance_vs_sigma(2.0, 1.3)])
cut, tab = Z.dominance_table(100.0); print("cutoff sqrt(t/2pi) =", round(cut, 3)); [print(" ", k, v) for k, v in tab.items()]
for n in (2, 3, 5): print("n=%d: t=%.1f  zero spacing %.3f  wheel period %.3f  ratio %.3f  | measured spacing %.3f" % ((n,) + Z.wavelength_pair(n) + (Z.measured_spacing(n)[0],)))
""")
M("""
**Findings.** (1) the Screw's explicit formula reproduces ψ from true zeros and its jumps are exactly ln p; (2) at n ≥ 51 its tabulated zero heights fall back to a Lambert approximation and drift by ≈ 2.3 — logged;
(3) Σ_RB's `balance` is (Σ_p p^{−σ})·(E_Red − E_Blue), so it vanishes on a σ-independent set and cannot select σ = ½ — logged;
(4) in each hemisphere the forward (Red) and mirror (Blue) waves swap dominance at n = √(t/2π), in opposite senses north and south (the polarisation ratio of notebook 04);
(5) at the resonance t = 2πn² the mean zero spacing is exactly half the n-th wheel's period, so λ_prime = λ_factor/2 there.

## 7  Anomaly log (Chase Every Anomaly)
| # | Anomaly | Verdict |
|---|---|---|
| 1 | Writing `rhp/flow.py` overwrote the existing notebook-02 module | CODE (mine): restored from git; new module is `rhp/extinction.py` |
| 2 | Finite-N survivors/smooth points do not sit on ω/ρ | MATHS: O(1/ln y) known; shrinks with N |
| 3 | Reciprocity residual ~1e-4 | METHOD: first-order grid error at the jump ω(1⁻)=0→1; scales with h; closed forms agree to 1e-12 |
| 4 | Screw zeros after n = 50 off by ≈ 2.3 | CODE (ValaQuenta, not edited): Lambert approximation fallback |
| 5 | Σ_RB balance independent of σ | CODE/METHOD (ValaQuenta, not edited): σ enters only as a common weight |
""")
nb["metadata"]["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nbf.write(nb, "05_extinction_flow_nodal_ladder_zerorb.ipynb")
