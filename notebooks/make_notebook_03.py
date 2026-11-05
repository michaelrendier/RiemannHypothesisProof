# SPDX-License-Identifier: GPL-3.0-only
"""Builds notebooks/03_clock_echo_harness_margins_primegame.ipynb."""
import nbformat as nbf
nb = nbf.v4.new_notebook()
M = lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))
C = lambda s: nb.cells.append(nbf.v4.new_code_cell(s.strip("\n")))

M("""
# 03 — The θ-clock, the echo spectrum, an operator harness, margins, balance, and PRIMEGAME

**Assistant note.** Written by an AI assistant (Claude) with the author, from `rhp/`. Labels: ESTABLISHED, OURS, THEORETICAL, CONJECTURE, OPEN;
`:CALCULATED` = this code reproduces it. Floating-point checks are not proofs. **Standing protocol (Chase Every Anomaly):** every result that was legal
but unexpected was chased to a verdict; they are collected in §9.

**Vantage point.** Orientation: the θ-clock (phase variable) beside t; log-polar for the pathway. Aperture: σ = ½ continuous for t ≤ 3600
(3066 zeros), plus σ ∈ (0, ½) for the balance field at t ≤ 400.
""")

C("""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd().parent))
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"svg.fonttype": "none", "font.family": "monospace", "font.size": 8})
from rhp import bigpath, clock, gram, echo, harness, lehmer, balance, primegame
zeros = np.loadtxt("../data/zeros_1_to_3600.csv")                 # γ_1 … γ_3066 < 3600
zz = zeros[zeros > 20]
p = bigpath.pathway_big(3600.0, 0.02)
print(len(zeros), "zeros")
""")

M("""
## 1. The θ-clock: one circle carries the zero density, the rest only shift it

ζ ≈ A + χB and the n = 1 pair is 1 + e^{−2iθ}: a **circle of radius 1 centred at 1 through the origin**, reached when θ ≡ π/2 (mod π). Its hits are the *carrier*
E_k with θ(E_k) = (k − 3/2)π. The wheels n ≥ 2 move the kth zero by S(γ_k): θ(γ_k) = π(k − 3/2 − S). Distances below are in mean local spacings.
""")
C("""
E = clock.carrier_zeros(len(zeros))
shift = clock.shift_in_spacings(zeros)
S = clock.jitter(zeros)
print(f"carrier vs zeros: mean |γ−E| = {np.abs(shift).mean():.3f} spacings, sd {shift.std():.3f}, max {np.abs(shift).max():.3f}")
print(f"S(γ_k): mean {S.mean():+.4f}, sd {S.std():.4f}, range [{S.min():+.3f}, {S.max():+.3f}];   max |S + shift·(θ'/π·spacing)| check: {np.max(np.abs(S - (-shift))):.2e}")
print("zeros are the carrier's hits, displaced by < 1 spacing everywhere:", bool(np.abs(shift).max() < 1.0))
""")

M("""
## 2. Where the pathway crosses the real axis — and the negative side

Im ζ = 0 only at Gram points θ(gₙ) = nπ (or zeros), where ζ = (−1)ⁿ Z(gₙ): positive when Gram's law holds, **negative when it fails**.
""")
C("""
r = gram.negative_crossings(3600.0)
n_, g_, v_, neg = r["n"], r["g"], r["value"], r["neg"]
print("Gram points below 3600:", len(g_), "| negative real-axis crossings (Gram's law fails):", len(neg), f"({len(neg)/len(g_):.3f})")
print("first ten failures (n, gₙ):", [(int(n_[i]), round(float(g_[i]), 3)) for i in neg[:10]])
print("smallest |ζ(½+i gₙ)| over all Gram points:", float(np.abs(v_).min()), "(refined with mpmath where < 0.02)")
dec = np.array_split(np.arange(len(g_)), 6)
print("failure fraction by height (six equal blocks of Gram points):", [f"{np.mean(v_[d] < 0):.3f}" for d in dec])
cnt = gram.gram_interval_counts(n_, g_, zeros)
bad = np.nonzero(cnt != 1)[0]
fail = set(int(i) for i in neg)
adj = set()
for i in fail: adj |= {i - 1, i}
print("Gram intervals with ≠ 1 zero:", len(bad), "| intervals adjacent to a failing point:", len(adj & set(bad.tolist())), "of", len(adj))
print("intervals with ≠1 zero that are NOT adjacent to a failure:", sorted(set(bad.tolist()) - adj)[:10])
""")

M("""
## 3. The echo spectrum: the zeros' Fourier sum is the explicit formula

F_Γ(t) = Σ_{γ>0} e^{−(γ/Γ)²/2} cos(γt) is the (windowed) trace of cos(Ht) for a hypothetical H with spectrum {γ}. Weil's explicit formula predicts it exactly:
the pole term cosh(t/2) (the **path of least primes**), a smooth term from ψ(¼ + ir/2), and **negative Gaussian bumps of area ½·Λ(n)/√n at t = ln n**.
The comparison below uses all 3066 zeros, Γ = 300.
""")
C("""
Gam = 300.0
t_e = np.linspace(0.3, 4.1, 3800)
F = echo.zero_sum(t_e, zeros, Gam)
P, prime_part, smooth_part = echo.explicit_prediction(t_e, Gam)
print(f"max |F − prediction| = {np.abs(F - P).max():.2e}   (max |F| = {np.abs(F).max():.1f})")
lam = echo.lam_table(80)
print("\\n  t=ln n   n    Λ(n)/√n   predicted bump height −½Λ/√n·Γ/√(2π)   measured (F − smooth part) at t")
for n in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, 49):
    x = np.log(n); i = np.argmin(abs(t_e - x))
    print(f"  {x:6.3f}  {n:3d}  {lam[n]/np.sqrt(n):8.4f}   {-0.5*lam[n]/np.sqrt(n)*Gam/np.sqrt(2*np.pi):12.2f}              {(F - smooth_part)[i]:12.2f}")
""")

M("""
## 4. The operator harness

A candidate H is a spectrum. Four scores against the zeros: **weyl** (max deviation from the smooth count θ/π + 1 − ½ + (k−1)… = 0 for the carrier), **spacing** (variance of unfolded
neighbour spacings: GUE ≈ 0.18, Poisson 1, a ladder 0), **numvar** (number variance of 5 spacings) and **locate** (mean distance to the actual zeros, in spacings).
Controls: the carrier; a ladder (n+½)ħ — the spectrum of ValaQuenta's H_NN as coded, and of xp on a periodic u-box; Poisson; a GUE matrix unfolded onto the zeros' density.
A real candidate (a boundary condition on xp, or a Σ_RB construction) must beat the carrier on `locate` **and** the GUE on `locate` while keeping the statistics.
""")
C("""
cands = harness.candidates(zeros)
rows = {k: harness.score(v, zeros) for k, v in cands.items()}
print(f"{'candidate':52s} {'weyl':>8s} {'spacing':>8s} {'numvar5':>8s} {'locate':>8s}")
for k, s in rows.items():
    print(f"{k:52s} {s['weyl']:8.3f} {s['spacing_var']:8.3f} {s['numvar5']:8.3f} {s['locate']:8.3f}")
""")

M("""
## 5. Bifurcation margins (Lehmer)

For each gap between consecutive zeros, the height of Z's extremum relative to the median extremum of the 100 gaps around it. A margin near 0 is a **near-double zero**: a small real
perturbation could merge the pair and split it off the line.
""")
C("""
mg = lehmer.gap_margins(p, zz)
order = np.argsort(mg["rel"])[:8]
print("   γ_n        gap (spacings)   |Z| extremum   relative height")
for j in order:
    print(f"{mg['gamma'][j]:10.4f}   {mg['gap_norm'][j]:10.3f}      {mg['peak'][j]:10.4f}   {mg['rel'][j]:10.4f}")
j = order[0]
from mpmath import zetazero
k_idx = int(np.searchsorted(zeros, mg["gamma"][j]))            # 0-based index of γ_n in the full list
g1, g2 = float(zetazero(k_idx + 1).imag), float(zetazero(k_idx + 2).imag)
t_star, zpk = lehmer.refine_peak(g1, g2)
print(f"\\nchase of the smallest margin: zeros n = {k_idx+1}, {k_idx+2} (mpmath zetazero): {g1:.6f}, {g2:.6f}; gap = {g2-g1:.5f}")
print(f"mpmath golden-section extremum INSIDE the gap: t* = {t_star:.5f}, |Z(t*)| = {zpk:.5f}   (Riemann–Siegel grid: {mg['peak'][j]:.5f})")
""")

M("""
## 6. The modulus balance of the two halves, off the line

D(σ,t) = ln|χ| + ln|B| − ln|A| on the **truncated** halves (N = ⌊√(t/2π)⌋ terms). On σ = ½ it vanishes identically (B = Ā). The ratio M = D/(½−σ) = ln(t/2π) + c(σ,t) is positive
unless c ≈ −ln(t/2π). The crude bound |c| ≤ 2 ln N would make that impossible, but 2 ln N = ln(t/2π) − ln(t/2πN²), so **the crude bound has zero slack at t = 2πN²**,
the heights at which a new wheel switches on. We search for sign changes of M for σ ∈ (0.02, 0.49), t ∈ [20, 400].
""")
C("""
sig = np.linspace(0.02, 0.49, 48)
ts = np.arange(20, 400, 0.02)
cr = balance.offline_crossings(ts, sig)
arr = np.array(cr)
print("off-line sign changes of M in [20,400] × (0.02,0.49):", len(cr))
# group into branches: successive crossings continuous in t
br = []
cur = [arr[0]]
for row in arr[1:]:
    if row[0] - cur[-1][0] <= 0.021 and abs(row[1] - cur[-1][1]) < 0.05: cur.append(row)
    else: br.append(np.array(cur)); cur = [row]
br.append(np.array(cur))
print("distinct branches:", len(br))
for b in br[:10]:
    print(f"  t {b[0,0]:7.2f}–{b[-1,0]:7.2f}  σ* {b[:,1].min():.3f}–{b[:,1].max():.3f}  (points {len(b)})")
print("heights 2πN²:", [round(2*np.pi*N*N, 1) for N in range(3, 8)])
from mpmath import mp, mpc, zeta
dz = np.array([min(abs(zeros - b[len(b)//2, 0])) for b in br])
loc_sp = np.array([2*np.pi/np.log(b[len(b)//2, 0]/(2*np.pi)) for b in br])
rng = np.random.default_rng(5)
tr_ = rng.uniform(20, 400, 20000)
null = np.array([np.min(np.abs(zeros - x)) for x in tr_]) / (2*np.pi/np.log(tr_/(2*np.pi)))
print(f"median distance (local spacings) from a branch midpoint to the nearest zero: {float(np.median(dz/loc_sp)):.3f}   |  random t in [20,400]: {float(np.median(null)):.3f}  (110 branches; null sd of the median ≈ {float(null.std()/np.sqrt(len(dz))*1.25):.3f})")
zeta_vals = [abs(complex(zeta(mpc(b[len(b)//2, 1], b[len(b)//2, 0])))) for b in br[:12]]
print("|ζ(σ*+it)| (exact, mpmath) on the first 12 branch midpoints:", [round(x, 3) for x in zeta_vals])
""")

M("""
## 7. PRIMEGAME and the fret picture

Conway's 14 fractions: start at 2; multiply by the first fraction that gives an integer. **Listing 1** is the generator (extracted from `rhp/primegame.py`).
""")
C("""
import inspect
print(inspect.getsource(primegame.primegame))
g = primegame.primegame(); print([next(g) for _ in range(14)])
""")
C("""
ev = primegame.events(25)
print("p, step of 2^p, steps+uses(15/14), cumulative uses of 15/14:")
claimed = {5: 281, 7: 710, 11: 2375}
for pr, st, use in ev[:9]:
    print(f"{pr:3d} {st:8d} {st+use[11]:8d} {use[11]:5d}", "   <- video:", claimed.get(pr, ""))
slope, coef = primegame.fit_cost(ev)
print(f"\\ncost fit over p ≥ 11 (first 25 primes): steps ≈ {coef:.2f} · p^{slope:.3f}")
tr = primegame.trace(2400)
print("registers (2,3,5,7 | 11,13,17,19,23,29): max exponent over 2400 steps:", tr.max(axis=0).tolist(), "| flags only ever take values", sorted(set(tr[:, 4:].ravel().tolist())))
print("\\nrotation by d on a circle of n positions closes early (period < n) iff gcd(n,d) > 1:")
for n in (12, 13, 15, 17):
    ct = primegame.closing_table(n)
    print(f"  n={n:2d}: early-closing d = {[d for d,per in ct.items() if per < n]}")
""")

M("""
## 8. Figures
""")
C("""
fig, ax = plt.subplots(2, 3, figsize=(16.8, 9.4))
# (a) Z against θ/π: the carrier 2cos θ and the modulation
m = (p["t"] > 20) & (p["t"] < 200)
ax[0, 0].plot(p["theta"][m] / np.pi, p["Z"][m], lw=0.7, color="#3b528b", label="Z(t) against θ/π")
ax[0, 0].plot(p["theta"][m] / np.pi, 2 * np.cos(p["theta"][m]), lw=0.5, color="#c0392b", label="n=1 carrier 2cos θ")
ax[0, 0].legend(frameon=False); ax[0, 0].set_title("a  the θ-clock: zeros at half-integers of θ/π, jittered by S")
# (b) shift distribution
ax[0, 1].hist(shift, bins=60, color="#21918c"); ax[0, 1].set_title("b  zero − carrier, in spacings")
# (c) Gram failures vs height
ax[0, 2].plot(g_, np.where(v_ > 0, 0, 1), ".", ms=2, color="k"); ax[0, 2].set_title("c  Gram points: 1 = negative real-axis crossing")
# (d) echo spectrum
ax[1, 0].plot(t_e, F, lw=0.6, color="k", label="Σ cos(γ t) (3066 zeros, Γ=300)"); ax[1, 0].plot(t_e, P, lw=0.5, color="#c0392b", ls="--", label="explicit formula")
for n in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, 49): ax[1, 0].axvline(np.log(n), lw=0.3, color="#3b528b")
ax[1, 0].set_ylim(-50, 60); ax[1, 0].legend(frameon=False, fontsize=7); ax[1, 0].set_title("d  echoes at t = ln(prime power)")
# (e) spacing histograms
for k, c in (("zeros (reference)", "k"), ("GUE (unfolded)", "#3b528b"), ("Poisson", "#c0392b")):
    u = harness.unfold(cands[k]); ax[1, 1].hist(np.diff(u), bins=50, range=(0, 4), histtype="step", color=c, label=k, density=True)
ax[1, 1].legend(frameon=False, fontsize=7); ax[1, 1].set_title("e  unfolded spacings")
# (f) PRIMEGAME registers
ax[1, 2].imshow(tr[:600].T, aspect="auto", cmap="viridis", interpolation="nearest")
ax[1, 2].set_yticks(range(10)); ax[1, 2].set_yticklabels(primegame.REGISTERS); ax[1, 2].set_title("f  PRIMEGAME registers, first 600 steps (rows = prime)")
fig.tight_layout(); fig.savefig("../figures/clock_echo_harness.svg"); plt.close(fig)
fig, ax = plt.subplots(1, 3, figsize=(16.8, 5))
for b in br[:40]: ax[0].plot(b[:, 1], b[:, 0], lw=0.8)
ax[0].plot(np.full(int((zeros < 400).sum()), 0.5), zeros[zeros < 400], "k.", ms=2); ax[0].set_xlabel("σ"); ax[0].set_ylabel("t"); ax[0].set_title("off-line balance crossings of the truncated halves (curves), zeros (black)")
ax[1].semilogy(mg["gamma"], mg["rel"], ".", ms=2, color="k"); ax[1].set_title("relative extremum height of Z in each gap (Lehmer margins)")
ax[2].loglog([e[0] for e in ev], [e[1] for e in ev], "o-", ms=3); ax[2].set_title(f"PRIMEGAME steps vs p: slope {slope:.2f}")
fig.tight_layout(); fig.savefig("../figures/balance_margins_primegame.svg"); plt.close(fig)
print("wrote figures/clock_echo_harness.svg and balance_margins_primegame.svg")
""")

M("""
## 9. Anomaly log (Chase Every Anomaly)
""")
C("""
# A5: are the zeros' spacing variance 0.149 < GUE 0.18 and Σ²(5) = 0.34 < 0.49 a finite-height effect? Compare low and high halves of the list.
u = harness.unfold(zeros)
lo, hi = u[:1500], u[1500:]
print(f"spacing variance: zeros 1–1500: {np.diff(lo).var():.3f}   zeros 1501–3066: {np.diff(hi).var():.3f}   (GUE: 0.180)")
print(f"number variance (L=5): zeros 1–1500: {harness.numvar(lo):.3f}   zeros 1501–3066: {harness.numvar(hi):.3f}   (GUE: ~0.51)")
""")
M("""
| # | what was anomalous | verdict | evidence |
|---|---|---|---|
| A1 | A first quick check of the zeros' Fourier sum gave "peaks" of inconsistent size and sign at ln p | **METHOD** (grid vs window) | the Gaussian window Γ = 2500 has width 4·10⁻⁴ in t, finer than the 10⁻³ grid. With Γ = 300 the sum equals the explicit formula to ~10⁻¹⁰ (§3); the prime bumps are **negative**, as the formula says. |
| A2 | PRIMEGAME step counts disagree with the video (281, 710, 2375 against my 280, 707, 2363) | **COUNTING CONVENTION** | the outputs and their order are identical (first 20 terms match the published sequence). The video's figures equal my steps plus the number of times 15/14 has fired — 1, 3, 12 (and 17, 34 for 2¹³, 2¹⁷) — exactly. Their counter evidently charges 15/14 twice; I could not confirm why. |
| A3 | The video says the flags are "five" registers | **VIDEO IMPRECISION** | six registers (11, 13, 17, 19, 23, 29) reach 1 and only ever take the values 0, 1 (§7). |
| A4 | 393 Gram intervals with ≠ 1 zero against 197 failing Gram points | **EXPECTED** | each failing point spoils the two intervals that meet at it (§2 counts the adjacency) |
| A5 | zeros' spacing variance (0.15) and Σ²(5) (0.34) are below GUE (0.18, 0.51) | **FINITE-HEIGHT EFFECT (expected)** | the higher half is only slightly closer to GUE (0.147→0.152, 0.321→0.360, §9 cell above), consistent with a slow approach, with corrections of order 1/ln²(E) (Bogomolny–Keating), and Σ² saturates at L ≳ ln(E/2π) spacings |
| A6 | a near-double zero at n = 1496, 1497 (t = 1977.1739, 1977.2714): gap 0.089 spacings, |Z| extremum 0.022 (typical ~1.4) | **LEGAL, GENUINE** (after a **METHOD** error) | my first golden-section search used a ±0.2 window around the grid peak and ran into the neighbouring lobe (|Z| = 0.315); restricted to the gap between the two mpmath zeros it gives t* = 1977.2224, |Z| = 0.0221, agreeing with the Riemann–Siegel grid (0.0219) |
| A7 | 664 off-line sign changes of the balance M in [20,400] | **ARTEFACT OF THE TRUNCATED HALVES; NOT ZEROS; association with zeros UNRESOLVED** | they lie on smooth branches (§6) with |ζ| not small at their midpoints (0.4–1.5); their distance to the nearest zero, in local spacings, is slightly **smaller** than for random t (§6 prints both): a real association in the truncated model whose mechanism is not identified (**UNRESOLVED**, next test: follow a branch to its end and see whether it terminates on a zero). The crude bound on c has zero slack at t = 2πN²: the reason the margin is *marginal*, not a violation. |
| A8 | S(γₖ) at the zeros reaches ±0.81 but S(t) between zeros only ±0.54 | **CONVENTION** | S jumps by 1 at each zero; at the zero it is the mean of the two sides |
| A9 | My earlier suggestion to map "the nodal line of the modulus balance" | **CORRECTION (MATHS)** | for the exact functions ζ(s) = χ(s)ζ(1−s) gives |ζ(s)| = |χ||ζ(1−s)| identically, so the exact balance is not a line; only the truncated halves' balance is informative (§6) |
| A10 | H_NN (the ValaQuenta ladder) has weyl 203, locate 137 | **EXPECTED FAILURE** | the control fails every score, as it should; it is a harmonic ladder, not a spectrum of the zeros |

### What this says about mapping the space
1. **One circle carries the density, the rest carry the information.** The n = 1 circle through the origin reproduces every zero to within one spacing (mean 0.22). The zeros' individuality is the shift S, bounded below 0.81 for 3066 zeros.
2. **The two sides of the duality are both visible and both exact.** Primes → zeros (notebook 02) and zeros → primes (§3, identity to 10⁻¹⁰): the zeros' trace has the primes as fixed, negative echoes at ln p, and the pole term cosh(t/2) is the path of least primes.
3. **No candidate operator has been found; a test exists.** H_NN and the xp-in-a-box ladder fail; GUE has the statistics without the positions; the carrier has the density without the statistics. A real candidate must pass all four columns of §4.
4. **The margin is marginal.** The balance of the two halves leaves exactly zero slack at the heights where wheels switch on; that is where a zero *could* leave the line in the truncated model, and they do not (A7). It is a statement about the approximation, not about ζ.
5. **Near-misses are where the structure is.** Margins of 0.016 (A6) are the places to watch for any bifurcation; none occurs below 3600 (all 3066 zeros are on the line, Platt–Trudgian rigorously to 3·10¹²).
6. **PRIMEGAME is the discrete dual of the same structure:** primes as the survivors of the divisor test; cost ≈ 2.2·p^2.9 steps against ≈ p ln p zeros needed to resolve a prime in the explicit formula.

| item | label |
|---|---|
| carrier circle; zero − carrier < 1 spacing; θ(γ) = π(k − 3/2 − S) | `:CALCULATED` |
| negative real-axis crossings = Gram-law failures; 197 of 3066 below 3600 | `:CALCULATED` |
| echo spectrum = explicit formula (10⁻¹⁰) | ESTABLISHED formula, `:CALCULATED` |
| harness scores and controls | `:CALCULATED`; **no operator candidate passes** |
| Lehmer margins; smallest 0.016 near 1977 | `:CALCULATED` |
| balance crossings; crude-bound slack zero at 2πN² | `:CALCULATED` / ESTABLISHED bound |
| PRIMEGAME generator, register roles, cost p^2.9 | ESTABLISHED program, `:CALCULATED` |
""")
nbf.write(nb, "03_clock_echo_harness_margins_primegame.ipynb")
