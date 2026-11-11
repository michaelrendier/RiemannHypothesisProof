"""Beurling systems and the non-UFD level (rhp.beurling).  Writes beurling_probe.json."""
import json, sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import beurling as B
from rhp.sedenion_slots import primes_upto

R = {}; X = 100_000_000; t0 = time.time()
P = primes_upto(X); print("primes", len(P), f"{time.time()-t0:.0f}s", flush=True)
xs = np.exp(np.linspace(np.log(1e4), np.log(X), 600))

# 1  true primes: RMS of (ψ − x)/√x, and the prediction from the zeros
psiT = B.psi_beurling(P, np.ones(len(P), np.int64), xs)
g = np.array([float(l.split()[-1]) for l in open("data/zeros_25000_hp.tsv") if l[0].isdigit()])
pred = float(np.sqrt(2 * np.sum(1 / (0.25 + g ** 2))))
wins = [(1e4, 1e5), (1e5, 1e6), (1e6, 1e7), (1e7, 1e8)]
R["true"] = {"rms_windows": [B.rms_over(psiT, xs, *w) for w in wins], "zero_sum_prediction": pred}
print("true primes RMS (ψ−x)/√x per decade:", [round(v, 3) for v in R["true"]["rms_windows"]], " predicted from zeros", round(pred, 3), flush=True)

# 2  gated primes: PNT fails
R["gated"] = {}
for tau in (1.0, 3.0):
    m = B.gated_multiplicity(P, tau); psiG = B.psi_beurling(P, m, xs); ratio = psiG / xs; pr = B.gated_prediction(xs, tau)
    sel = xs >= 1e5
    R["gated"][tau] = dict(ratio_min=float(ratio[sel].min()), ratio_max=float(ratio[sel].max()), pred_min=float(pr[sel].min()), pred_max=float(pr[sel].max()),
                           max_abs_err_vs_prediction=float(np.abs(ratio[sel] - pr[sel]).max()), true_ratio_range=[float((psiT / xs)[sel].min()), float((psiT / xs)[sel].max())])
    print(f"gated τ={tau}: ψ_P/x ∈ [{ratio[sel].min():.3f}, {ratio[sel].max():.3f}] (prediction [{pr[sel].min():.3f}, {pr[sel].max():.3f}], max |err| {np.abs(ratio[sel]-pr[sel]).max():.4f});  true ψ/x ∈ [{(psiT/xs)[sel].min():.4f}, {(psiT/xs)[sel].max():.4f}]", flush=True)
# generalized integers N(x)/x of the gated system, τ = 1
XN = 2_000_000; m = B.gated_multiplicity(P[P <= XN], 1.0); a = B.generalized_integers(P[P <= XN], m, XN); Nc = np.cumsum(a)
grid = [10 ** e for e in (3, 4, 5, 6)] + [2_000_000]; R["gated_N_over_x"] = {int(x): float(Nc[x] / x) for x in grid}
lg = np.exp(np.linspace(np.log(1e3), np.log(XN), 400)).astype(int); r = Nc[lg] / lg
R["gated_N_osc_amplitude_by_decade"] = {f"{lo:.0e}-{hi:.0e}": float((r[(lg >= lo) & (lg <= hi)].max() - r[(lg >= lo) & (lg <= hi)].min()) / 2) for lo, hi in [(1e3, 1e4), (1e4, 1e5), (1e5, 1e6), (1e6, 2e6)]}
print("gated N(x)/x:", {k: round(v, 4) for k, v in R["gated_N_over_x"].items()}, " half-range per decade:", {k: round(v, 4) for k, v in R["gated_N_osc_amplitude_by_decade"].items()}, flush=True)

# 3  random (Cramér) primes
R["random"] = {"seeds": []}
for seed in range(4):
    Q = B.random_system(X, seed); psiR = B.psi_beurling(Q, np.ones(len(Q), np.int64), xs); R["random"]["seeds"].append([B.rms_over(psiR, xs, *w) for w in wins])
R["random"]["predicted_sqrt_ln_x_minus_2"] = [float(np.sqrt(np.log(np.sqrt(lo * hi)) - 2)) for lo, hi in wins]
print("random RMS per decade (4 seeds):", [[round(v, 2) for v in s] for s in R["random"]["seeds"]], " √(ln x − 2):", [round(v, 2) for v in R["random"]["predicted_sqrt_ln_x_minus_2"]], flush=True)

# 4  Z[√−5]
R["z5"] = {"brute_vs_formula": {}}
for Xb in (300, 1000):
    R["z5"]["brute_vs_formula"][Xb] = [B.z5_irreducibles_brute(Xb), B.z5_irreducibles(Xb)[2]]
print("Z[√−5] brute vs formula (irreducibles, norm ≤ X):", R["z5"]["brute_vs_formula"], flush=True)
R["z5"]["table"] = {}
for e in range(3, 9):
    Xz = 10 ** e; n_all, n_pr, n_irr = B.z5_irreducibles(Xz); L = np.log(Xz)
    R["z5"]["table"][Xz] = dict(prime_ideals=n_all, principal=n_pr, irreducibles=n_irr, ideals_x_lnx_over_x=n_all * L / Xz, irr_x_lnx_over_x=n_irr * L / Xz)
    print(f"  X=1e{e}: prime ideals {n_all}  principal {n_pr}  irreducible elements {n_irr}   π_ideal·lnX/X = {n_all*L/Xz:.3f}   π_irr·lnX/X = {n_irr*L/Xz:.3f}   (ln ln X = {np.log(L):.3f})", flush=True)
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
print(f"done {time.time()-t0:.0f}s")
