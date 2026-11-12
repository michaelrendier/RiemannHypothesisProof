"""Which pairs of curves built from Cody's 'natural terms' (π, e, ln, 2/π, 2/e, 2/ln, α, golden) cross at d* (and at Ω)?  And is that better than chance?
Grammar (from the 2026-04-11 Gemini chat): basis functions of x times a natural multiplier.  Crossings found on x ∈ (0.02, 3) by sign change + brentq; matches counted to 1e-6 relative.  Writes natural_crossings.json."""
import json, math, sys, pathlib, itertools
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
OM = 0.5671432904097838; DS = OM / math.log(10); PHI = (1 + 5 ** .5) / 2; AL = 1 / 137.035999084
F = {"x": lambda x: x, "x^2": lambda x: x * x, "sqrt x": np.sqrt, "ln x": np.log, "ln(1+x)": np.log1p, "1/x": lambda x: 1 / x, "e^x": np.exp, "e^-x": lambda x: np.exp(-x),
     "10^-x": lambda x: 10.0 ** (-x), "1-e^-x": lambda x: 1 - np.exp(-x), "x e^x": lambda x: x * np.exp(x), "ln x /x": lambda x: np.log(x) / x}
A = {"1": 1, "2/pi": 2 / math.pi, "2/e": 2 / math.e, "pi/2": math.pi / 2, "e/2": math.e / 2, "1/pi": 1 / math.pi, "1/e": 1 / math.e, "pi": math.pi, "e": math.e, "ln10": math.log(10),
     "1/ln10": 1 / math.log(10), "alpha": AL, "1/alpha": 1 / AL, "phi": PHI, "1/phi": 1 / PHI, "2pi": 2 * math.pi, "2": 2.0}
curves = [(f"{an}*{fn}", (lambda a, f: (lambda x: a * f(x)))(av, ff)) for (fn, ff) in F.items() for (an, av) in A.items()]
xs = np.linspace(0.02, 3.0, 6001); cross = []                      # (x*, nameA, nameB)
vals = [np.where(np.isfinite(c(xs)), c(xs), np.nan) for _, c in curves]
with np.errstate(all="ignore"):
    for i, j in itertools.combinations(range(len(curves)), 2):
        d = vals[i] - vals[j]; idx = np.nonzero(d[:-1] * d[1:] < 0)[0]
        for k in idx:
            try: cross.append((brentq(lambda x: curves[i][1](x) - curves[j][1](x), xs[k], xs[k + 1], xtol=1e-14), curves[i][0], curves[j][0]))
            except Exception: pass
xsol = np.array([c[0] for c in cross]); print("curves", len(curves), "pairs", len(curves) * (len(curves) - 1) // 2, "crossings", len(cross), flush=True)
def hits(T, tol=1e-6): return [c for c in cross if abs(c[0] - T) <= tol * T]
rng = np.random.default_rng(0); T_rand = np.exp(rng.uniform(math.log(0.05), math.log(2.5), 3000))
null = np.array([len(hits(t)) for t in T_rand]); R = {"n_curves": len(curves), "n_crossings": len(cross), "null_mean_hits": float(null.mean()), "null_p99": float(np.percentile(null, 99)), "targets": {}}
for name, T in (("d* = Ω/ln10", DS), ("Ω", OM), ("0.24600", 0.246), ("1/phi^2", PHI ** -2), ("1/e", 1 / math.e), ("2/pi", 2 / math.pi)):
    h = hits(T); R["targets"][name] = dict(T=T, n_hits=len(h), examples=[f"{a}  =  {b}   (x* = {x:.12f})" for x, a, b in h])
    print(f"{name:12s} T={T:.9f}  hits within 1e-6: {len(h):4d}  (null mean {null.mean():.2f}, p99 {np.percentile(null,99):.0f})", flush=True)
    for x, a, b in h[:6]: print("     ", a, "=", b, f"  x*={x:.12f}")
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
