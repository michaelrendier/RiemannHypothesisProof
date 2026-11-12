"""Curves built from the spacing of the primes (the Fermat domain: composites between primes) in the digit coordinate D (x = 10^D): which pairs cross at d*?  Writes gap_curves_scan.json."""
import json, math, sys, pathlib, itertools
import numpy as np
from scipy.optimize import brentq
OM = 0.5671432904097838; LN10 = math.log(10); DS = OM / LN10; LOGE = 1 / LN10; PHI = (1 + 5 ** .5) / 2; AL = 1 / 137.035999084
g = lambda D: LN10 * D                                           # mean prime gap at x = 10^D  (ln x)
F = {"D (digits of x)": lambda D: D, "g = ln x (mean gap)": g, "1/g (prime density)": lambda D: 1 / g(D), "1-1/g (composite density)": lambda D: 1 - 1 / g(D),
     "g-1 (composites per gap)": lambda D: g(D) - 1, "D(g-1) (digit weight / gap)": lambda D: D * (g(D) - 1), "(D+1)(g-1)": lambda D: (D + 1) * (g(D) - 1),
     "log10 g (digits of the gap)": lambda D: np.log10(g(D)), "(1+ln g)/ln10 (gap entropy, digits)": lambda D: (1 + np.log(g(D))) / LN10,
     "x = 10^D": lambda D: 10.0 ** D, "1/x = 10^-D": lambda D: 10.0 ** (-D), "D^2": lambda D: D * D, "1/D": lambda D: 1 / D, "e^-D": lambda D: np.exp(-D), "g^2": lambda D: g(D) ** 2}
A = {"1": 1, "2/pi": 2 / math.pi, "2/e": 2 / math.e, "pi/2": math.pi / 2, "e/2": math.e / 2, "1/pi": 1 / math.pi, "1/e": 1 / math.e, "pi": math.pi, "e": math.e, "ln10": LN10, "1/ln10": LOGE,
     "alpha": AL, "1/alpha": 1 / AL, "phi": PHI, "1/phi": 1 / PHI, "2pi": 2 * math.pi, "2": 2.0}
curves = [(f"{an}*[{fn}]", (lambda a, f: (lambda D: a * f(D)))(av, ff)) for fn, ff in F.items() for an, av in A.items()]
Ds = np.linspace(0.02, 3.0, 6001); vals = []
with np.errstate(all="ignore"):
    for _, c in curves: v = c(Ds); vals.append(np.where(np.isfinite(v), v, np.nan))
    cross = []
    for i, j in itertools.combinations(range(len(curves)), 2):
        d = vals[i] - vals[j]; idx = np.nonzero(d[:-1] * d[1:] < 0)[0]
        for k in idx:
            try: cross.append((brentq(lambda D: curves[i][1](D) - curves[j][1](D), Ds[k], Ds[k + 1], xtol=1e-14), curves[i][0], curves[j][0]))
            except Exception: pass
print("curves", len(curves), "pairs", len(curves) * (len(curves) - 1) // 2, "crossings", len(cross), flush=True)
def hits(T, tol=1e-6): return [c for c in cross if abs(c[0] - T) <= tol * T]
rng = np.random.default_rng(0); T_rand = np.exp(rng.uniform(math.log(0.05), math.log(2.5), 3000)); null = np.array([len(hits(t)) for t in T_rand])
R = {"n_curves": len(curves), "n_crossings": len(cross), "null_mean": float(null.mean()), "null_p99": float(np.percentile(null, 99)), "targets": {}}
import re
def parse(s): m = re.match(r"(\S+?)\*\[(.+)\]$", s); return m.group(1), m.group(2)
for name, T in (("d* = Ω/ln10", DS), ("log10 e = 1/ln10", LOGE), ("0.24600", 0.246), ("Ω", OM), ("1/ln10 * Ω? (0.2463 x 1)", DS)):
    h = hits(T); rows = []
    for x, a, b in h:
        (ma, fa), (mb, fb) = parse(a), parse(b); rows.append((a, b, ma == mb))
    R["targets"][name] = dict(T=T, n_hits=len(h), n_equal_multiplier=sum(1 for r in rows if r[2]), distinct_curve_pairs=sorted({tuple(sorted((parse(a)[1], parse(b)[1]))) for a, b, _ in rows}), different_multiplier=[f"{a} = {b}" for a, b, e in rows if not e][:20])
    print(f"{name:26s} T={T:.9f} hits {len(h):4d} (equal multiplier {sum(1 for r in rows if r[2])}); null mean {null.mean():.2f}")
    for pr in R["targets"][name]["distinct_curve_pairs"][:14]: print("      ", pr)
    if R["targets"][name]["different_multiplier"]: print("   different multipliers:", R["targets"][name]["different_multiplier"][:6])
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
