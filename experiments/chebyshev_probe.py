import sys, pathlib, json; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import chebyshev_bias as CB
g = CB.zeros(75.0); print("L(χ-4) zeros", len(g), g[:6].round(4), flush=True)
c = CB.c_chi(); print("c =", c)
X = 200000; D = CB.D_exact(X); R = {"zeros": len(g), "c": c, "first": g[:6].tolist()}
rng = np.random.default_rng(3)
for lo, hi in ((50, 200), (200, 2000), (2000, 20000), (20000, 200000)):
    x = np.floor(rng.uniform(lo, hi, 300)) + .5; ex = D[np.floor(x).astype(int)]
    for sgn in (1, -1):
        pass
    r = {K: float(np.sqrt(np.mean(((CB.D_zeros(x, g, K, c) - ex) / np.sqrt(x)) ** 2))) for K in (5, 20, len(g))}
    r["rms_exact/sqrtx"] = float(np.sqrt(np.mean((ex / np.sqrt(x)) ** 2))); R[f"{lo}-{hi}"] = r; print(lo, hi, {k: round(v, 3) for k, v in r.items()}, flush=True)
n = np.arange(3, X + 1); neg = float((D[3:] < 0).mean()); R["fraction_D_negative_to_2e5"] = neg; print("fraction of x≤2e5 with D<0 (inert primes lead):", round(neg, 4))
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
