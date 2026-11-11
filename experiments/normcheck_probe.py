"""Norm check at the interfaces of the CD tower (rhp.normcheck).  Writes normcheck_probe.json."""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import normcheck as nc

res = {"levels": {}, "basis_pair_sums": {}}
for N in (1, 2, 4, 8, 16, 32, 64, 128, 256):
    r = nc.level_report(N, m=1000 if N >= 128 else 2000); res["levels"][N] = r
    row = f"T_{N:<4} all max|D|={r['all']['max_abs']:.3e}  mean={r['all']['mean_abs']:.3e}  D∈[{r['all']['min']:+.2f},{r['all']['max']:+.2f}]"
    if N >= 2: row += "   " + " ".join(f"{k}:{r[k]['max_abs']:.1e}" for k in ("LL", "LU", "UL", "UU"))
    print(row)
for N in (8, 16, 32):
    v = nc.basis_pair_sums(N); res["basis_pair_sums"][N] = {str(k): c for k, c in sorted(v.items())}
    tot = sum(v.values()); print(f"T_{N}: |(e_a±e_b)(e_c±e_d)|² values", dict(sorted(v.items())), f"  zero share {v.get(0,0)/tot:.4f}")
json.dump(res, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
