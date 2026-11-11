# SPDX-License-Identifier: GPL-3.0-only
"""zdtree_probe.py — the Zero-Divisor Tree T_16 … T_256 rebuilt in RHP (rhp/zdtree.py) and checked: counts, descent rule, clean corridors, survivors, wiki-88 closed forms, equality with the anticommutation graph.   :CALCULATED
Writes experiments/zdtree_probe.json."""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
from rhp import zdtree as T, sedenion_slots as SS
out = {}; prev = None
for k in (4, 5, 6, 7, 8):
    t0 = time.time(); n = 2 ** k; Z, per, tot = T.zd_incidence(k); carry, clean = T.classes(per, tot); iu = np.triu_indices(n, 1)
    zd = int(Z[iu].sum()); pairs = (n - 1) * (n - 2) // 2
    row = dict(k=k, dim=n, zd_pairs=zd, imaginary_pairs=pairs, fraction=zd / pairs, survivors=pairs - zd, carrying_delta=len(carry), clean_delta=clean,
               closed_form_carrying=2 ** k - k - 5, wiki88_zd_ordered=n * (n + 1.5 - 3 * k) - 4, ordered_from_count=2 * zd)
    row["clean_rule_holds"] = (clean == list(range(1, 8)) + [2 ** j for j in range(3, k)])
    if prev is not None:
        h = n // 2; old = Z[:h, :h]
        row.update(new_zd=zd - prev["zd_pairs"], old_pairs_status_unchanged=bool((old[1:, 1:] == prev["Z"][1:, 1:]).all()),
                   new_pairs=pairs - prev["imaginary_pairs"], new_fraction=(zd - prev["zd_pairs"]) / (pairs - prev["imaginary_pairs"]),
                   old_zd_classes_stay=set(prev["carry"]) <= set(carry), old_clean_stay=set(prev["clean"]) <= set(clean), new_clean=sorted(set(clean) - set(prev["clean"])),
                   new_carrying=len(set(carry) - set(prev["carry"])), survivor_recursion_ok=(pairs - zd) == 2 * prev["survivors"] + 3 * (2 ** (k - 1) - 1))
    if k <= 7:
        D = SS.deform_matrix(n); row["equals_anticommutation_graph"] = bool((Z[1:, 1:] == D[1:, 1:]).all())
    out[str(k)] = {kk: vv for kk, vv in row.items()}
    prev = dict(row, Z=Z, carry=carry, clean=clean)
    print(k, zd, row.get("new_zd"), len(carry), round(time.time() - t0, 1), "s"); sys.stdout.flush()
json.dump(out, open(os.path.join(HERE, "zdtree_probe.json"), "w"), indent=1)
