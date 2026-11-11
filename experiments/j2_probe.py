"""The J₂ definition tested (rhp.j2).  Writes j2_probe.json."""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import j2, sedenion_slots as SS

R = {}
# A — operators
R["operators"] = {}
for N in (8, 16, 32, 64):
    Jl, Jp = j2.J_left(N), j2.J_perm(N)
    r = {"L4": dict(order=j2.order(Jl), proj_order=j2.projective_order(Jl), swaps=j2.swaps_colour(Jl, N), square_is_minus_I=bool(np.allclose(Jl @ Jl, -np.eye(N))), auto_defect=j2.automorphism_defects(Jl, N)),
         "P4": dict(order=j2.order(Jp), proj_order=j2.projective_order(Jp), swaps=j2.swaps_colour(Jp, N), auto_defect=j2.automorphism_defects(Jp, N)),
         "P4_moves_e0_to": int(np.nonzero(Jp[:, 0])[0][0]),
         "all_Lk_square_minus_I": all(np.allclose(SS.Alg(N).Lmat(k) @ SS.Alg(N).Lmat(k), -np.eye(N)) for k in range(1, N))}
    R["operators"][N] = r
    print(f"T_{N}: L4 order {r['L4']['order']} (rays {r['L4']['proj_order']}) swaps {r['L4']['swaps']} L4²=−I {r['L4']['square_is_minus_I']} auto/anti defect {r['L4']['auto_defect'][0]:.2f}/{r['L4']['auto_defect'][1]:.2f}"
          f" | P4 order {r['P4']['order']} swaps {r['P4']['swaps']} auto/anti {r['P4']['auto_defect'][0]:.2f}/{r['P4']['auto_defect'][1]:.2f} | all L_k²=−I {r['all_Lk_square_minus_I']}")
R["stroke_cycle_e0"] = j2.orbit_of_identity(16); print("J₂^n e₀:", R["stroke_cycle_e0"])

# B — Klein group of the functional equation
rng = np.random.default_rng(1); pts = [complex(rng.uniform(-0.5, 1.5), rng.uniform(1, 40)) for _ in range(12)]
R["klein"] = j2.klein_check(pts); print("Klein:", R["klein"])
g = np.loadtxt("data/zeros_25000_hp.tsv", usecols=-1, dtype=float) if False else None
import csv
zs = []
for line in open("data/zeros_25000_hp.tsv"):
    t = line.split()
    try: zs.append(float(t[-1]))
    except ValueError: pass
orb = [j2.klein_orbit_size(complex(0.5, z)) for z in zs]; R["klein_orbits"] = {"zeros": len(zs), "size2": int(sum(o == 2 for o in orb)), "size4": int(sum(o == 4 for o in orb))}
print("Klein orbits of the zero table:", R["klein_orbits"])

# C — arithmetic J₂
P = SS.primes_upto(1_000_000)[1:]
inert, none = j2.hole_shares(P[:20000]); R["hole_share_inert"] = {str(d): v for d, v in inert.items()}; R["share_trivial_frobenius"] = {str(m): v for m, v in none.items()}
print("inert (hole) share per level:", {d: round(v, 3) for d, v in inert.items()}); print("trivial-Frobenius share (m = 1…7):", {m: round(v, 4) for m, v in none.items()})
R["decomposition_check"] = {str(m): j2.decomposition_check(P[:20000], SS._DISC[:m]) for m in (1, 2, 3, 4)}; print("ζ_K decomposition check (checked, failures):", R["decomposition_check"])

# D — the slot rules against the J₂ axis
pr = SS.primes_upto(40000)[:2000]; R["rule_shares"] = {}
for N in (16, 32, 64, 128):
    R["rule_shares"][N] = {}
    tab = j2.relation_table(N); R["rule_shares"][N]["slot_classes"] = {c: tab.count(c) for c in ("zero", "axis", "anti", "comm", "deform")}
    print(f"T_{N} slot classes vs e_4:", R["rule_shares"][N]["slot_classes"])
    for rule in ("seq1", "mod", "char", "charx", "random"):
        sh = j2.rule_class_shares(N, rule, pr, rng=np.random.default_rng(3)); R["rule_shares"][N][rule] = sh
        print(f"   {rule:6s}", {k: round(v, 3) for k, v in sh.items()})
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
