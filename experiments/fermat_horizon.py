"""The Fermat side of the horizon: exceptional primes of xⁿ + yⁿ = zⁿ mod p, certified complete up to the Weil threshold (rhp.fermat_curve).  Writes fermat_horizon.json."""
import json, sys, pathlib, math, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from rhp import fermat_curve as F
t0 = time.time(); R = {}
for n in range(2, 21):
    thr = F.weil_threshold(n); pmax = int(math.ceil(thr)) + 1; ex = F.exceptional_primes(n, pmax)
    big = max(ex) if ex else None
    R[n] = dict(genus=F.genus(n), tones=(n - 1) * (n - 2), weil_threshold=thr, certified_to=pmax, exceptional=ex, count=len(ex), largest=big, ratio_largest_to_threshold=(big / thr if big else None))
    print(f"n={n:2d} g={F.genus(n):3d} 2g={(n-1)*(n-2):3d} Weil threshold {thr:10.1f}  exceptional primes ({len(ex)}): {ex if len(ex) <= 12 else ex[:6] + ['…'] + ex[-4:]}   largest/threshold {big / thr if big else float('nan'):.3f}   [{time.time()-t0:.0f}s]", flush=True)
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
