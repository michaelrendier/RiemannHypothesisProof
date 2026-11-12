"""Spacing of the primes by decade (rhp.gapdigits): mean gap, composites per gap, digit weight, gap entropy, composite density; digit-transition gaps at 10^m.  Primes to 2·10^8.  Writes gap_decades.json."""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from rhp import gapdigits as G
from rhp.sedenion_slots import primes_upto
P = primes_upto(200_000_000)
R = {"decades": [G.decade_stats(P, m) for m in range(1, 8)], "digit_transition_gaps": G.digit_transition_gaps(P, 8)}
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
for r in R["decades"]: print(r["m"], r["primes"], round(r["mean_gap"], 3), round(r["pnt_mean_gap"], 3), round(r["mean_composites"], 3), round(r["gap_entropy_nats"], 3))
