# zeros25k — in-house recreation and re-run on the 25,000 exact zeros

Script: `../zeros25k_rerun.py`. Started 2026-11-09 14:22 as a disowned process (`setsid nohup nice -n 10 … & disown`, 4 workers).

| stage | what | output |
|---|---|---|
| A | every γ_k re-found at 32 digits (mpmath `siegelz` root); sign change across γ_k ± 1e-9; strictly increasing; table error | `../../data/zeros_25000_hp.tsv`, `results.json["A"]` |
| B | Backlund exact count N(T) = θ/π + 1 + S(T) at the mid-gap above every 50th zero and above the 25,000th (next zero from `mpmath.zetazero(25001)`) | `results.json["B"]`, `stageB.pkl` |
| C | wave reconstruction with K = 10…25,000 tones; prime-resolution crossover p* against T/2π; per-integer ping (ψ_K jump across [n−½, n+½] against Λ(n)) for K = 3066 and 25,000 | `stageC.json`, `jump_K*.npy` |

Check: `tail run.log`, `cat status.json`, `pgrep -f zeros25k_rerun`. Resume after a stop: rerun `python3 experiments/zeros25k_rerun.py A B C` (stages A and B checkpoint to `*.pkl`).
Labels: every number here is :CALCULATED (floating point or mpmath at the stated dps); nothing is evidence about RH.
