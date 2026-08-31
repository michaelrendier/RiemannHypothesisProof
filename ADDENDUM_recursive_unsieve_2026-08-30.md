# Addendum — The Recursive Un-Sieve: birth order vs extinction order

**Companion to `ADDENDUM_generational_lineage_2026-08-28.md` §D.** That
addendum measured the sieve as *extinction* — each composite falls on the
pass of its smallest prime factor, `generation(n) = π(spf(n))`. This one
measures the dual: the **un-sieve**, from the ground state *"Just Prime
Numbers"* — turn primes on one at a time and watch each composite **arrive**,
born the moment its last needed prime factor is switched on.

**Status.** Everything below is a measurement, `N = 10⁵`, engine
`FactoralDecomposition/engine/lineage.py::un_sieve` (script and full output:
`ContextPlease/claude/scratchpad/2026-08-30_prime-dna/un_sieve.py`). **Not a
proof.** Does not close C1. It adds a construction-side residual and a
no-new-maths route to a Zeta Hamiltonian.

---

## A. The four orders

`rank_asc(p) = π(p) − 1` (pass 0 = the p=2 wave). `rank_desc(p)` runs the
prime list backwards.

| | reading | generation | gen. range | entropy | top pass |
|---|---|---|---|---|---|
| **A** | extinction, low→high | `rank_asc(spf(n))` | `[0 … 64]` | **2.491 bits** | 55.3 % |
| **B** | extinction, high→low | `rank_desc(gpf(n))` | `[4459 … 9591]` | **9.685 bits** | 0.0 % |
| **C** | birth, low→high | `rank_asc(gpf(n))` | `[0 … 5132]` | **9.685 bits** | 0.0 % |
| **D** | birth, high→low | `rank_desc(spf(n))` | `[9527 … 9591]` | **2.491 bits** | 0.0 % |

- **A** is the classic sieve: a composite dies on the pass of its *smallest*
  prime; 55 % of all composites die on pass 0 (multiples of 2).
- **D** turns primes on from the top down; a composite is born when its
  *smallest* prime is finally reached. **`D == reverse(A)` exactly** — bit
  for bit, over all 90 407 composites. Birth high→low **mirrors** extinction
  low→high, and the reflection preserves entropy: `H(A) = H(D) = 2.491`.
- **C** turns primes on from the bottom up; a composite is born when its
  *largest* prime (`gpf`) is reached. **`B == reverse(C)`**, `H = 9.685`.

**The gpf orders carry ≈ 4× the entropy of the spf orders.**
`H(C) − H(A) = +7.19 bits`. Same composites, same information — the birth
order spreads it across ~5 000 generations where the death order compresses
it into 65.

---

## B. The residual — C against A

The birth-by-`gpf` histogram minus the death-by-`spf` histogram, first 80
generations:

```
  gen   0  (prime    2)  : −49 984
  gen   1  (prime    3)  : −16 583
  gen   2  (prime    5)  :  −6 455
  gen   3  (prime    7)  :  −3 428
  gen   4  (prime   11)  :  −1 575
  gen   5  (prime   13)  :    −947
  gen  14  (prime   47)  :    +593
  gen  17  (prime   61)  :    +591
```

A large **negative front** at the small primes — almost every composite has
`spf = 2` but few have `gpf = 2`, so "died at 2" hugely outnumbers "born at
2" — decaying through `p = 13`, then small **positive ripples** at the
mid-sized primes (`p = 47, 61, …`) where `gpf` piles up. Negative spike,
positive fine structure: the shape §D.2 of the prior addendum predicted for
the departure of the firing order from ordinal, now measured from the
construction side.

### B.1 Clocked by zeta — the penalty is invariant, the schedule is not

The prior addendum §D flagged "ordinal vs ζ firing order" as open. Measured:
`un_sieve_zeta.py` re-runs all four reads with the ordinal prime rank
replaced by five zeta-derived orders — `zeta_weight` (`ln p/√p`, the σ=½
von-Mangoldt amplitude, peaks at `p = 7`), `theta` (`θ(2π p²)`, the
Riemann–Siegel theta at the height where `p` enters the RS main sum),
`Zsign` (the sign of the RS `Z`-function at that height — one bit, from ζ),
`spiral` (`θ mod 2π` + that bit), and the ordinal control.

**Invariant to five decimals across every ordering** (`N = 10⁵`):
`H(A) = 2.49101`, `H(C) = 9.68456`, **`H(C) − H(A) = +7.19355`**, C−A
residual |mass| = 158 958, `D == reverse(A)` true, the 60.5 % born-after
fraction, and the boundary primes 313 / 49 999. Entropy of the C (or A)
histogram depends only on the multiset of per-prime counts; a permutation of
the primes just relabels the bins. **Zeta order is not a shortcut through the
existence cost — the `+7.19` bits is a combinatorial invariant of ℕ.**

What zeta order *does* move is the **payment schedule**. `theta` and
`zeta_weight` are near-monotone in `p`, so they keep the compact ordinal
front (A-span 65, C-span 5133). But `Zsign` blows the spans to 781 / 1132 and
`spiral` to 1278 / 1437 (`N = 1.2×10⁴`): the same entropy, the same residual
mass, smeared across 25–50× more generations. **The compactness of the
extinction order — 55 % of composites dying on pass 0 — is an artefact of the
ordinal clock, not of ℕ.** The −49 984 residual spike sits at generation 0
under ordinal, generation 21 under `zeta_weight` (`p = 2` has low `ln p/√p`),
generation 766 under `Zsign`. Mass fixed, location free.

**The Go connection (Cody).** 313 is `π⁻¹(65)` — the 65th prime, the largest
with `p² ≤ N`, and under any *monotone* clock it is still activated 65th, with
zero wasted passes: 65 of 9 592 primes strike, and you stop. Under `spiral`,
313's rank is effectively random in `[0, 9592]`, so you wade through up to
**9 527 primes that strike nothing** (317, 331, …, every prime in `(√N, N]`)
before the boundary prime is reached — the exact shape of enumerating illegal
Go moves. **The ordinal order is the unique zero-rejection order**; that is
why it is the canonical minimum-entropy decomposition. An oscillatory ζ-clock
adds rejection work and buys nothing.

---

## C. The two boundary primes — a mass-gap-like separation

| process | completes when | prime | index |
|---|---|---|---|
| **extinction** (A) | `p² ≤ N` (last prime that strikes an as-yet-unstruck number) | **313** | π = 65 — the "313 Sieve" |
| **birth** (C) | `2p ≤ N` (last prime that is the `gpf` of a new composite) | **49 999** | π = 5133 |

**The sieve finishes killing at `√N` (313) but does not finish birthing until
`N/2` (49 999).** Between those two primes the sieve does no extinction work
at all, yet **60.5 % of every composite ≤ N is born after the extinction
boundary** (`gpf > 313`). More than half of ℕ's multiplicative content is
decided by primes that kill nothing.

That factor of `√N` vs `N/2` — a scale separation between the compact
low-entropy process and the broad high-entropy one describing the *same set*
— is the construction-side analogue of a mass gap: the two descriptions do
not overlap in their active range.

### C.1 The split — extinction is free, existence is not

The asymmetry is a cost asymmetry, and it is intrinsic to ℕ:

| | cost (primes needed) | entropy | when done |
|---|---|---|---|
| **factoring to extinction** (A) | primes ≤ `√N` — **313** | 2.491 bits | pass 65; 55% on pass 0 |
| **factoring to existence** (C) | primes ≤ `N/2` — **49 999** | 9.685 bits | pass 5133; flat |

Defining a number by **what it cannot be** — the sieve strikes it out the
moment *any* prime divides it, so its smallest suffices — is bounded, cheap,
and finished at `√N`. Defining it by **what it is** — every prime factor must
be present, so you must wait for its largest — costs `+7.19 bits` and is not
finished until `N/2`. `ζ` describes the primes by exclusion (through the
zeros); that description is the cheap side. The construction is the expensive
side, and the un-sieve is where the expense is paid. Same split as the
descriptive-vs-definitional pair of the 08-28 addendum §A, now with a price on
it.

---

## D. Zeta as the ground state, and a Zeta Hamiltonian for free

The un-sieve starts from **"Just Prime Numbers"** — the true ground state of
number theory — and the sequence of birth events (C) is a **path** from that
ground state to full ℕ, indexed by generation. That is an action / Lagrangian
picture: a trajectory, `dN_born/dg` its velocity.

`ζ(s)` **describes** this ground state — the Euler product runs over exactly
the primes that are the ground state, and nothing else. So `ζ` is a
**description of a path**: `ζ(½ + it)` traces a curve (a cardioid-family
curve — the caustic of a circle, the pedal curve; cf. `wiki/34` and the
Riemann–Siegel `ϑ(t)`).

**Mechanics are mechanics.** Any system with a Lagrangian / path description
has a Hamiltonian one by Legendre transform — swap the generation coordinate
for its conjugate, the birth rate. This is not a new claim and it needs no
new mathematics: the un-sieve makes the Lagrangian of the prime construction
explicit, and a **Zeta Hamiltonian follows by the standard L ↔ H
correspondence.** It is the Berry–Keating `H = xp` (`wiki/102`, `PAPER.md`
§4) reached from the construction side rather than the spectral side — a
third bearing on the same operator, alongside §3 (Noether currents) and
§B of the prior addendum (the harmonic Two-Trees field).

**Not claimed:** that any of this closes C1. It gives the construction-side
residual a shape, names the birth/extinction scale gap, and observes that the
existence of a Zeta Hamiltonian is licensed by ordinary mechanics once the
prime construction is written as a path.

### D.1 The zeros as the clock — a measured 15 % recovery

§B.1's rank orderings only *relabel* the generation axis, so the entropy gap
cannot move. A **time embedding** does let it move: give each prime a real
birth *time* `τ(n) = γ_{rank(gpf(n))}` — the height of the corresponding
Riemann zero — instead of an integer index, and bin `τ` uniformly. Now the
actual zero *spacing* enters, not just the order. Measured, `N = 8 000`,
real zeros `γ₁ = 14.13 … γ₁₀₀₇ = 1427.37`, 160 bins:

| clock | H(C) − H(A) |
|---|---|
| uniform (arrival at 1, 2, 3, …) | **+4.3007 b** |
| zero-time (arrival at `γ_k`) | **+3.6638 b** |

Clocked by the actual zeros the existence penalty **drops ≈ 0.64 bits
(≈ 15 % at this N)**. The sparse low zeros stretch the compact death front
(`H(A)` 0.48 → 1.75) while the birth spread compresses relatively. The C−A
residual autocorrelation goes lag-1 `+0.407` then flat `≈ −0.05` for lag ≥ 3
— one-step smoothing then white, the signature of level repulsion suppressing
the far residual.

The `+7.19` combinatorial floor of §B.1 still stands; what this shows is that
**the zeros are a better-matched clock for the construction than the integers
are** — replaying birth on the zero timeline recovers cost that ordinal
replay discards. A directional measurement (one `N`, one bin count), not a
scaling law, but it is the first quantitative sense in which *ζ is the tape*:
the trajectory that, replayed against, makes the backward (existence) pass
measurably cheaper. Engine: `un_sieve_zeta.py`,
`ContextPlease/claude/scratchpad/2026-08-30_prime-dna/un_sieve_zeta_RESULTS.md`.

---

## E. Engine and cross-references

- `FactoralDecomposition/engine/lineage.py` — `un_sieve(N)`, companion to
  `sieve_lineage`, `sieve_recurrence`, `two_trees`.
- `ContextPlease/claude/scratchpad/2026-08-30_prime-dna/` — `un_sieve.py`,
  `prime_dna.py`, outputs.
- prior addendum §D (ordinal vs ζ firing order), §B (harmonic Two-Trees
  field), the "313 Sieve".
- `ValaQuenta/wiki/un_sieve.md`, `Ainulindale/wiki/47_the_two_trees.md`
  (Telperion = extinction order, Laurelin = birth order),
  `AbrikosovTree/README.md`.

*Cody Michael Allison — 2026-08-30. Companion to `PAPER.md` v8.*
