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
