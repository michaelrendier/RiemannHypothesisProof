# Build plan — what follows from `paper/PAPER.md` (draft 0.1)

## Software to copy into `rhp/` (vendored, provenance header + source commit + GPL; `sync_check.py` flags drift)

| need (paper §) | source | status of the source |
|---|---|---|
| Noether currents, amplitude Lagrangian (§5) | `archive/.../engines/noether_derivation.py` | runs; imports `fermat_monster_engine` optionally |
| Wiles–Noether product e^(−E) (§5) | `FourthAgePapers/FermatMonster/engine/fermat_monster_engine.py` | exists; vendor only the functions used |
| π-factoring, equator, envelope, origin (§3–4) | `ContextPlease/claude/scratchpad/2026-08-13_apex_path/half_is_the_equator.py` | runs; T1–T4 |
| balance, σ expansion (§4, §10) | `ValaQuenta/modules/{noether,sigma_expansion}` | to check |
| forced σ, Newton form (§10) | `ValaQuenta/noether.py` (`forced_sigma`) | fixed 2026-08-28 |
| Berry–Keating xp candidate (§9) | `ValaQuenta/modules/berry_keating` (+ `NOTES_2026-08-31.md`) | to check |
| spin / wobble (§6, §8) | `ValaQuenta/modules/spectral_primes`, `GenerationalLineage/engine/toolsets/spectral_primes.py` | VQ part needs mpmath |
| Σ_RB candidate, ∅_RB (§5.2, §10) | `ValaQuenta/modules/h_rb_hat` | THEORETICAL |
| prime gauge field (§8) | `ValaQuenta/modules/prime_gauge_field` | one connection provably flat |
| Oblique Gear angles (§10.4) | `ValaQuenta/modules/oblique_gear` | to check |
| constants, four d\* (§10.5) | `ValaQuenta/modules/constants` | to check |
| unit circle catalogue (§4.3) | `GenerationalLineage/engine/toolsets/circle.py` | 29th toolset |
| spring coordinate u = ln x, Mellin (§7) | `GenerationalLineage/engine/toolsets/logperiodic.py` | numerical fit |
| tests, listings, doc generators | `ValaQuenta/tests/{test_known_results,test_listings,test_readme}.py`, `docs/gen_listings.py` | pattern to copy |

## To be WRITTEN (no code exists today)

1. `rhp/limacon.py` — the loop-fit sifting test (numbers are in the project record only; no longer on the critical path, §7 is rewritten).
2. `rhp/zero_counting.py` — §6 table from `mpmath`, with Gram's relation and spacing ratio.
3. `rhp/orthogonality.py` — the test of §8 (inner product, null model); currently only a conjecture.
4. ~~`rhp/lift.py` — the cusp-lift test~~ withdrawn: the pathway has no cusp. DONE instead (2026-11-03): `rhp/uns.py`, `rhp/flatten.py`, `rhp/render.py`, `rhp/viewport.py`/`camera.py` (vendored skill assets), `tests/test_uns.py` (11), notebook 01, figure + `data/pathway.csv` (Blender import). DONE 2026-11-04: `rhp/bigpath.py`, `flow.py`, `epicycles.py`, `render_big.py`, notebook 02, `tests/test_bigpath.py`. DONE 2026-11-04: `rhp/clock.py`, `gram.py`, `echo.py`, `harness.py`, `lehmer.py`, `balance.py`, `primegame.py`, notebook 03, `tests/test_notebook03.py` (29 tests total). NEXT: a real operator candidate for the harness (boundary conditions on xp), follow the balance branches to their ends (UNRESOLVED), sphere-chart picture, Blender render of `data/pathway.csv`.
5. `paper/make_paper.py` with `--check`.
6. `.venv` (numpy, mpmath, sympy, jupyter, nbconvert, pytest) and `requirements.txt`.

## Documents to create

`README.md` (informal introduction + contents, as GenerationalLineage), `docs/` code-reference Manual (API per module, Listings, transcripts),
`notebooks/` (one per paper section: 01 π-factoring · 02 emergent-½ catalogue · 03 amplitude Lagrangian · 04 zero-counting phase ·
05 cardioid sifting · 06 prime & factor wavelength · 07 the open step), `tests/`, `archive/` (done), `CONTRIBUTING.md`, `CITATION.cff`, `references/INDEX.md`.


## TODO — The Zipper (added 2026-11-16, Cody)

A mechanism, to be built and measured, that **zips and unzips the two domains along the positive integer line**: the primes (Telperion, the sieve) and the factor spacings (Laurelin, the un-sieve; Fermat's side), i.e. T + L = ζ − 1 read as two interlocking rows. Name: *the Zipper*.
Standing facts the Zipper starts from (all already in this repository, `:CALCULATED`): T + L = ζ − 1 exact (item 19); W_L = −W_T to ≤ 2/√x, antiphase corr −0.82 (items 16, 39); every integer n ≥ 2 sits in exactly one row (prime or composite); the Möbius inversion 1/ζ is the unzip; the wheels (frequencies ln p) are incommensurate and never lock exactly (item 22).
Engineering model to test (INTERPRETATION until measured): teeth = integers; two staggered rows; a registry condition (one pitch); a head/neck lock; a slider that enforces local order one pair at a time. See `paper/PAPER.md` §11.5 for the stub and the measurements it needs. Not built; no code yet.
