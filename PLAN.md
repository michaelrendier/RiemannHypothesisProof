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
4. ~~`rhp/lift.py` — the cusp-lift test~~ withdrawn: the pathway has no cusp. DONE instead (2026-11-03): `rhp/uns.py`, `rhp/flatten.py`, `rhp/render.py`, `rhp/viewport.py`/`camera.py` (vendored skill assets), `tests/test_uns.py` (11), notebook 01, figure + `data/pathway.csv` (Blender import). DONE 2026-11-04: `rhp/bigpath.py`, `flow.py`, `epicycles.py`, `render_big.py`, notebook 02, `tests/test_bigpath.py`. NEXT: sphere-chart picture, double-zero margin map (Lehmer pairs), Blender render of `data/pathway.csv`.
5. `paper/make_paper.py` with `--check`.
6. `.venv` (numpy, mpmath, sympy, jupyter, nbconvert, pytest) and `requirements.txt`.

## Documents to create

`README.md` (informal introduction + contents, as GenerationalLineage), `docs/` code-reference Manual (API per module, Listings, transcripts),
`notebooks/` (one per paper section: 01 π-factoring · 02 emergent-½ catalogue · 03 amplitude Lagrangian · 04 zero-counting phase ·
05 cardioid sifting · 06 prime & factor wavelength · 07 the open step), `tests/`, `archive/` (done), `CONTRIBUTING.md`, `CITATION.cff`, `references/INDEX.md`.
