# Zero-phase windings against Archimedean (Ulam, polygonal) windings — 2026-11-11   :CALCULATED

Code `rhp/zerophase.py`, `rhp/atlas_riemann_fermat/zerophase.py` → `zerophase.svg`; run `experiments/zerophase_probe.py` → `.json`, `zerophase_scans.npz`; test `test_zerophase`.

## Notation (the code names the variation)
- **winding** `'zero'`: θ_k(n) = γ_k ln n. A zero-tone x^ρ = e^{(½+iγ)u} IS the logarithmic spiral r = √x, θ = γ ln x (identity), pitch arctan(1/2γ) = 2.03°, 1.36°, 1.15° for zeros 1–3.
  `'arch'`: θ_h(n) = 2πh√n. h = 1 is Sacks/Ulam (squares on one ray). The m-gon spiral is the SAME family at h = √(2/(m−2)) (polygonal numbers P_m(k) lie on the ray θ = 2πk at n = (m−2)k²/2): triangle 1.414, square 1, pentagon 0.816, hexagon 0.707. 'Which polygon' is one axis, h.
- **lift** `'stereo'`: inverse stereographic projection, equator at |z| = 1 (n = 1, or σ = ½ for the normalised z = x^{ρ−½}); conformal. `'linear'`: colatitude linear in ln n (what `AbrikosovTree/render/lattice_planes.py` actually draws; its comment block described the stereographic lift, never implemented). Now named in that file (`LIFT_KIND = 'linear'`, AbrikosovTree `1998b8d`).

## Results
1. **Loxodrome.** Under the stereographic lift each zero-tone has constant bearing arctan(2γ): 87.9740° (zero 1), 88.6374°, 88.8546°, 89.4241° (zero 10), 89.8772° (zero 100), spread ≤ 5e−4° and equal to theory to 1e−3°. Under the linear lift the bearing varies by ±1.41°, 0.95°, 0.80°, 0.40°, 0.09° respectively. The critical line is the equator; each zero is a rotation of rate γ_k about the polar axis, and the real standing tone is its shadow.
2. **One statistic for both windings:** A(θ) = |Σ_{n≤x}(Λ(n) − 1) n^{−½} e^{iθ(n)}| / ln x, x = 5·10⁴ (ln x = 10.8). Random-phase level √(Σw²)/ln x = 0.64, Rayleigh median 0.53.
3. **Zero scan, t ∈ [10, 100].** 29 zeros; 29 lines, mean height 0.919, minimum 0.834, largest offset from the ordinate 0.069 (< 1/ln x = 0.092); background: median 0.18, midpoint median 0.152, maximum between zeros 0.376. Contrast ≈ 6 to 1. (Landau–Gonek: the line height tends to ln x − O(1).)
4. **Archimedean scan, h ∈ [0.02, 8].** Median 0.425 (random-phase level), maximum 2.24 at h = 4.83 (the extreme of ~8000 noise samples); polygon values: 3-gon 0.27, 4-gon 0.21, 5-gon 1.12, 6-gon 0.15, 8-gon 0.57 at j = 1; the larger values at j·h are within the noise extremes. No polygon stands out.
5. **Caveat.** A is a GLOBAL coherent sum. It cannot see the local structure of the Ulam diagonals (quadratic families of length ~√x); those are a property of the square lattice and the local residue structure, not a global phase alignment of the primes. The statement is: no global alignment at any Archimedean winding, at every zero ordinate for the logarithmic one.

## Reading (not a result)
The arrangement the primes follow, as a global phase structure, is the zero-phase one: windings at constant rates γ_k in u = ln x, i.e. a flow on a torus with the incommensurate frequencies γ_k, whose line spectrum IS the zeros. The Ulam / polygonal spirals are projections whose rays are local quadratic families. The polygon the lattice picks (square ℤ[i], hexagon ℤ[ω]) is a question about discriminants (−4, −3), the same slot-bit data as the `char` rule, not about a global winding.
