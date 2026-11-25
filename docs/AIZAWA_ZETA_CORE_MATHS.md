# The core maths of the Aizawa pathway, and what it says about the zeta pathway

Written 2026-11-24 for Cody. Every number below was measured in this repository (experiments or tests named in each row); nothing here is evidence about the Riemann Hypothesis.
Code: `ValaQuenta/modules/aizawa_attractor`, `ValaQuenta/modules/zeta_uns_path`, `experiments/diagonal_quadratics.py`, viewer `visualizations/javascript/webgl_aizawa_vs_zeta_uns.html`.

## 1. What the Aizawa path is

The system dx/dt = (z−b)x − dy, dy/dt = dx + (z−b)y, dz/dt = c + az − z³/3 − (x²+y²)(1+ez) + f z x³ looks like a tangle. In cylindrical coordinates (r, φ, z) the first two equations are exact:

    dr/dt = (z − b) r        dφ/dt = d

The orbit **turns about the z axis at a constant rate d**, z sets the **growth rate of the radius**, and r can never reach the axis. Only the last term (f z x³) breaks the rotation symmetry; at f = 0 the 3-D orbit is the 2-D (r, z) flow swept round the axis (checked to 1.7×10⁻⁸). The fixed points sit on the axis with eigenvalues (z*−b) ± id and a − z*². The largest Lyapunov exponent is +0.085. A consequence worth keeping: ln r stays bounded on the attractor, so the long-time average of d(ln r)/dt = z − b is exactly zero, i.e. **⟨z⟩ = b** (measured 0.701 against b = 0.7 over 60,000 points; the orbit has z > b, radius growing, 48% of the time and z < b, radius shrinking, 52%). r and z are uncorrelated (−0.003) and mean z is the same in the tube, the middle and the shell, so the tube (a third of the time at r < 0.35) is not picked out by a value of z. The zeta analogue is the radial law ⟨Z′/Z⟩ ≈ 0 in the same bounded-log sense.

## 2. The same skeleton in zeta, and where it differs

Hardy: ζ(½+it) = e^(−iθ(t)) Z(t), Z real. So the path **turns at −θ′(t) = −½ ln(t/2π)** with a **signed radius Z(t)**.

| | Aizawa | ζ on the critical line | measured here |
|---|---|---|---|
| rotation rate | d, constant | −θ′(t), slowly rising | identity, residual 10⁻¹⁰ |
| radial law | (dr/dt)/r = z − b | (d/dt) ln\|Z\| = Z′/Z | Z′/Z = Σ over zeros 2t/(t²−γ²) + smooth: spread 0.573 → 0.007 away from the zeros (25,000 zeros); the remainder is a near-constant 0.75 |
| the axis | never reached | passed through at every zero (Z changes sign) | sign change at each of the first 30 zeros |
| driver | z, a third variable | the zeros themselves (poles of Z′/Z) | as above |

So "the same animation from the pole" is true of the skeleton and not of the picture: down the axis Aizawa is concentric spirals around a hole; ζ in the value plane is lobes that go through the centre, because its radius is signed. The spring is the UNS cylinder (ρ = log₁₀|ζ| + 3, φ = arg ζ, height = t) with the zeros on its axis.

**A lesson from the first render.** Sampling by equal steps of θ puts every point on 32 fixed angles of the flat point, which draws 32 rigid spokes. That was a sampling artefact. A golden-ratio jitter on the phase grid removes it (tested).

## 3. Where ½ comes from as a phase

The functional equation ζ(s) = χ(s) ζ(1−s) has |χ(½+it)| = 1 and χ = e^(−2iθ(t)): on the line it is a **pure phase**, and Z is what remains when that phase is taken out. σ does not rotate anything: in each term p^(−s) = p^(−σ) e^(−it ln p) the **t part is the phase** and the σ part is the **amplitude**; ½ is where the amplitude p^(−σ) of the series and p^(σ−1) of its mirror agree. The axis of rotation is therefore the **t (imaginary) direction**, not the σ axis; σ moves the radius.

## 4. Where the primes are

| view | what it shows | measured |
|---|---|---|
| prime spring at height t | one step per prime, length p^(−½), turned by −t ln p (loops uniform in ln p, spacing 2π/t, as the boundary engine's `boundary_uclock` finds) | steps equal p^(−s) to 0.02 for p > 70 |
| resonance heights | t_k = 2πk/ln p: prime p's phasor back at angle 0 | exact |
| zero spectrum Σ cos(γx) | spikes at x = ln p^k, the explicit formula on its zero side | 25,000 zeros: every prime to 47 and the powers 4, 8, 9, 25 spike 28× to 750× above the off-grid background |

## 5. The diagonals, as phase offsets, and the orthogonal wave

A diagonal of the Ulam spiral is a quadratic f(n). Shifting along it, n → n+δ, is a **phase offset in n**; it changes b and c but **keeps the discriminant** Δ = b² − 4c. The prime density depends on Δ through the **character (Δ/p)**:

| claim | result |
|---|---|
| the number of roots of f mod p is 1 + (Δ/p) (p odd, p ∤ Δ) | 2,267 of 2,267 cases |
| shifts of Euler's n²+n+41 (all Δ = −163) give the same density | counts 2,261 … 2,234 over 12 shifts: spread 9, against 47 for random noise (the slow drift is 1/ln f(n) as the window moves out) |
| Bateman–Horn predicts the counts across 120 random quadratics | r = 0.9994; observed/predicted median 0.996, 5–95% range 0.94–1.05 over the 89 with a nonzero prediction; the 31 with a fixed prime divisor are predicted 0 and produced 2 primes; Euler's 41: observed 2,261, predicted 2,266.8 |

So the "other wave orthogonal to the numbers" is the **Dirichlet character χ_Δ(p) = (Δ/p)**. It is orthogonal to the trivial character that gives ζ, it has its own L-function L(s, χ_Δ) with its own critical line, and a diagonal's prime density is L(1, χ_Δ) read off through Bateman–Horn. Polynomials with a large L(1, χ_Δ) (Δ = −163 for Euler's) are the dense diagonals. Bateman–Horn is a conjecture (1962); 0.9994 is agreement with it, not a proof.

## 6. The torus, and the gradient mesh

The CSS "Gradient Mesh" (#1074 in the collection) is nine radial gradients on a 3 × 3 grid of control points, summed with soft edges: a smooth field interpolated from control values (a radial-basis mixture; the vector-graphics version is a Coons or bicubic patch). The torus that does appear in this work is the **prime-phase torus**: the phases (t ln p₁, t ln p₂, …) mod 2π. For two primes the path (t ln 2, t ln 3) covers 100% of a 16 × 16, 32 × 32 and 64 × 64 grid of that torus (Kronecker: ln 2/ln 3 is irrational), and the two-prime Euler factors are a function on it, which is exactly a mesh. The UNS cylinder is around the t axis; it becomes a torus only when t is wrapped for finitely many primes. The strip-to-sphere chart (θ = πσ) is a sphere, not a torus.

## 7. The core tools, and where each helps

| tool | helps with | checked |
|---|---|---|
| rotation + radial driver (polar form) | reading any attractor or oscillating path; the ζ pathway | §1–2 |
| phasor sums and their Fourier dual | primes ↔ zeros (explicit formula); the spring | §4 |
| Kronecker–Weyl flow on a torus | several primes at once; the mesh picture | §6 |
| characters and Bateman–Horn | diagonals, Ulam spaces, "other waves" | §5 |
| A₂ root system (the trig hexagon) | all six trig functions from rotation; `rhp/trig_hexagon.py` | its tests |
| RK4 and Lyapunov exponents | any path you integrate | engine tests |

## 8. Open

- The Aizawa paper reference is not yet verified against a title page; the origin of the animation page is unknown.
- Asked, not yet done: `AbrikosovTree` Ulam spaces read as phase offsets of a shared Δ, and the L-function zeros of one χ_Δ drawn beside ζ's.
- The boundary engine's finding that different spirals end on the same object was read from its poster only, not analysed here.
