# Ω = W(1) and d*, reproduced from the First Age experiment — 2026-11-11   :CALCULATED

Code `rhp/omega_dstar.py` (mpmath, 50 digits), test `test_omega_dstar`.
Sources read: `Ainulindale/AgeFirst/original_2026/ResearchPaper/SMNNIP_Omega_Constant_Paper.txt` and `SMNNIP_Maths_Notes_Log.txt` (2026-04-09), `conjecture_proof.py` (constants only, no Ω step); wikis 17, 50, 73 §4, 96; ValaQuenta `modules/constants/maths.py` (`derive_omega_zs`: iteration + Banach contraction; `derive_d_star`: the four faces, with d*_taut = Ω/ln 10 labelled 'tautological'); the old `AbrikosovTree/render/fermat_riemann_entropy.py` is a different entropy (zeros against gaps), not this experiment.

## The experiment as written
x = T/T_Planck. Forward entropy ceiling S_f/k_B = e^{−x} (falls); inverse entropy floor S_i/k_B = x (rises); in decimal digits D_f = −x/ln 10, D_i = log₁₀ x. They cross once, at x = e^{−x}.

## Reproduced (blind: the crossing is found by bisection on the two digit curves, with no W and no e^x solved)
- **Ω** = 0.567143290409783872999968662210355549753815787 by four routes: closed form W(1); the contraction x → e^{−x} (192 steps to 50 digits); Newton on x e^x − 1; bisection on D_f = D_i. Largest disagreement 3.5e−49.
- **d*** = −log₁₀ Ω = Ω/ln 10 = 0.246307201473422578…; closure d*·ln 10 = Ω exactly (the Translator); W branches W₀(0) = 0, W₀(1) = Ω, W₀(−1/e) = W₋₁(−1/e) = −1 (the fold), ÷ ln 10: 0, 0.246307, −0.434294 (twice).
- **T\*** = Ω·T_Planck = 8.0352·10³¹ K (paper quotes 8.0349·10³¹: a 4·10⁻⁵ slip); log₁₀ T\* = 31.90500 for both curves. The paper's 'match 3.05·10⁻¹⁶' is float64 rounding at an exact identity (here 0.0).

## What is definitional, what is not
- d* = Ω/ln 10 is a change of base: 'zero gap by construction' (wiki 17 says so). The crossing is exact given S_f = e^{−x}.
- **The step S_f = e^{−x} is asserted, not derived** (paper §III 'substituting the Bekenstein entropy'). Orientation check: the Boltzmann weight of a Planck-energy excitation at temperature T is e^{−T_P/T} = e^{−1/x}, and x = e^{−1/x} has NO real root (x ln x ≥ −1/e > −1). So the paper's substitution is not the Boltzmann weight; it is a different, assumed, law.
- The independent number is d*_spec = 0.24600 (measured, SPARC, Berry–Keating): gap 3.07·10⁻⁴ in d, 7.07·10⁻⁴ in Ω units. Not reproducible here (no SPARC data); still open.
- 'Hagedorn start 140 quadrillion °F' (`conjecture_proof.py` OMEGA_LIMIT = 140e15 °F ≈ 7.8·10¹⁶ K) shares its leading digits with T_Planck = 1.4168·10³² K but is 10¹⁵ lower: unit/scale slip or a different quantity; noted, not chased.

## Entropy and inertia as one pair of maps (the part that connects)
f(x) = e^{−x} (the entropy map, falling onto the ceiling) and its inverse f⁻¹(x) = −ln x (the inertia map, running away from it) share the fixed point Ω. Multipliers f′(Ω) = −Ω = −0.56714 and (f⁻¹)′(Ω) = −1/Ω = −1.76322: product 1 (measured 0.9999999979 from iterates), both alternate sides every step. Forward: contraction by Ω per step, the iterates 0.6, 0.549, 0.578, … alternate above and below. Backward: the same orbit read the other way, growing by 1/Ω per step. This is the saddle-pair (μ, 1/μ).
**The N-shape closes on this**: J_f = e^{−σE}, J_b = e^{−(1−σ)E} (wiki 73 §4) cross at σ = ½ at height e^{−E/2}; that height equals E/2 exactly when E = 2Ω (J = Ω, J_f·J_b = e^{−2Ω} = Ω²). Reading, not a theorem: the crossing of the two currents is the fixed point x = e^{−x} with x = E/2.
**Same statement as the Klein group** (docs/j2_definition): with E = ln n, J_f = n^{−σ} and J_b = n^{−(1−σ)} are the two terms n^{−s}, n^{−(1−s)} of the functional equation; their product is n^{−1} for every σ and they are equal only at σ = ½ (AM–GM). Exact, and the same identity as |x^{σ−½}|·|x^{½−σ}| = 1.

## The fractal that applies (UFformulary, 2026-11-11)
Of the 507 files in `PtolemyDesktop/Archimedes/Maths/Formula/UFformulary/` (Ultra Fractal) none computes the Riemann ζ or the critical-line zeros; the 'zeta' names are other things (Luschny's ZetaPolynomials in `jam3.ufm`; Toby Marshall's `mandalaZeta`; the 'Zeta-Like' type z = Σ aₖ^z in `mch.ufm`, a Dirichlet-polynomial iteration). Two are useful: (1) **Kerry Mitchell's `cardioid-julia` (`lkm.ufm`, comment at line 1259, formula at 1327)**: for z = c·e^z the main 'cardioid' boundary is c = k/e^k, k = e^{2πi m/n}. That is the multiplier circle, and the Lambert values are its real points: c = −1 (our entropy map, k = −Ω, inside), c = 1/e (k = +1, the cusp = the fold W(−1/e) = −1), c = −e (k = −1, period doubling: x → e^{−λx} flips at λ = e). `rhp.omega_dstar.exp_cardioid`, `test_exp_cardioid`. (2) The Riemann-sphere projection transforms (`pdf.uxf`, `reb.uxf`, `sam.uxf`, 'view from the South Pole, North Pole = ∞'), the same picture as our stereographic lift. `pdf.uxf`'s 'Lambert' is a map projection, unrelated to W.

