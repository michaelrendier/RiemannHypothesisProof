# The pathway of the first K primes (2026-11-08)

Code `rhp/primepath.py`; figures `rhp/atlas/primepath_flat.svg`, `primepath_3axis.svg`, `primepath_uns.svg` (`python -m rhp.atlas.primepath`); tests `test_primepath_*`.
Path: log P_K(s) = Σ_{k≤K} −ln(1 − p_k^{−s}), s = ½ + it, K up to 10⁶ primes (p_K = 15,485,863). Labels: ESTABLISHED / OURS-CALCULATED / READING.

## Corrections to the premises (measured)
- **After the first 65 primes the survivors are not only primes.** Density of integers with no prime factor ≤ 313: Π(1−1/p) = 0.0965 (Mertens e^{−γ}/ln 313 = 0.0977). In [10⁶, 2·10⁶) 9.62 % survive,
  of which 73.2 % are prime and 26.8 % are composites with every factor > 313 (e.g. 317²). The 65 primes bring down every composite only below 317² = 100,489 (exhaustively verified, log item 15).
  A number N needs the primes up to √N, so the set of primes used as factors is infinite.
- **Both sets are countable (ℵ₀).** What differs is natural density: primes 0, composites 1. The 'larger infinity' is a density statement.

## The path as K grows (does it keep growing?)
- Flat: at t = γ₁ the first 65 primes give a jagged tangle (|log P| = 2.6); afterwards the path is a smooth OUTWARD spiral that never closes: |log P| = 3.8 / 8.9 / 13.2 at K = 10⁴ / 10⁵ / 10⁶.
  At t = 3 it is 9.5 / 27.3 / 78.5. The prime sum tracks the prime-number-theorem part Ei((1−s) ln x): ratio 1.001 at t = 3, K = 10⁵. Growth ~ x^{½}/(|1−s| ln x), turning at rate −t ln x: a logarithmic spiral.
- 3 axes (Re, Im, ln p): a cone/funnel opening with ln p (the sieve clock) — a conical helix; the first 65 primes are the jagged tip.
- Why 65 matters for the SHAPE: below x ≈ (t·ln x)² the smooth part is smaller than the random-phase scatter, so the tangle dominates; the outward spiral takes over later (~10⁴ primes at t = 14).
- The torus: K wheels are circles in K orthogonal complex planes, a point of a flat torus T^K ⊂ C^K; the path is the sum of coordinates (a projection). The sphere S^{2K−1} holding it (rms sense) has
  radius R_K = √Σ1/p = 1.422 / 1.646 / 1.705 / 1.752 at K = 65 / 10⁴ / 10⁵ / 10⁶: grows like √ln ln p, without bound but absurdly slowly.
- On the Riemann sphere (|P| = 1 the equator): at t = γ₁, ln|P_K| stays near −2.5 (north, product near 0) for ~10⁴ primes, then oscillates with growing amplitude (±19 by 10⁶), crossing |P| = 1 17 times.
  The product does NOT converge to ζ(½+iγ₁) = 0 within 10⁶ primes (|P_K| ranges from 6·10⁻⁹ to 5·10⁵): the Euler product on the critical line is not convergent in the usual sense (ESTABLISHED).
- The s-point itself never leaves the equator: s = ½ + it has |w| = 1 for every K. The primes act through the VALUE P_K, not through s.

## Precision (what resolution matters)
Reconstruct W_T (u ∈ [6, 14.5], x ≤ 2·10⁶) from the first 1000 zero heights rounded to d decimals — correlation with the sieve wave: d=0 0.234, 1 0.941, 2 0.9834, 3 0.9841, 8 0.9841.
Rule: a height error δγ shifts a tone's phase by δγ·u; coherence needs δγ·u ≲ 0.1, so ~2 decimals up to u ≈ 15 and one more digit per factor of 10 in x… (δγ ≲ 0.1/ln x). Finite precision, finite tones: no uncountable wall appears.
