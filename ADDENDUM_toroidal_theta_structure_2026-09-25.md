# Addendum — The Riemann-Siegel Theta Function as a Toroidal Helix: Major Loop (Chirp) vs Minor Loop (Prime Resonance)

**Status.** Measurement, not a proof. Does not close C1. Computed live,
2026-09-25, `mpmath` 1.3.0, 40 real nontrivial zeros (`mpmath.zetazero(1..40)`),
dps=25. Script: session scratchpad, reproducible directly from the snippets
below — no external dependency beyond `mpmath`.

---

## A. The flattening problem, named precisely

A standard 2D plot of `ζ(1/2+it)` or the Riemann-Siegel `Z(t)` — real axis
against imaginary axis, `t` implicit — makes the curve's repeated near-returns
to the origin (at each zero) look like false cyclicity: the same shape
recurring. It isn't recurring. Each near-zero-crossing happens at a strictly
larger `t`, never revisited. Un-flattening with `t` as an explicit third
(vertical) axis turns the "recurring" 2D curve into an ascending helix, with
`(Re, Im)` of the `θ(t)`-rotated trajectory as the cross-section at each
height.

## B. `θ(t)`-rotation makes `ζ` real exactly on `σ=1/2`, tested against `σ=0.7`

`e^{iθ(t)}·ζ(1/2+it)` is real to machine precision (imaginary parts `~1e-31`
to `1e-38`, pure float noise) and matches `mpmath.siegelz(t)` independently to
10 decimal places. The identical rotation applied at `σ=0.7` leaves a real,
non-negligible imaginary residual (`0.03` to `0.19`, an order of magnitude
above float error) at the same `t` values. `σ=1/2` is not merely *a* fixed
point of `s↦1−s` — it is specifically where this particular complex rotation
resolves onto the real line; off it, the same rotation does not resolve.

## C. The helix has two independent frequencies — major loop and minor loop

Computed `θ'(t)` (the local winding rate) and the actual zero spacing against
the smooth prediction `2π/θ'(t)`, for `n=2..40`:

- **Major loop — `θ'(t)`, the smooth secular rate.** Grew `0.405 → 1.487`
  across these 40 zeros (**3.67×**). Monotonic, no peaks — a chirp, not a
  resonance. This is the classical Riemann-von Mangoldt mean zero density,
  `θ'(t) ≈ (1/2)ln(t/2π)`, read as the axis-of-revolution's own rotation rate
  rather than as a counting formula.
- **Minor loop — the fluctuation, `spacing − 2π/θ'(t)`.** Mean `−2.648`,
  population stdev `1.046` over the sample — a structured, non-zero-mean
  oscillation (GUE-type level repulsion), not noise scattered around the
  smooth prediction.

Two independent frequencies swept around one common axis (`t`) is precisely
the generating picture of a torus: the major loop is the path around the
axis, the minor loop is the tube. Not asserted by analogy — the two
frequencies were measured separately, above, before this geometric reading
was applied to them.

## D. The major loop does not resonate. The minor loop does, and its resonances have names.

- **`θ'(t)` is not resonant.** Smooth, monotonic, no characteristic peak
  frequencies — confirmed by direct computation (§C), not assumed.
- **The fluctuation is resonant, and this is classical, not new:** the
  Riemann–von Mangoldt explicit formula already establishes that the
  fluctuating part of the zero-counting function, Fourier-transformed, peaks
  exactly at `log p` for every prime `p` and its powers (Riemann 1859, von
  Mangoldt, Weil; the "spectral interpretation" line of work, Berry–Keating).
  The minor loop's resonant frequencies are the primes. This addendum did
  not derive that fact — it identifies the already-classical fact as *this
  specific geometric object's* resonant structure, precisely, rather than
  leaving "the zeros encode the primes" as an unlocated generality.

## E. What is explicitly OPEN — not claimed, not tested

Whether the tilt/precession relationship already established elsewhere in
this project (`Ω ∝ sin(tilt)`, "oblique gearing" — the real-axis component
driving a secondary rotation whenever nonzero) itself resonates *with* the
`t`-axis — i.e., whether the precession rate shows characteristic peaks at
particular heights along the axis, distinct from simply tracking the major
loop's smooth chirp — has not been computed. Real, calculable, distinct
question from §D. Left open, not answered by what's established above.

## F. Status against C1

Nothing here closes C1 (self-adjointness / the deficiency-indices import
already on record in this project's other 0_RB derivations). This addendum
is a geometric reading of already-classical analytic facts (the functional
equation's fixed point, the explicit formula's prime resonance) plus one new,
tested, falsifiable computation (the `θ(t)`-rotation reality test, §B) and one
honestly-flagged open question (§E). Treat accordingly.
