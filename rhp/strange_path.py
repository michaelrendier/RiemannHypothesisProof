# SPDX-License-Identifier: GPL-3.0-only
"""rhp.strange_path: a port of the 14-line `generateGeometry()` in the 'Holographic Animation Using JavaScript' short (Code with Bhurtel;
screenshots downloads/riemannsphere/javascript-path-1.jpeg and -2.jpeg), and the same geometry for the zeta function.

THE JAVASCRIPT, line by line (the screenshot shows only generateGeometry; rk4Step, CONFIG and the ODE are NOT shown):
  state = {x: 0.1, y: 0, z: 0}              one point in a 3-D phase space, started near the origin
  for i < CONFIG.burnIn: rk4Step(state)     run the ODE with 4th-order Runge-Kutta and throw the output away, so the point reaches the attractor
  for i < CONFIG.points: rk4Step(state)     then record the orbit
      positions[3i]   = state.x * scale     x stays x
      positions[3i+1] = state.z * scale     the SECOND slot gets z and the THIRD gets y (y and z swapped, so z is drawn upward)
      positions[3i+2] = state.y * scale
      lineDistances[i] = i / points         0 -> 1 along the path: the renderer uses it to colour and to dash the line (the rainbow fade in the images)
So the picture is ONE trajectory of a 3-D ODE, drawn as a line.  It is a strange attractor.  The screenshots look like the Aizawa attractor (a sphere with a tube through
its axis); the code shown does not name the system, so `aizawa` below is an IDENTIFICATION BY SHAPE, checked in experiments/strange_path.py against a render.
IT IS NOT THE ZETA FUNCTION.  Nothing in the shown code evaluates ζ; an RK4 orbit of Aizawa's equations has no zeros, no critical line, no primes.

THE ZETA VERSION (built here, same layout): the critical-line path t ↦ ζ(½ + it), put on the Riemann sphere by inverse stereographic projection
      (X, Y, Z) = (2 Re w, 2 Im w, |w|² − 1) / (|w|² + 1),   w = ζ(½ + it)
Every zero is w = 0, which lands on the south pole (0, 0, −1): the path visits the pole once per zero.  w = ∞ (a pole of ζ) is the north pole; on the critical line ζ has none.
The positions array uses the same slot order and the same lineDistances as the JavaScript, so one renderer can draw either.
Order of work for readers: `rk4_step`, `generate_geometry` (the port), `aizawa` (one system), `zeta_sphere_path` (the zeta geometry).
"""
from __future__ import annotations
import numpy as np
from mpmath import mp, mpf, mpc, zeta

# ── the port ──
def rk4_step(f, state, h):
    """One classical Runge-Kutta step of x' = f(x) (4th order: local error O(h^5))."""
    k1 = f(state); k2 = f(state + 0.5 * h * k1); k3 = f(state + 0.5 * h * k2); k4 = f(state + h * k3)
    return state + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

def aizawa(p=(0.95, 0.7, 0.6, 3.5, 0.25, 0.1)):
    """Aizawa attractor, standard parameters a b c d e f."""
    a, b, c, d, e, f = p
    def rhs(s):
        x, y, z = s
        return np.array([(z - b) * x - d * y, d * x + (z - b) * y, c + a * z - z ** 3 / 3 - (x * x + y * y) * (1 + e * z) + f * z * x ** 3])
    return rhs

def generate_geometry(rhs, points=20000, burn_in=2000, scale=1.0, h=0.01, start=(0.1, 0.0, 0.0)):
    """The port of generateGeometry(): returns (positions float32 [points*3] in the JS slot order x, z, y; line_distances float32 [points])."""
    state = np.array(start, float)
    for _ in range(burn_in): state = rk4_step(rhs, state, h)
    positions = np.zeros(points * 3, np.float32); line = np.zeros(points, np.float32)
    for i in range(points):
        state = rk4_step(rhs, state, h)
        positions[3 * i] = state[0] * scale; positions[3 * i + 1] = state[2] * scale; positions[3 * i + 2] = state[1] * scale; line[i] = i / points
    return positions, line

# ── the zeta version ──
def zeta_critical(t):
    return complex(zeta(mpc(0.5, t)))

def sphere_point(w):
    """Inverse stereographic projection of w onto the unit sphere; w = 0 -> south pole (0,0,-1), w = inf -> north pole."""
    m = abs(w) ** 2; return np.array([2 * w.real, 2 * w.imag, m - 1.0]) / (m + 1.0)

def zeta_sphere_path(points=6000, t0=0.0, t1=60.0, scale=1.0):
    """Positions in the JS slot order (X, Z, Y) and lineDistances, for the critical-line image of zeta on the Riemann sphere."""
    mp.dps = 15; ts = np.linspace(t0, t1, points)
    positions = np.zeros(points * 3, np.float32); line = np.zeros(points, np.float32); w = np.zeros(points, complex)
    for i, t in enumerate(ts):
        w[i] = zeta_critical(float(t)); X, Y, Z = sphere_point(w[i]) * scale
        positions[3 * i], positions[3 * i + 1], positions[3 * i + 2] = X, Z, Y; line[i] = i / points
    return positions, line, ts, w
