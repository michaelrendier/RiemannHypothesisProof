# SPDX-License-Identifier: GPL-3.0-only
"""The trig hexagon (rhp.trig_hexagon): six functions on a ring, 1 at the centre; what the rotations of that ring preserve.

The picture (MATH Talk short, saved as downloads/trig-identity-1.jpeg and -2.jpeg): vertices around the hexagon in the order
    sin, cos, cot, csc, sec, tan        and 1 in the centre.
ESTABLISHED (textbook identities; checked below by residuals at random angles):
  opposite vertices are reciprocals          sin·csc = cos·sec = tan·cot = 1
  every vertex is the product of its two ring neighbours      tan = sin·sec,  sin = tan·cos,  cos = sin·cot, ...
  the three shaded triangles (centre + two adjacent vertices) are the Pythagorean identities
      sin² + cos² = 1,    1 + cot² = csc²,    1 + tan² = sec²
  the quotient identities tan = sin/cos, cot = cos/sin follow from the product rule.
CALCULATED here (the question 'can all of trig be represented with some rotations?'):
  1. In log space the six vertices are ±u, ±v, ±(u−v) with u = ln sin, v = ln cos: the six roots of A2, and the product rule is vector addition
     (two neighbours 60° apart add to the vertex between them).  The ring's rotation by 60° is the Weyl-group-with-sign action on A2.
  2. Which rotations/reflections of the ring leave each family invariant (`invariance()`, MEASURED): reciprocal and product rules: all 12 symmetries of the hexagon;
     Pythagorean forms: only the identity and the co-function reflection θ ↦ π/2 − θ (the hypotenuse role differs per triangle, so no rotation preserves them).
  3. The sum, difference, double-angle and De Moivre formulas are NOT on the hexagon; they are the group law of the rotation R(θ) of the plane
     (R(a)R(b) = R(a+b)); sin and cos are its entries and the other four functions are the hexagon's products and quotients of those two.
     So 'all of trig' = the rotation group (two entries) + the hexagon (the other four).  Half-angle needs a square root and is outside both.
"""
from __future__ import annotations
import itertools
import numpy as np

RING = ('sin', 'cos', 'cot', 'csc', 'sec', 'tan')
FN = {'sin': np.sin, 'cos': np.cos, 'tan': np.tan, 'cot': lambda x: 1 / np.tan(x), 'sec': lambda x: 1 / np.cos(x), 'csc': lambda x: 1 / np.sin(x)}
OPPOSITE = [(RING[i], RING[i + 3]) for i in range(3)]
PYTHAG = [('sin', 'cos', None), ('tan', None, 'sec'), ('cot', None, 'csc')]     # (a, b, c): a² + b² = c², None stands for the centre value 1; sin²+cos²=1: c is the centre

def values(theta):
    return {k: FN[k](theta) for k in RING}

def neighbours(k):
    i = RING.index(k); return RING[(i - 1) % 6], RING[(i + 1) % 6]

# identities as predicates over a labelling v: label -> value (a labelling = a permutation of the ring names)
def reciprocal_ok(v, pairs=OPPOSITE, tol=1e-9):
    return all(abs(v[a] * v[b] - 1) < tol for a, b in pairs)
def product_ok(v, tol=1e-9):
    return all(abs(v[k] - np.prod([v[n] for n in neighbours(k)])) < tol for k in RING)
def pythag_ok(v, tol=1e-9):
    return all(abs(v['sin'] ** 2 + v['cos'] ** 2 - 1) < tol for _ in [0]) and abs(1 + v['tan'] ** 2 - v['sec'] ** 2) < tol and abs(1 + v['cot'] ** 2 - v['csc'] ** 2) < tol

def relabel(v, k, reflect=False):
    """Rotate the ring labels by k steps of 60° (optionally reflect first): the value that sat at position i is now read as the label at position i+k."""
    ring = RING[::-1] if reflect else RING
    out = {}
    for i, name in enumerate(RING):
        out[ring[(i + k) % 6] if reflect else RING[(i + k) % 6]] = v[name]
    return out

def invariance(n_angles=50, seed=3):
    """For each of the 12 symmetries of the hexagon (6 rotations, 6 reflections): do the three families of identities still hold when the values are moved by it?
       Returns {family: [(symmetry, holds_for_all_angles), ...]}."""
    rng = np.random.default_rng(seed); angles = rng.uniform(0.2, 1.3, n_angles)
    fams = {'reciprocal': reciprocal_ok, 'product': product_ok, 'pythagorean': pythag_ok}
    out = {f: [] for f in fams}
    for refl, k in itertools.product((False, True), range(6)):
        for f, test in fams.items():
            ok = all(test(relabel(values(t), k, refl)) for t in angles)
            out[f].append((('reflect+' if refl else 'rotate ') + f'{60 * k}°', ok))
    return out

def log_roots(theta):
    """The six vertices in log space as vectors in the (u, v) plane, u = ln sin, v = ln cos (angle in (0, π/2)); returns {name: (x, y)} in the A2 picture."""
    u, v = np.log(np.sin(theta)), np.log(np.cos(theta))
    return {'sin': (u, 0.0), 'cos': (v, 0.0), 'tan': (u - v, 0.0), 'cot': (v - u, 0.0), 'sec': (-v, 0.0), 'csc': (-u, 0.0)}

def a2_coordinates(theta):
    """ln of each vertex = a·(ln sin) + b·(ln cos) with integer (a, b): the six vectors ±(1,0), ±(0,1), ±(1,−1).  Rotation by 60° acts on these as a linear map of order 6."""
    return {'sin': (1, 0), 'cos': (0, 1), 'cot': (-1, 1), 'csc': (-1, 0), 'sec': (0, -1), 'tan': (1, -1)}

def rotation_matrix_action():
    """The 6×6 permutation of the ring by 60°, and the 2×2 integer matrix R on (a, b) that realises it on the A2 coordinates; checks R^6 = I and the neighbour sum rule."""
    A = a2_coordinates(0.5); M = np.array([[0, -1], [1, 1]])          # columns = images of (1,0) and (0,1): sin->cos, cos->cot
    out = {}
    for i, k in enumerate(RING):
        nxt = RING[(i + 1) % 6]; out[k] = tuple(M @ np.array(A[k]))
    ok_map = all(out[k] == A[RING[(i + 1) % 6]] for i, k in enumerate(RING))
    ok_order = np.array_equal(np.linalg.matrix_power(M, 6), np.eye(2, dtype=int))
    ok_sum = all(tuple(np.array(A[neighbours(k)[0]]) + np.array(A[neighbours(k)[1]])) == A[k] for k in RING)
    return M, ok_map, ok_order, ok_sum

# ── the rotation group is the rest of trig ──
def R(theta):
    c, s = np.cos(theta), np.sin(theta); return np.array([[c, -s], [s, c]])
def six_from_rotation(theta):
    """All six functions from the matrix R(θ): sin and cos are its entries, the other four are the hexagon's products/quotients."""
    r = R(theta); s, c = r[1, 0], r[0, 0]
    return {'sin': s, 'cos': c, 'tan': s / c, 'cot': c / s, 'sec': 1 / c, 'csc': 1 / s}
def sum_formulas_from_rotation(a, b):
    """R(a)R(b) = R(a+b): read sin(a+b), cos(a+b) off the matrix product, then tan(a+b) from the hexagon quotient."""
    r = R(a) @ R(b); s, c = r[1, 0], r[0, 0]; return {'sin': s, 'cos': c, 'tan': s / c}
def double_angle_from_rotation(a):
    return sum_formulas_from_rotation(a, a)
def de_moivre(theta, n):
    return np.linalg.matrix_power(R(theta), n)
