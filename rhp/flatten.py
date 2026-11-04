# SPDX-License-Identifier: GPL-3.0-only
"""Artefact classifier: what the flat plot of ζ(½+it) does to the pathway.

Each function takes the dict from rhp.uns.pathway() and returns plain rows.
Artefact types (see notebooks/01): hub, false crossing, petal sign, petal size,
lost clock, missed petal (Gram's law), mirror half.
"""
from __future__ import annotations

import numpy as np
from mpmath import mp, mpf, mpc, zeta

from . import uns
from .viewport import Viewport, assert_viewport_sound


def hub_passes(p, zeros):
    """One row per zero: the direction the flat curve crosses the origin vs −θ(γ) mod π.

    direction = arg(i·ζ′(½+iγ)) mod π  (exact derivative from mpmath).  CALCULATED: equal to −θ(γ) mod π.
    """
    rows = []
    for n, g in enumerate(zeros, 1):
        s = mpc(0.5, g)
        d = 1j * complex(zeta(s, derivative=1))
        direction = float(np.angle(d) % np.pi)
        pred = float((-uns.theta(g)) % np.pi)
        # sign of Z just before/after the zero (petal sign flip)
        zb = float(np.real(np.exp(1j * uns.theta(g - 1e-3)) * uns.zeta_line(g - 1e-3)))
        za = float(np.real(np.exp(1j * uns.theta(g + 1e-3)) * uns.zeta_line(g + 1e-3)))
        rows.append(dict(n=n, gamma=g, direction=direction, minus_theta_mod_pi=pred,
                         agree=abs(direction - pred) < 1e-6, speed=abs(d),
                         sign_before=int(np.sign(zb)), sign_after=int(np.sign(za))))
    return rows


def _segments_intersect(P, Q, i0, i1):
    """Vectorised proper-intersection test of segment rows i0:i1 against all later segments."""
    A, B = P[i0:i1], P[i0 + 1:i1 + 1]
    C, D = P[:-1], P[1:]
    out = []
    r = B - A
    s = D - C
    # solve A + a r = C + b s for each (i, j)
    rxs = r[:, None, 0] * s[None, :, 1] - r[:, None, 1] * s[None, :, 0]
    q = C[None, :, :] - A[:, None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        a = (q[..., 0] * s[None, :, 1] - q[..., 1] * s[None, :, 0]) / rxs
        b = (q[..., 0] * r[:, None, 1] - q[..., 1] * r[:, None, 0]) / rxs
    hit = (rxs != 0) & (a >= 0) & (a < 1) & (b >= 0) & (b < 1)
    ii, jj = np.nonzero(hit)
    for i, j in zip(ii, jj):
        gi = i0 + i
        if j > gi + 1:                              # not adjacent, count each pair once
            out.append((gi, j, a[i, j]))
    return out


def false_crossings(p, hub_eps=0.05, chunk=250):
    """Self-intersections of the flat polyline. Rows: t1, t2, point, hub (|ζ| < hub_eps).

    A crossing is a 'false crossing' (artefact) when t1 ≠ t2 — in the 3D spring the two parameters
    are always a distance |t2 − t1| apart; the flat plot erases it.
    """
    P = np.column_stack([p["zeta"].real, p["zeta"].imag])
    t = p["t"]
    rows = []
    for i0 in range(0, len(P) - 1, chunk):
        for gi, gj, a in _segments_intersect(P, None, i0, min(i0 + chunk, len(P) - 1)):
            pt = P[gi] + a * (P[gi + 1] - P[gi])
            rows.append(dict(t1=float(t[gi]), t2=float(t[gj]), x=float(pt[0]), y=float(pt[1]),
                             hub=bool(np.hypot(*pt) < hub_eps), dt=float(t[gj] - t[gi])))
    return rows


def petals(p, zeros):
    """One row per petal (between consecutive zeros, plus the lead-in): linear vs log size."""
    edges = np.concatenate([[p["t"][0]], zeros, [p["t"][-1]]])
    rows = []
    for k in range(len(edges) - 1):
        m = (p["t"] >= edges[k]) & (p["t"] <= edges[k + 1])
        if not m.any():
            continue
        a = np.abs(p["zeta"][m])
        rows.append(dict(petal=k, t_from=float(edges[k]), t_to=float(edges[k + 1]),
                         max_abs=float(a.max()), max_log10=float(np.log10(a.max()))))
    return rows


def lost_clock(p, zeros, guard=0.2):
    """Winding rate: numeric dφ/dt of the continuous angle vs exact −θ′(t), away from zeros.

    φ is the arg ζ unwrapped with the π jump at each zero removed (signed-radius form).
    """
    t = p["t"]
    ang = np.unwrap(np.angle(p["zeta"] * np.where(p["Z"] < 0, -1, 1)))   # remove sign-of-Z jump
    rate = np.gradient(ang, t)
    far = np.all(np.abs(t[:, None] - zeros[None, :]) > guard, axis=1) if len(zeros) else np.ones_like(t, bool)
    far[:2] = far[-2:] = False
    return dict(max_abs_dev=float(np.max(np.abs(rate[far] + p["thetap"][far]))),
                min_rate=float(np.min(-p["thetap"])), sign_change_t=float(t[np.argmax(p["thetap"] > 0)]))


def gram_law(T=300.0):
    """Zero counts in each Gram interval [gₙ, gₙ₊₁), n = 0..; Gram's law fails where the count ≠ 1."""
    zs = uns.zeros_below(T)
    gp = uns.gram_points(T)
    counts = [int(np.sum((zs >= gp[i]) & (zs < gp[i + 1]))) for i in range(len(gp) - 1)]
    bad = [dict(n=i, g_lo=float(gp[i]), g_hi=float(gp[i + 1]), zeros=c)
           for i, c in enumerate(counts) if c != 1]
    return dict(T=T, intervals=len(counts), failures=bad)


def mirror_check(ts=(3.0, 14.1, 27.7, 41.3)):
    """ζ(½−it) = conj ζ(½+it): the t<0 half is the conjugate image."""
    return max(abs(complex(zeta(mpc(0.5, -mpf(x)))) - uns.zeta_line(x).conjugate()) for x in ts)


def window_map(p, zeros, crossings, vh=500):
    """Scroll a t-window over the samples with the nes-viewport model; one summary row per window.

    Sweep soundness (clamp/flush/thumb/sweep) is asserted; the union of windows must be every sample.
    """
    n = len(p["t"])
    vp = Viewport(cw=1, ch=n, vw=1, vh=vh)
    assert_viewport_sound(vp, track_len=20)
    rows, seen = [], set()
    vp.move_to_top()
    while True:
        r = vp.visible_rows()
        seen |= set(r)
        ta, tb = p["t"][r.start], p["t"][r.stop - 1]
        rows.append(dict(oy=vp.oy, t_from=float(ta), t_to=float(tb),
                         zeros=int(np.sum((zeros >= ta) & (zeros <= tb))),
                         false_crossings=int(sum(1 for c in crossings if not c["hub"] and ta <= c["t1"] <= tb)),
                         hub_crossings=int(sum(1 for c in crossings if c["hub"] and ta <= c["t1"] <= tb))))
        if vp.oy >= vp.max_oy:
            break
        vp.page_down()
    assert seen == set(range(n))
    return rows
