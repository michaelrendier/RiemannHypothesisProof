"""Nodal extraction with both sides mapped (rhp.two_sided).  R3: ψ and M node locking in 1-D to 10^7; R4: ξ and 1/ξ nodal lines in 2-D (reuses experiments/nodal_extract/xi_grid.npz).  Writes nodal_two_sided.json."""
import json, math, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
from rhp import two_sided as T, mertens_wave as W, nodal_extract as NX
rng = np.random.default_rng(0); R = {}
X = 10_000_000; mu = W.mobius_table(X); M = np.cumsum(mu.astype(np.int64)); psi, _ = NX.psi_table(X); n = np.arange(X + 1)
wp = (psi - n); wm = M
R["R3"] = {}
for lo, hi in ((10 ** 3, 10 ** 4), (10 ** 4, 10 ** 5), (10 ** 5, 10 ** 6), (10 ** 6, 10 ** 7)):
    sl = slice(lo, hi); idx = np.arange(lo, hi)
    nodes_p = T.sign_changes(wp[sl]) + lo; nodes_m = T.sign_changes(wm[sl]) + lo
    sp = np.sign(wp[sl]); sm = np.sign(wm[sl]); both = (sp != 0) & (sm != 0)
    m_in_p, _ = T.locking(nodes_p, nodes_m); p_in_m, _ = T.locking(nodes_m, nodes_p)
    nullA = T.null_locking(nodes_p, len(nodes_m), lo, hi, rng); nullB = T.null_locking(nodes_m, len(nodes_p), lo, hi, rng)
    R["R3"][f"{lo:.0e}-{hi:.0e}"] = dict(psi_nodes=int(len(nodes_p)), M_nodes=int(len(nodes_m)), nodes_per_M_node=float(len(nodes_p) / max(1, len(nodes_m))), M_in_psi_gaps_mean=m_in_p, M_in_psi_null=nullA,
                                         psi_in_M_gaps_mean=p_in_m, psi_in_M_null=nullB, opposite_sign_fraction=float(np.mean(sp[both] != sm[both])))
    print(f"{lo:.0e}-{hi:.0e}: ψ nodes {len(nodes_p)}, M nodes {len(nodes_m)}; M-node position in ψ gap {m_in_p:.3f} (null {nullA:.3f}); ψ-node in M gap {p_in_m:.3f} (null {nullB:.3f}); opposite sign {R['R3'][f'{lo:.0e}-{hi:.0e}']['opposite_sign_fraction']:.3f}", flush=True)
# R4: 2-D
d = np.load("experiments/nodal_extract/xi_grid.npz"); Zh, sig, ts = d["Z"], d["sig"], d["ts"]
sigfull = np.concatenate([1 - sig[:0:-1], sig]); Zfull = np.concatenate([np.conj(Zh[:, :0:-1]), Zh], axis=1)
inv = 1 / Zfull
agree_re = float(np.mean(np.sign(Zfull.real) == np.sign(inv.real))); agree_im = float(np.mean(np.sign(Zfull.imag) == -np.sign(inv.imag)))
Rl, Rc = NX.zero_set_2d(sigfull, ts, Zfull.real); Rl2, _ = NX.zero_set_2d(sigfull, ts, inv.real)
len1 = sum(NX.polyline_length(l) for l in Rl); len2 = sum(NX.polyline_length(l) for l in Rl2)
# angle of each Re = 0 rung where it crosses the equator σ = ½
eq = np.array([[0.5, ts[0]], [0.5, ts[-1]]]); angs = []; gam = []
for L in Rl:
    a0, a1 = L[:-1], L[1:]; s0 = a0[:, 0] - 0.5; s1 = a1[:, 0] - 0.5; k = np.nonzero(s0 * s1 <= 0)[0]
    for j in k:
        v = a1[j] - a0[j]
        if abs(v[0]) < 1e-14 and abs(v[1]) < 1e-14: continue
        t_cross = a0[j, 1] + (0.5 - a0[j, 0]) * v[1] / v[0] if v[0] != 0 else a0[j, 1]
        angs.append(math.degrees(math.atan2(v[1], v[0]))); gam.append(t_cross)
angs = np.array(angs); angs = np.where(angs > 90, angs - 180, np.where(angs < -90, angs + 180, angs))          # angle of the rung from the σ axis (horizontal), deg
R["R4"] = dict(grid=[len(ts), len(sigfull)], sign_agree_Re=agree_re, sign_agree_Im=agree_im, re_lines_xi=len(Rl), re_lines_inv=len(Rl2), re_length_xi=float(len1), re_length_inv=float(len2),
               rung_crossings=int(len(angs)), rung_angle_from_horizontal_deg=dict(median_abs=float(np.median(np.abs(angs))), max_abs=float(np.max(np.abs(angs))), mean=float(angs.mean())))
print("R4:", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in R["R4"].items() if k != "rung_angle_from_horizontal_deg"}, R["R4"]["rung_angle_from_horizontal_deg"], flush=True)
json.dump(R, open(pathlib.Path(__file__).with_suffix(".json"), "w"), indent=1)
