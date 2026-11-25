"""Pole-view density of the zeta path, the pseudo-cardioid decomposition, and the Ulam-spiral diagonals.  Writes pole_density_and_ulam.json.
Expectations written before the run:
  P1  seen from the pole the path is NOT uniform in azimuth: the mean value of zeta(1/2+it) over t is 1 (the n = 1 term), so the path leans to Re > 0; the less dense region is
      the Re < 0 side; it does not shrink as more zeros are added (so it is not a truncation artefact of uncalculated primes).
  P2  zeta = S_N + exp(-2 i theta) conj(S_N) + R with S_N = sum_{n<=N} n^(-1/2-it), N = floor(sqrt(t/2pi)) and |R| <= (t/2pi)^(-1/4) |C0|: a circle through 0 for N = 1, then a pseudo-cardioid without a cusp
      (at a simple zero the speed |Z'| is not 0).
  U1  every straight diagonal piece of the square spiral has second difference 8 (f(k) = 4k^2 + bk + c), so every diagonal is 4(k + b/8)^2 - D/16, D = b^2 - 16c: the perfect squares shifted by a constant.
  U2  in the Sacks picture (r = sqrt(m), angle = 2 pi sqrt(m)) the diagonal f(k) = 4k^2 + bk + c sits at the fixed angle pi b / 2 (mod 2 pi), bending by D / (32 k)."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.zeta_uns_path import maths as M
R = {}
g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()])
ts = M.sample_by_phase(10.0, 25000, 16, jitter=True); z = M.zeta_line(ts); th = M.theta(ts)
dt = np.gradient(ts); P1 = {}
for nz in (100, 1000, 5000, 25000):
    k = nz * 16 + 1; zz = z[:k]; w_t = dt[:k] / dt[:k].sum()
    mean_t = complex(np.sum(zz * w_t)); ang = np.angle(zz)
    H, _ = np.histogram(ang, bins=72, range=(-math.pi, math.pi)); dens = H / H.mean(); less = dens < 0.5
    Ht, _ = np.histogram(ang, bins=72, range=(-math.pi, math.pi), weights=w_t); dt_ = Ht / Ht.mean()
    P1[nz] = dict(mean_value_t_weighted=[mean_t.real, mean_t.imag], frac_Re_neg_phase_uniform=float(np.mean(zz.real < 0)), frac_Re_neg_t_weighted=float(np.sum(w_t[zz.real < 0])),
                  azimuth_frac_below_half_density=float(less.mean()), azimuth_frac_below_half_density_t_weighted=float((dt_ < 0.5).mean()),
                  min_density_azimuth_deg=float(np.degrees(-math.pi + (np.argmin(dens) + 0.5) * 2 * math.pi / 72)))
R['P1'] = P1
# the n = 1 term alone: where would the path be?
c1 = 1 + np.exp(-2j * th); R['P1_n1_circle'] = dict(frac_Re_neg=float(np.mean(c1.real < -1e-9)), min_Re=float(c1.real.min()))
# P2
tt = np.linspace(200, 20000, 600); zt = M.zeta_line(tt); th2 = M.theta(tt); Nn = np.floor(np.sqrt(tt / (2 * math.pi))).astype(int); resid = []; bound = []
for i, t in enumerate(tt):
    n = np.arange(1, Nn[i] + 1); S = np.sum(n ** (-0.5 - 1j * t)); rec = S + np.exp(-2j * th2[i]) * np.conj(S); resid.append(abs(zt[i] - rec)); bound.append((t / (2 * math.pi)) ** -0.25 * 1.0)
resid = np.array(resid); bound = np.array(bound); R['P2'] = dict(max_residual=float(resid.max()), max_ratio_to_bound=float((resid / bound).max()), n=len(tt))
zer = g[:1500]; h = 1e-4; zp = np.array([(float(M.Z_rs(x + h)[0]) - float(M.Z_rs(x - h)[0])) / (2 * h) for x in zer]); R['P2']['min_abs_Zprime_at_first_1500_zeros'] = float(np.abs(zp).min()); R['P2']['median_abs_Zprime'] = float(np.median(np.abs(zp)))
# Ulam square spiral
SZ = 301; X = np.zeros(SZ * SZ + 1, int); Y = np.zeros(SZ * SZ + 1, int); x = y = 0; dx, dy = 1, 0; step = 1; m = 1; X[1] = Y[1] = 0
mm = 1
while mm < SZ * SZ:
    for _ in range(2):
        for _ in range(step):
            if mm >= SZ * SZ: break
            x += dx; y += dy; mm += 1; X[mm] = x; Y[mm] = y
        dx, dy = -dy, dx
    step += 1
half = SZ // 2; G = np.zeros((SZ, SZ), np.int64)
for m_ in range(1, SZ * SZ + 1):
    if abs(X[m_]) <= half and abs(Y[m_]) <= half: G[X[m_] + half, Y[m_] + half] = m_
S_ = np.ones(SZ * SZ + 1, bool); S_[:2] = False
for i in range(2, int((SZ * SZ) ** 0.5) + 1):
    if S_[i]: S_[i * i::i] = False
pieces = []
for sgn in (1, -1):
    for c in range(-2 * half + 20, 2 * half - 20):
        seq = []
        for xx in range(-half, half + 1):
            yy = sgn * xx + c if sgn == 1 else -xx + c
            yy = xx + c if sgn == 1 else -xx + c
            if -half <= yy <= half: seq.append(G[xx + half, yy + half])
        seq = np.array([v for v in seq if v > 0]); 
        if len(seq) < 8: continue
        d2 = np.diff(seq, 2); i0 = 0
        for i in range(1, len(d2) + 1):
            if i == len(d2) or d2[i] != d2[i0]:
                if i - i0 >= 6: pieces.append((sgn, c, int(d2[i0]), seq[i0:i + 2]))
                i0 = i
d2s = [p[2] for p in pieces]; vals, cnts = np.unique(d2s, return_counts=True); R['U1'] = dict(pieces=len(pieces), second_difference_counts={int(v): int(n) for v, n in zip(vals, cnts)})
# fit pieces with d2 = +-8 : f(k) = (d2/2) k^2 + b k + c ; k = 0.. ; discriminant and prime fraction
rows = []
for sgn, c, d2, seq in pieces:
    if abs(d2) != 8: continue
    a = d2 // 2; k = np.arange(len(seq)); A = np.vstack([k ** 2, k, np.ones_like(k)]).T; coef = np.linalg.lstsq(A, seq.astype(float), rcond=None)[0]; b, c0 = coef[1], coef[2]
    rows.append(dict(a=a, b=round(b), c=round(c0), D=round(b * b - 4 * a * c0), n=len(seq), primes=int(S_[seq].sum()), x_identity_err=float(np.max(np.abs(seq - (a * (k + b / (2 * a)) ** 2 - (b * b - 4 * a * c0) / (4 * a)))))))
R['U1']['fitted_pieces'] = len(rows); R['U1']['identity_max_err'] = max(r['x_identity_err'] for r in rows) if rows else None
big = sorted(rows, key=lambda r: -r['primes'])[:6]; R['U1']['most_prime_rich'] = big
frac = np.array([r['primes'] / r['n'] for r in rows]); Dv = np.array([r['D'] for r in rows]); R['U1']['prime_fraction_mean'] = float(frac.mean())
# U2 Sacks angle of central diagonals
U2 = []
for b, c in ((2, 1), (4, 1), (6, 3), (0, 1), (8, 5), (-2, 1)):
    k = np.arange(20, 2000); f = 4.0 * k * k + b * k + c; ang = np.mod(2 * math.pi * np.sqrt(f), 2 * math.pi); D = b * b - 16 * c
    pred = (math.pi * b / 2) % (2 * math.pi); err = np.angle(np.exp(1j * (ang - pred)))
    U2.append(dict(b=b, c=c, D=D, predicted_angle=pred, mean_abs_err_k100_2000=float(np.mean(np.abs(err[80:]))), err_times_k_median=float(np.median(err * k)), predicted_err_times_k=float(-2 * math.pi * D / 32)))
R['U2'] = U2
json.dump(R, open(root / 'experiments/pole_density_and_ulam.json', 'w'), indent=1, default=float)
print(json.dumps({k: R[k] for k in ('P1', 'P1_n1_circle', 'P2')}, indent=1)); print(json.dumps({'U1': {k: v for k, v in R['U1'].items() if k != 'most_prime_rich'}, 'U2': R['U2']}, indent=1)); print(R['U1']['most_prime_rich'][:4])
