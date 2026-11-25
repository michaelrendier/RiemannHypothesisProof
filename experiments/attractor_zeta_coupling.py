"""The zeta function and the Aizawa attractor run through each other; the two precessions and the twin primes; shadows of lines.  Writes attractor_zeta_coupling.json and
visualizations/javascript/attractor_zeta_data.js (small data for the page).  Which part to run: python attractor_zeta_coupling.py A B C D
Expectations written before the run:
  A  The attractor has two frequencies, the rotation d about z and the (r, z) oscillation w_rz.  At f = 0 the (r, z) system is a planar flow, so it settles to a limit cycle (Poincare-Bendixson) and the 3-D orbit is a torus
     with rotation number rho = d / w_rz.  With f = 0.1 the torus is broken into the chaotic attractor.  Sampling the orbit at the twin primes (n and n+2 prime) is expected to show NO phase alignment (Rayleigh p large)
     and no difference of z from matched primes; the ratios of the two frequencies are NOT expected to equal the twin-prime constants to 1% (reported with how many constants were tried).
  B  Reading the orbit as s = x + i K y (K = 19) puts the x = 1/2 plane at the critical line; the orbit crosses it many times; the heights K y at the crossings are smooth (not clustered on the zeros).
  C  The rotation rate of the zeta path, theta'(t), equals the Aizawa d = 3.5 at t = 2 pi e^7 = 6889.5, about zero number 6,5xx, step 66 of the animation.
  D  Of the vertical lines sigma = const, only sigma = 1/2 has an image that reaches 0 (the others stay away from 0); horizontal lines Im s = c (c = 1/2, 1, pi/2, pi, 2 pi) have images that miss 0."""
import base64, json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.aizawa_attractor import maths as AZ
from ValaQuenta.modules.zeta_uns_path import maths as ZM
which = sys.argv[1:] or ['A', 'B', 'C', 'D']
OUT = root / 'experiments/attractor_zeta_coupling.json'; R = json.load(open(OUT)) if OUT.exists() else {}
P = AZ.STANDARD; d = P['d']; TW = 2 * math.pi

def primes_to(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return s

if 'A' in which:
    h = 0.0045
    # f = 0: limit cycle of the (r, z) system
    p0 = dict(P, f=0.0); s = np.array([0.1, 0.0]); step = lambda s: (lambda k1: (lambda k2: (lambda k3: (lambda k4: s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4))(AZ.axisymmetric_rhs(s + h * k3, p0)))(AZ.axisymmetric_rhs(s + h / 2 * k2, p0)))(AZ.axisymmetric_rhs(s + h / 2 * k1, p0)))(AZ.axisymmetric_rhs(s, p0))
    for _ in range(40000): s = step(s)
    zs = []; ts = []
    for i in range(80000):
        s_new = step(s)
        if s[1] > 0 and AZ.axisymmetric_rhs(s, p0)[1] > 0 and AZ.axisymmetric_rhs(s_new, p0)[1] <= 0: ts.append(i * h); zs.append(s_new[1])
        s = s_new
    per = np.diff(ts); T_rz = float(per.mean()); R['A_f0'] = dict(periods=len(per), T_rz=T_rz, T_rz_sd=float(per.std()), z_max=float(np.mean(zs)), z_max_sd=float(np.std(zs)), w_rz=TW / T_rz, d=d, rotation_number=d * T_rz / TW, frac_part=(d * T_rz / TW) % 1)
    # f = 0.1: chaotic; return map at the z maxima, and spectra
    o = AZ.orbit(120000, 3000, h); x, y, z = o[:, 0], o[:, 1], o[:, 2]; rr = np.hypot(x, y); phi = np.unwrap(np.arctan2(y, x))
    idx = [i for i in range(1, len(z) - 1) if z[i] > z[i - 1] and z[i] >= z[i + 1] and z[i] > 0.9]
    dphi = np.diff(phi[idx]); Tm = np.diff(np.array(idx)) * h
    fft = lambda v: (np.fft.rfftfreq(len(v), h) * TW, np.abs(np.fft.rfft(v - v.mean())))
    wz, az = fft(z); wr, ar = fft(rr); top = lambda w, a, k=5: [round(float(w[j]), 3) for j in np.argsort(-a)[:k] if w[j] > 0.05]
    R['A_f01'] = dict(n_cycles=len(idx), mean_T_between_zmax=float(Tm.mean()), sd_T=float(Tm.std()), mean_dphi_per_cycle=float(dphi.mean()), dphi_mod_2pi_mean=float(np.mean(np.mod(dphi, TW))), dphi_mod_2pi_sd=float(np.std(np.mod(dphi, TW))),
                      z_spectrum_peaks=top(wz, az, 8), r_spectrum_peaks=top(wr, ar, 8), w_z=TW / float(Tm.mean()), d=d, ratio_d_over_wz=d / (TW / float(Tm.mean())))
    w1, w2 = d, R['A_f0']['w_rz']; C2 = 0.6601618158468696; B2 = 1.9021605831; consts = {'C2': C2, '2C2': 2 * C2, 'C2/2': C2 / 2, 'B2': B2, '1/C2': 1 / C2, 'ln2': math.log(2), 'ln3': math.log(3), '1/2': 0.5, 'pi/e': math.pi / math.e}
    rats = {'d/w_rz': w1 / w2, 'w_rz/d': w2 / w1, 'd/w_z(f=.1)': R['A_f01']['ratio_d_over_wz'], 'w_z/d': 1 / R['A_f01']['ratio_d_over_wz'], 'frac(d/w_rz)': (w1 / w2) % 1}
    hits = [(a, round(v, 5), b, round(c, 5)) for a, v in rats.items() for b, c in consts.items() if abs(v - c) / c < 0.01]
    R['A_twin_constants'] = dict(ratios={k: float(v) for k, v in rats.items()}, constants_tried=len(consts), ratios_tried=len(rats), within_1pct=hits, expected_by_chance=len(consts) * len(rats) * 0.02 * 0.5)
    # twin primes: sample the chaotic orbit at times n * h2
    h2 = 0.0451; N = 60000; sub = int(round(h2 / h)); orb = AZ.orbit(N * sub, 3000, h)[::sub]; xs, ys, zs_ = orb[:N, 0], orb[:N, 1], orb[:N, 2]; sv = primes_to(N + 3); n = np.arange(N)
    isp = sv[:N]; twin = isp & (sv[2:N + 2] | (n - 2 >= 0) & np.concatenate([[False, False], sv[:N - 2]])); lone = isp & ~twin; twin = twin & (n > 20); lone = lone & (n > 20)
    ph = np.mod(np.arctan2(ys, xs), TW)
    def rayleigh(a):
        C = np.mean(np.exp(1j * a)); Rn = abs(C) * len(a); p = math.exp(math.sqrt(1 + 4 * len(a) + 4 * (len(a) ** 2 - Rn ** 2)) - (1 + 2 * len(a))); return float(abs(C)), float(min(1.0, p))
    from scipy.stats import ks_2samp, mannwhitneyu
    t2 = {}
    for name, m in (('twin_primes', twin), ('lone_primes', lone), ('all_integers', (n > 20))): t2[name] = dict(n=int(m.sum()), rayleigh_R_p=rayleigh(ph[m]), rayleigh_2phi_R_p=rayleigh(2 * ph[m]))
    t2['phase_KS_twin_vs_lone'] = [float(v) for v in ks_2samp(ph[twin], ph[lone])]; t2['z_MWU_twin_vs_lone_p'] = float(mannwhitneyu(zs_[twin], zs_[lone]).pvalue); t2['z_KS_twin_vs_lone'] = [float(v) for v in ks_2samp(zs_[twin], zs_[lone])]
    t2['sampling'] = dict(h2=h2, alpha=d * h2 % TW, note='alpha = rotation per sample mod 2 pi: irrational-looking, so integer-indexed samples equidistribute for any subsequence of density > 0 (Weyl/Vinogradov)')
    R['A_twin_sampling'] = t2
    print(json.dumps({k: R[k] for k in ('A_f0', 'A_f01', 'A_twin_constants', 'A_twin_sampling')}, indent=0, default=float)[:3800])

if 'B' in which:
    from mpmath import mp, zeta, mpc; mp.dps = 15; K = 19.0
    o = AZ.orbit(400000, 3000, 0.0045)[::20]; x, y, z = o[:, 0], o[:, 1], o[:, 2]; s = x + 1j * K * y
    zg = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()][:6]); zer = np.array([0.5 + 1j * g for g in zg] + [0.5 - 1j * g for g in zg])
    dmin = np.min(np.abs(s[:, None] - zer[None, :]), axis=0); R['B_orbit_vs_zeros'] = dict(K=K, points=len(s), x_range=[float(x.min()), float(x.max())], t_range=[float(K * y.min()), float(K * y.max())], zero_in_range=[float(g) for g in zg if abs(g) <= K * max(abs(y.min()), y.max())],
        min_distance_s_plane_to_first_zeros=[round(float(v), 3) for v in dmin])
    cr = np.where(np.diff(np.sign(x - 0.5)) != 0)[0]; th = K * y[cr]; R['B_crossings_x_half'] = dict(n=int(len(cr)), t_min=float(th.min()), t_max=float(th.max()), frac_with_t_positive=float(np.mean(th > 0)))
    near = np.array([np.min(np.abs(np.abs(t) - zg)) for t in th if abs(t) < 22]); zz = np.array([14.134725, 21.022040, 25.010858]); rng = np.random.default_rng(4); inr = th[np.abs(th) < 26]
    nd = lambda v: np.min(np.abs(np.abs(v)[:, None] - zz[None, :]), axis=1); null = [float(nd(rng.uniform(-26, 26, len(inr))).mean()) for _ in range(500)]
    R['B_crossings_x_half'].update(dict(n_within_26=int(len(inr)), mean_dist_to_nearest_of_first3_zero_heights=float(nd(inr).mean()), null_mean=float(np.mean(null)), null_sd=float(np.std(null)), z_score=float((nd(inr).mean() - np.mean(null)) / np.std(null))))
    sub = np.arange(0, len(s), 8); vals = []
    for i in sub: vals.append(complex(zeta(mpc(float(x[i]), K * float(y[i])))))
    vals = np.array(vals); lg = np.log10(np.abs(vals) + 1e-9)
    R['B_paint'] = dict(points=len(sub), logabs_min=float(lg.min()), logabs_max=float(lg.max()), corr_logabs_with_z=float(np.corrcoef(lg, z[sub])[0, 1]), corr_logabs_with_r=float(np.corrcoef(lg, np.hypot(x[sub], y[sub]))[0, 1]))
    pos = np.stack([x[sub], z[sub], y[sub]], 1).astype(np.float32); col = np.clip((lg + 1.5) / 3.0, 0, 1).astype(np.float32); b64 = lambda a: base64.b64encode(np.asarray(a, np.float32).tobytes()).decode()
    R['_paint_blob'] = dict(pos=b64(pos.ravel()), col=b64(col), n=int(len(sub)), K=K, zeros=[float(g) for g in zg[:2]])
    print(json.dumps({k: R[k] for k in ('B_orbit_vs_zeros', 'B_crossings_x_half', 'B_paint')}, indent=0, default=float))

if 'C' in which:
    from scipy.optimize import brentq
    tstar = brentq(lambda t: float(ZM.theta_prime(t)) - d, 100, 50000); g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()]); nz = int(np.sum(g < tstar))
    R['C_theta_prime_equals_d'] = dict(t_star=tstar, two_pi_e7=2 * math.pi * math.e ** 7, zeros_below=nz, step=int(math.ceil(nz / 100)), block_zeros=[int((math.ceil(nz / 100) - 1) * 100 + 1), int(math.ceil(nz / 100) * 100)])
    rows = []
    for tt in (14.13, 236.5, 1000.0, tstar, 21942.0):
        dd = float(ZM.theta_prime(tt)); pp = dict(P, d=dd); lam = AZ.lyapunov_max(steps=30000, p=pp); o = AZ.orbit(30000, 3000, 0.0045, pp); rows.append(dict(t=tt, d=dd, lyapunov=float(lam), r_mean=float(np.hypot(o[:, 0], o[:, 1]).mean()), z_mean=float(o[:, 2].mean())))
    R['C_attractor_with_zeta_rotation_rate'] = rows; print(json.dumps({k: R[k] for k in ('C_theta_prime_equals_d', 'C_attractor_with_zeta_rotation_rate')}, indent=0, default=float))

if 'D' in which:
    from mpmath import mp, zeta, mpc; mp.dps = 15
    vert = {}; ts_ = np.linspace(0.1, 60.0, 600)
    for sg in (0.3, 0.4, 0.5, 0.6, 0.7, 1.0): vert[sg] = np.array([complex(zeta(mpc(sg, t))) for t in ts_])
    hor = {}; sx = np.linspace(-3.0, 4.0, 701)
    for c in (0.5, 1.0, math.pi / 2, math.pi, 2 * math.pi):
        vals = []
        for sgm in sx:
            if abs(sgm - 1.0) < 1e-9 and c == 0: continue
            vals.append(complex(zeta(mpc(float(sgm), c))))
        hor[c] = np.array(vals)
    R['D_vertical_lines_min_abs_zeta_t0_60'] = {str(k): float(np.abs(v).min()) for k, v in vert.items()}; R['D_horizontal_lines_min_abs_zeta_sigma_-3_4'] = {f'{k:.4f}': float(np.abs(v).min()) for k, v in hor.items()}
    pack = lambda dct: {f'{k:.4f}': [[round(v.real, 4), round(v.imag, 4)] for v in arr[::3]] for k, arr in dct.items()}
    R['_shadows'] = dict(vertical=pack(vert), horizontal=pack(hor)); print(json.dumps({k: R[k] for k in ('D_vertical_lines_min_abs_zeta_t0_60', 'D_horizontal_lines_min_abs_zeta_sigma_-3_4')}, indent=0, default=float))

js = {k: R[k] for k in ('_paint_blob', '_shadows') if k in R}
if js:
    with open(root / 'visualizations/javascript/attractor_zeta_data.js', 'w') as f: f.write('// generated by experiments/attractor_zeta_coupling.py (mpmath zeta along the Aizawa orbit read as s = x + 19 i y; images of the lines sigma = const and Im s = const). Not hand-edited.\nwindow.ATTR_ZETA = ' + json.dumps({'paint': js.get('_paint_blob'), 'shadows': js.get('_shadows')}) + ';\n')
json.dump({k: v for k, v in R.items() if not k.startswith('_')}, open(OUT, 'w'), indent=1, default=float)
