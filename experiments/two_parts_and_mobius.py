"""Two parts of the zero counting function (discrete = zeros, smooth = theta'/pi) and the Mobius-band reading of the path.  Writes two_parts_and_mobius.json.
Expectations: (1) the smooth part carries NO prime-power lines: its spectrum in x = ln n is flat at ln p^k; the discrete part has lines at ln p^k of height -Lambda(n)/sqrt(n) (Landau);
(2) the fluctuation S(T) is uncorrelated with the smooth part N_smooth(T) and its slope; (3) Z changes sign at every zero, and theta advances by ~pi between consecutive zeros, so the path is a section of the tautological line bundle over RP^1 (a Mobius band): core crossings = zeros, half-twists = theta/pi."""
import json, math, pathlib, sys, numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.zeta_uns_path import maths as ZM
g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()]); T = g[-1]; R = {}
def Lam(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0
ns = np.arange(2, 61); x = np.log(ns); lam = np.array([Lam(int(n)) for n in ns])
# discrete part: sum over zeros of cos(gamma x) / (T/2pi)
disc = np.array([np.sum(np.cos(g * xx)) for xx in x]) / (T / (2 * math.pi))
# smooth part: integral of (theta'(t)/pi) cos(t x) dt over [0,T], by quadrature
tt = np.linspace(14.0, T, 4_000_001); rho = ZM.theta_prime(tt) / math.pi
smooth = np.array([np.trapezoid(rho * np.cos(tt * xx), tt) for xx in x]) / (T / (2 * math.pi))
pred = -lam / np.sqrt(ns)
pp = lam > 0
R['spectra'] = dict(corr_discrete_vs_minus_Lambda_over_sqrt_n=float(np.corrcoef(disc, pred)[0, 1]), corr_smooth_vs_minus_Lambda_over_sqrt_n=float(np.corrcoef(smooth, pred)[0, 1]),
                    rms_discrete_at_prime_powers=float(np.sqrt(np.mean(disc[pp] ** 2))), rms_discrete_at_other_n=float(np.sqrt(np.mean(disc[~pp] ** 2))),
                    rms_smooth_at_prime_powers=float(np.sqrt(np.mean(smooth[pp] ** 2))), rms_smooth_at_other_n=float(np.sqrt(np.mean(smooth[~pp] ** 2))),
                    mean_discrete_at_prime_powers=float(disc[pp].mean()), mean_pred=float(pred[pp].mean()))
# (2) S(T) vs smooth
th = ZM.theta(g); k = np.arange(1, len(g) + 1); S = k - th / math.pi - 1 + 0.0
Nsm = th / math.pi + 1; R['S_vs_smooth'] = dict(corr_S_Nsmooth=float(np.corrcoef(S, Nsm)[0, 1]), corr_S_logslope=float(np.corrcoef(S, ZM.theta_prime(g))[0, 1]), S_mean=float(S.mean()), S_std=float(S.std()))
# (3) Mobius structure
d_theta = np.diff(th) / math.pi; tm = 0.5 * (g[:-1] + g[1:]); Zm = ZM.Z_rs(tm); Zq = ZM.Z_rs(g + 1e-3)
sgn_alt = float(np.mean(np.sign(Zm[1:]) != np.sign(Zm[:-1])))
R['mobius'] = dict(share_adjacent_lobes_with_opposite_sign_of_Z=sgn_alt, theta_advance_between_zeros_over_pi_mean=float(d_theta.mean()), p5_p95=[float(np.percentile(d_theta, 5)), float(np.percentile(d_theta, 95))], total_theta_over_pi_to_T=float(th[-1] / math.pi), zeros=len(g), half_twists_minus_zeros=float(th[-1] / math.pi - len(g)))
json.dump(R, open(root / 'experiments/two_parts_and_mobius.json', 'w'), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
