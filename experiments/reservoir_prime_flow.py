"""Can the Aizawa attractor carry a prime-hashed word stream?  A reservoir-memory test, the nearest thing to an orbweaver test that can be run today.
(v2: additive input injection s += 0.4 (cos a, sin a, gamma) as in standard reservoirs, and an 8-way task on the 8 most frequent words; v1 rotated the state by the kick and scored 16 classes: zero gain everywhere, a flawed design, kept in the log.)
Each word is hashed to a prime (VAPMIP prime_hash: spelling code, next_prime); its phase kick is ln(prime) mod 2 pi (the 'prime flow') and its windspeed (gamma_radial of its first WordNet sense, the scaled context scalar) kicks z.
The state is then advanced by the attractor (RK4, tau = 0.2).  A nonlinear readout (extra trees) tries to name the word seen m steps ago (15 most frequent words + 'other') from the state alone.
Expectations written before the run: memory falls with the lag m; a linear 3-D rotation reservoir of the same size does worse at m >= 1 than the chaotic attractor OR is no worse (reported either way); shuffled words give the majority rate (leak check).
No claim about the orbweaver is made beyond what this measures.  Writes reservoir_prime_flow.json."""
import collections, json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent)); sys.path.insert(0, str(root.parent / 'VAPMIP'))
from ValaQuenta.modules.aizawa_attractor import maths as AZ
import prime_hash as PH, wordnet_boxkite as WB
from nltk.corpus import brown, wordnet as wn
from sklearn.ensemble import ExtraTreesClassifier
N = 30000; words = [w.lower() for w in brown.words()[:60000] if w.isalpha()][:N]
cnt = collections.Counter(words); top = [w for w, _ in cnt.most_common(15)]; cls = np.array([top.index(w) if w in top else 15 for w in words])
cache = {}
def hashp(w):
    if w not in cache:
        c = PH.spell(w); a = PH.next_prime(c) if c > 1 else 2; g = -0.5
        ss = wn.synsets(w)
        if ss:
            v = WB.gamma_radial(ss[0]); g = -0.5 if v is None else float(v)
        cache[w] = (math.log(a) % (2 * math.pi), g)
    return cache[w]
kick = np.array([hashp(w) for w in words]); TW = 2 * math.pi
def rot(s, a): c, sn = math.cos(a), math.sin(a); return np.array([c * s[0] - sn * s[1], sn * s[0] + c * s[1], s[2]])
def run_attractor(ext=False, tau_steps=45, h=0.0045):
    s = np.array([0.1, 0.0, 0.0])
    for _ in range(3000): s = AZ.rk4_step(s, h)
    F = []
    for k in range(N):
        s = s + 0.4 * np.array([math.cos(kick[k, 0]), math.sin(kick[k, 0]), kick[k, 1]]); fs = []
        for i in range(tau_steps):
            s = AZ.rk4_step(s, h)
            if ext and i in (14, 29): fs += list(s)
        F.append(fs + list(s))
        if not np.all(np.isfinite(s)) or np.abs(s).max() > 20: s = np.array([0.1, 0.0, 0.0])
    return np.array(F)
def run_linear():
    s = np.array([0.1, 0.0, 0.0]); F = []; c, sn = math.cos(0.5), math.sin(0.5)
    for k in range(N):
        s = s + 0.4 * np.array([math.cos(kick[k, 0]), math.sin(kick[k, 0]), kick[k, 1]]); s = np.array([0.95 * (c * s[0] - sn * s[1]), 0.95 * (sn * s[0] + c * s[1]), 0.9 * s[2]]); F.append(list(s))
    return np.array(F)
def score(F, y, lags=range(0, 9), seed=0):
    out = {}
    for m in lags:
        X = F[m:]; t = y[:len(y) - m]; keep = t < 8; X = X[keep]; t = t[keep]; n = len(X); cut = int(0.7 * n)      # 8-way task: only positions where the word m steps ago is one of the 8 most frequent
        clf = ExtraTreesClassifier(n_estimators=60, min_samples_leaf=3, n_jobs=-1, random_state=seed).fit(X[:cut], t[:cut]); acc = float(np.mean(clf.predict(X[cut:]) == t[cut:])); maj = float(np.mean(t[cut:] == np.bincount(t[:cut]).argmax()))
        out[m] = dict(acc=round(acc, 4), majority=round(maj, 4), gain=round(acc - maj, 4))
    return out
R = {}
Fa = run_attractor(False); R['attractor_3'] = score(Fa, cls); print('attractor 3 features', {m: v['gain'] for m, v in R['attractor_3'].items()}, flush=True)
Fe = run_attractor(True); R['attractor_9'] = score(Fe, cls); print('attractor 9 features', {m: v['gain'] for m, v in R['attractor_9'].items()}, flush=True)
Fl = run_linear(); R['linear_3'] = score(Fl, cls); print('linear 3 features', {m: v['gain'] for m, v in R['linear_3'].items()}, flush=True)
rng = np.random.default_rng(0); R['shuffled_attractor_9_control'] = score(Fe, rng.permutation(cls), lags=[0, 1, 4]); print('shuffled control', R['shuffled_attractor_9_control'], flush=True)
R['setup'] = dict(words=N, classes=16, top15=top, majority_class_share=float(np.mean(cls == np.bincount(cls).argmax())), phase_kick_std=float(np.std(kick[:, 0])), n_unique_words=len(cache))
json.dump(R, open(root / 'experiments/reservoir_prime_flow.json', 'w'), indent=1, default=float)
