"""ζ'(ρ_k) at the 25,000 tabulated zeros, from Hardy's Z:  ζ'(ρ) = −i e^{−iθ(γ)} Z'(γ)   (ζ(½+it) = e^{−iθ} Z(t), ρ a zero).
Chunks of 250 zeros are written to experiments/zetaprime25k/chunk_NNNN.tsv (restartable: finished chunks are skipped); merged into data/zetaprime_25000.tsv.
Run:  nohup setsid python3 -u experiments/zetaprime_at_zeros.py > experiments/zetaprime25k/run.log 2>&1 &
Checked against mpmath's direct zeta(s, derivative=1) at zeros 1, 2: agreement 4e−17."""
import sys, os, pathlib, time
import mpmath as mp
from multiprocessing import Pool
ROOT = pathlib.Path(__file__).resolve().parents[1]; OUT = ROOT / "experiments/zetaprime25k"; CH = 250
G = [mp.mpf(l.split()[-1]) for l in open(ROOT / "data/zeros_25000_hp.tsv") if l[0].isdigit()]

def work(c):
    f = OUT / f"chunk_{c:04d}.tsv"
    if f.exists(): return c
    mp.mp.dps = 20; rows = []
    for k in range(c * CH, min((c + 1) * CH, len(G))):
        g = G[k]; zp = -1j * mp.e ** (-1j * mp.siegeltheta(g)) * mp.siegelz(g, derivative=1)
        rows.append(f"{k+1}\t{mp.nstr(g, 25)}\t{mp.nstr(zp.real, 17)}\t{mp.nstr(zp.imag, 17)}")
    f.write_text("\n".join(rows) + "\n"); return c

if __name__ == "__main__":
    t0 = time.time(); nch = (len(G) + CH - 1) // CH
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as p:
        for c in p.imap_unordered(work, range(nch)): print(f"chunk {c} done {time.time()-t0:.0f}s", flush=True)
    with open(ROOT / "data/zetaprime_25000.tsv", "w") as out:
        out.write("k\tgamma\tRe_zeta_prime\tIm_zeta_prime  (zeta'(rho) = -i exp(-i theta) Z'(gamma); dps 20)\n")
        for c in range(nch): out.write((OUT / f"chunk_{c:04d}.tsv").read_text())
    print("merged", time.time() - t0)
