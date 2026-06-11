import argparse, os, sys
import numpy as np
from multiprocessing import Pool

RT = 0.59004
A0 = 0.357
B = 0.4
A_HALF = A0 / 2.0
K_RT = -0.2521201
SIGMA = 16.94800353
EXPLIMIT = -600.0
AVG = 29.6537135
STDV = 2.71997
_SQRT2 = 0.70710678118654752440
_SQRTPI = 0.564189583546

DBZED = np.array([
 [4.40,6.20,3.40,5.20,2.50,4.40,1.40,3.30,3.30,5.20,2.40,4.20,1.40,3.40,0.66,2.40,4.26],
 [4.40,2.50,3.30,1.40,6.20,4.40,5.20,3.40,3.40,1.40,2.40,0.66,5.20,3.30,4.20,2.40,4.26],
 [6.20,6.20,5.20,5.20,6.20,6.20,5.20,5.20,5.20,5.20,4.00,4.00,5.20,5.20,4.00,4.00,4.26],
 [6.20,6.20,5.20,5.20,6.20,6.20,5.20,5.20,5.20,5.20,4.00,4.00,5.20,5.20,4.00,4.00,4.26],
], dtype=np.float64)
EXPDBZED = np.exp(-DBZED / RT)

def bztwist_arr(n):
    return 0.8 + A0 * np.arange(1, n + 1, dtype=np.float64)

_LUT = np.full(256, 4, dtype=np.int8)
for _ch, _co in (("a", 0), ("t", 1), ("g", 2), ("c", 3)):
    _LUT[ord(_ch)] = _co
    _LUT[ord(_ch.upper())] = _co

def seq_to_dinuc_index(seq_bytes):
    codes = _LUT[np.frombuffer(seq_bytes, dtype=np.uint8)].astype(np.int64)
    b1 = codes[:-1]; b2 = codes[1:]
    idx = b1 * 4 + b2
    idx[(b1 >= 4) | (b2 >= 4)] = 16
    return idx

def assign_probability_vec(dl):
    z = np.abs(dl - AVG) / STDV
    x = z * _SQRT2
    y = _SQRTPI * np.exp(-x * x)
    z2 = z * z
    s = np.zeros_like(dl)
    xk = x.copy()
    k = 1.0
    for _ in range(600):
        s_new = s + xk
        k += 2.0
        xk = xk * (z2 / k)
        if not np.any(s_new > s):
            s = s_new
            break
        s = s_new
    tail = 0.5 - y * s
    return np.where(dl > AVG, tail, 1.0 / tail)

def zhunt_scores(dinuc_idx, nstart, mindin, maxdin):
    bestdl = np.full(nstart, 50.0, dtype=np.float64)
    for din in range(mindin, maxdin + 1):
        idx0 = dinuc_idx[0:nstart]
        vAS = DBZED[0][idx0].copy()
        vSA = DBZED[1][idx0].copy()
        bp_to_AS = np.empty((din, nstart), dtype=np.int8)
        bp_to_SA = np.empty((din, nstart), dtype=np.int8)
        for k in range(1, din):
            idxk = dinuc_idx[2 * k: 2 * k + nstart]
            cAS_fromAS = vAS + DBZED[0][idxk]
            cAS_fromSA = vSA + DBZED[3][idxk]
            from_AS = cAS_fromAS <= cAS_fromSA
            nvAS = np.where(from_AS, cAS_fromAS, cAS_fromSA)
            bp_to_AS[k] = np.where(from_AS, 0, 1)
            cSA_fromSA = vSA + DBZED[1][idxk]
            cSA_fromAS = vAS + DBZED[2][idxk]
            from_SA = cSA_fromSA <= cSA_fromAS
            nvSA = np.where(from_SA, cSA_fromSA, cSA_fromAS)
            bp_to_SA[k] = np.where(from_SA, 1, 0)
            vAS, vSA = nvAS, nvSA
        states = np.empty((din, nstart), dtype=np.int8)
        cur = np.where(vAS <= vSA, 0, 1).astype(np.int8)
        states[din - 1] = cur
        for k in range(din - 1, 0, -1):
            prev = np.where(cur == 0, bp_to_AS[k], bp_to_SA[k])
            states[k - 1] = prev
            cur = prev
        BE = np.empty((din, nstart), dtype=np.float64)
        BE[0] = EXPDBZED[np.where(states[0] == 0, 0, 1), dinuc_idx[0:nstart]]
        for k in range(1, din):
            idxk = dinuc_idx[2 * k: 2 * k + nstart]
            curk = states[k]; prevk = states[k - 1]
            row = np.where(curk == 0,
                           np.where(prevk == 0, 0, 3),
                           np.where(prevk == 1, 1, 2)).astype(np.int8)
            BE[k] = EXPDBZED[row, idxk]
        logcoef = np.empty((din, nstart), dtype=np.float64)
        prod = np.ones((din, nstart), dtype=np.float64)
        for i in range(din):
            sumv = np.zeros(nstart, dtype=np.float64)
            for j in range(din - i):
                prod[j] *= BE[i + j]
                sumv += prod[j]
            logcoef[i] = np.log(sumv)
        bzt = bztwist_arr(din)[:, None]
        deltatwist = A_HALF * din

        def delta_linking(dl_arr):
            z = dl_arr[None, :] - bzt
            expo = logcoef + K_RT * z * z
            mn = expo.min(axis=0)
            expmini = np.where(mn < EXPLIMIT, EXPLIMIT - mn, 0.0)
            E = np.exp(expo + expmini[None, :])
            sumq = E.sum(axis=0) + np.exp(K_RT * dl_arr * dl_arr + SIGMA + expmini)
            sump = (bzt * E).sum(axis=0)
            return deltatwist - sump / sumq

        f1 = delta_linking(np.full(nstart, 10.0))
        f2 = delta_linking(np.full(nstart, 50.0))
        bracketed = (f1 * f2 < 0.0)
        x = np.where(f1 < 0.0, 10.0, 50.0)
        dx = np.where(f1 < 0.0, 40.0, -40.0)
        for _ in range(16):
            dx *= 0.5
            xmid = x + dx
            fmid = delta_linking(xmid)
            x = np.where(fmid <= 0.0, xmid, x)
        dl_din = np.where(bracketed, x, 50.0)
        np.minimum(bestdl, dl_din, out=bestdl)
    return assign_probability_vec(bestdl)

def zhunt_scalar_one(idx, s, din):
    NEG = float("inf")
    vAS = DBZED[0][idx[s]]; vSA = DBZED[1][idx[s]]
    pAS = [0] * din; pSA = [0] * din
    for k in range(1, din):
        ik = idx[s + 2 * k]
        cAS = (vAS + DBZED[0][ik], 0) if (vAS + DBZED[0][ik]) <= (vSA + DBZED[3][ik]) else (vSA + DBZED[3][ik], 1)
        cSA = (vSA + DBZED[1][ik], 1) if (vSA + DBZED[1][ik]) <= (vAS + DBZED[2][ik]) else (vAS + DBZED[2][ik], 0)
        vAS, pAS[k] = cAS
        vSA, pSA[k] = cSA
    cur = 0 if vAS <= vSA else 1
    st = [0] * din; st[din - 1] = cur
    for k in range(din - 1, 0, -1):
        cur = pAS[k] if cur == 0 else pSA[k]
        st[k - 1] = cur
    BE = [0.0] * din
    BE[0] = EXPDBZED[0 if st[0] == 0 else 1][idx[s]]
    for k in range(1, din):
        ik = idx[s + 2 * k]
        if st[k] == 0:
            row = 0 if st[k - 1] == 0 else 3
        else:
            row = 1 if st[k - 1] == 1 else 2
        BE[k] = EXPDBZED[row][ik]
    logcoef = [0.0] * din
    prod = [1.0] * din
    for i in range(din):
        sm = 0.0
        for j in range(din - i):
            prod[j] *= BE[i + j]
            sm += prod[j]
        logcoef[i] = np.log(sm)
    bzt = bztwist_arr(din)
    deltatwist = A_HALF * din
    def dlf(dl):
        expo = [logcoef[i] + K_RT * (dl - bzt[i]) ** 2 for i in range(din)]
        mn = min(expo)
        expmini = (EXPLIMIT - mn) if mn < EXPLIMIT else 0.0
        E = [np.exp(e + expmini) for e in expo]
        sumq = sum(E) + np.exp(K_RT * dl * dl + SIGMA + expmini)
        sump = sum(bzt[i] * E[i] for i in range(din))
        return deltatwist - sump / sumq
    f1 = dlf(10.0); f2 = dlf(50.0)
    if f1 * f2 >= 0:
        return 50.0
    x = 10.0 if f1 < 0 else 50.0
    dx = 40.0 if f1 < 0 else -40.0
    while abs(dx) > 0.001:
        dx *= 0.5
        xm = x + dx
        if dlf(xm) <= 0.0:
            x = xm
    return x

def zhunt_scalar(idx, nstart, mindin, maxdin):
    out = np.empty(nstart)
    for s in range(nstart):
        best = 50.0
        for din in range(mindin, maxdin + 1):
            d = zhunt_scalar_one(idx, s, din)
            if d < best:
                best = d
        out[s] = best
    return assign_probability_vec(np.array([out]).ravel()) if False else out

def selftest():
    import random
    random.seed(12345)
    for trial in range(3):
        seq = "".join(random.choice("ACGT") for _ in range(400))
        seq = seq[:150] + "CG" * 20 + seq[190:]
        b = seq.encode()
        idx = seq_to_dinuc_index(b)
        mindin, maxdin = 6, 12
        nstart = len(seq) - 2 * maxdin
        zv = zhunt_scores(idx, nstart, mindin, maxdin)
        dl_scalar = zhunt_scalar(idx, nstart, mindin, maxdin)
        zs = assign_probability_vec(dl_scalar)
        diff = np.abs(zv - zs)
        rel = diff / np.maximum(np.abs(zs), 1e-9)
        print(f"trial {trial}: n={nstart} max|Δz|={diff.max():.3e} "
              f"max rel={rel.max():.3e} max z(vec)={zv.max():.1f} max z(scalar)={zs.max():.1f}")
        assert diff.max() < 1e-4, "vectorised and scalar disagree!"
    print("SELFTEST PASSED: vectorised engine matches scalar reference.")

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:].split()[0]; seq = []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)

_GLOBAL = {}

def _work(job):
    scaf, abs_start, sub_bytes, nstart, mindin, maxdin, thr = job
    idx = seq_to_dinuc_index(sub_bytes)
    z = zhunt_scores(idx, nstart, mindin, maxdin)
    hits = np.flatnonzero(z > thr)
    return scaf, abs_start, hits.astype(np.int64), z[hits].astype(np.float32)

def run_genome(genome, out_bed, mindin, maxdin, thr, procs, chunk, only_scaf=None):
    jobs = []
    seqs = {}
    for name, seq in read_fasta(genome):
        if only_scaf and name not in only_scaf:
            continue
        seqs[name] = seq
        L = len(seq)
        b = seq.encode()
        pad = 2 * maxdin
        pos = 0
        while pos < L - pad:
            core = min(chunk, L - pad - pos)
            sub = b[pos: pos + core + pad]
            nstart = core
            jobs.append((name, pos, sub, nstart, mindin, maxdin, thr))
            pos += core
    sys.stderr.write(f"Z-Hunt: {len(jobs)} chunks over {len(seqs)} scaffolds, "
                     f"{procs} procs, din {mindin}-{maxdin}, thr>{thr}\n")
    results = []
    with Pool(procs) as p:
        for i, r in enumerate(p.imap_unordered(_work, jobs)):
            results.append(r)
            sys.stderr.write(f"  chunk {i+1}/{len(jobs)} done\r")
    sys.stderr.write("\n")
    by_scaf = {}
    for scaf, abs_start, hits, zz in results:
        if hits.size == 0:
            continue
        by_scaf.setdefault(scaf, []).append((abs_start + hits, zz))
    n_regions = 0
    with open(out_bed, "w") as o:
        for scaf in sorted(by_scaf):
            allpos = np.concatenate([p for p, _ in by_scaf[scaf]])
            allz = np.concatenate([z for _, z in by_scaf[scaf]])
            order = np.argsort(allpos)
            allpos = allpos[order]; allz = allz[order]
            start = allpos[0]; prev = allpos[0]; maxz = allz[0]
            for p, z in zip(allpos[1:], allz[1:]):
                if p <= prev + 1:
                    prev = p; maxz = max(maxz, z)
                else:
                    o.write(f"{scaf}\t{start}\t{prev+1}\tZ\t{maxz:.1f}\t.\n")
                    n_regions += 1
                    start = p; prev = p; maxz = z
            o.write(f"{scaf}\t{start}\t{prev+1}\tZ\t{maxz:.1f}\t.\n")
            n_regions += 1
    print(f"Z-DNA regions (z>{thr}): {n_regions} -> {out_bed}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--genome", default="data/genome.fna")
    ap.add_argument("--out", default="results/zhunt.bed")
    ap.add_argument("--mindin", type=int, default=6)
    ap.add_argument("--maxdin", type=int, default=12)
    ap.add_argument("--threshold", type=float, default=400.0)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--scaffold", default=None, help="comma-separated scaffold ids to limit to")
    args = ap.parse_args()
    if args.selftest:
        selftest(); return
    only = set(args.scaffold.split(",")) if args.scaffold else None
    run_genome(args.genome, args.out, args.mindin, args.maxdin, args.threshold,
               args.procs, args.chunk, only)

if __name__ == "__main__":
    main()
