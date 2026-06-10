#!/usr/bin/env python3
import argparse, os, subprocess, tempfile, sys
import numpy as np
from multiprocessing import Pool

ZBIN = os.path.join(os.path.dirname(__file__), "..", "tools", "zhunt_c")

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

def _work(job):
    scaf, abs_start, sub, nstart, win, mindin, maxdin, thr = job
    with tempfile.NamedTemporaryFile("w", suffix=".seq", delete=False) as tf:
        tf.write(sub)
        tf.write("\n")
        path = tf.name
    try:
        subprocess.run([os.path.abspath(ZBIN), str(win), str(mindin), str(maxdin), path],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        zpath = path + ".Z-SCORE"
        z = np.empty(nstart, dtype=np.float64)
        with open(zpath) as fh:
            next(fh)
            for i in range(nstart):
                z[i] = float(fh.readline().split()[2])
        hits = np.flatnonzero(z > thr)
        return scaf, abs_start, (abs_start + hits).astype(np.int64), z[hits].astype(np.float32)
    finally:
        for p in (path, path + ".Z-SCORE"):
            try:
                os.remove(p)
            except OSError:
                pass

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default="data/genome.fna")
    ap.add_argument("--out", default="results/zhunt.bed")
    ap.add_argument("--mindin", type=int, default=6)
    ap.add_argument("--maxdin", type=int, default=12)
    ap.add_argument("--threshold", type=float, default=400.0)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--chunk", type=int, default=2_000_000)
    a = ap.parse_args()

    pad = 2 * a.maxdin
    win = a.maxdin
    jobs = []
    for name, seq in read_fasta(a.genome):
        L = len(seq)
        pos = 0
        while pos < L - pad:
            core = min(a.chunk, L - pad - pos)
            sub = seq[pos: pos + core + pad]
            jobs.append((name, pos, sub, core, win, a.mindin, a.maxdin, a.threshold))
            pos += core
    sys.stderr.write(f"zhunt_c: {len(jobs)} chunks, {a.procs} procs, "
                     f"win {win} din {a.mindin}-{a.maxdin}, thr>{a.threshold}\n")
    by_scaf = {}
    done = 0
    with Pool(a.procs) as p:
        for scaf, _, pos, zz in p.imap_unordered(_work, jobs):
            done += 1
            sys.stderr.write(f"  {done}/{len(jobs)}\r")
            if pos.size:
                by_scaf.setdefault(scaf, []).append((pos, zz))
    sys.stderr.write("\n")

    n_regions = 0
    with open(a.out, "w") as o:
        for scaf in sorted(by_scaf):
            allpos = np.concatenate([p for p, _ in by_scaf[scaf]])
            allz = np.concatenate([z for _, z in by_scaf[scaf]])
            order = np.argsort(allpos)
            allpos, allz = allpos[order], allz[order]
            start = prev = allpos[0]; maxz = allz[0]
            for pp, zz in zip(allpos[1:], allz[1:]):
                if pp <= prev + 1:
                    prev = pp; maxz = max(maxz, zz)
                else:
                    o.write(f"{scaf}\t{start}\t{prev+1}\tZ\t{maxz:.1f}\t.\n")
                    n_regions += 1
                    start = prev = pp; maxz = zz
            o.write(f"{scaf}\t{start}\t{prev+1}\tZ\t{maxz:.1f}\t.\n")
            n_regions += 1
    print(f"Z-DNA regions (z>{a.threshold}) from canonical Z-Hunt C: {n_regions} -> {a.out}")

if __name__ == "__main__":
    main()
