#!/usr/bin/env python3
import sys, os
import numpy as np

GFF = "data/annotation.gff"
GENOME = "data/genome.fna"
PROMOTER = 1000
DOWNSTREAM = 200
CATS = ["Exon", "Intron", "Promoter", "Downstream", "Intergenic"]
EXON, INTRON, PROMOTER_C, DOWNSTREAM_C, INTERGENIC = 0, 1, 2, 3, 4

def attrs(col):
    d = {}
    for kv in col.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d

def scaffold_lengths(path):
    lens, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    lens[name] = n
                name = line[1:].split()[0]; n = 0
            else:
                n += len(line.strip())
    if name:
        lens[name] = n
    return lens

def build_category_arrays(lens):
    genes = []
    exons = []
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t = f[2]
            if t == "gene":
                genes.append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
            elif t == "exon":
                exons.append((f[0], int(f[3]) - 1, int(f[4])))
    arr = {s: np.full(L, INTERGENIC, dtype=np.int8) for s, L in lens.items()}

    def paint(scaf, s, e, code, only_intergenic=False):
        L = lens.get(scaf)
        if L is None:
            return
        s = max(0, s); e = min(L, e)
        if s >= e:
            return
        if only_intergenic:
            seg = arr[scaf][s:e]
            seg[seg == INTERGENIC] = code
        else:
            arr[scaf][s:e] = code

    for scaf, s, e, strand in genes:
        if strand == "+":
            paint(scaf, e, e + DOWNSTREAM, DOWNSTREAM_C, only_intergenic=True)
        else:
            paint(scaf, s - DOWNSTREAM, s, DOWNSTREAM_C, only_intergenic=True)
    for scaf, s, e, strand in genes:
        ss, ee = (s - PROMOTER, s) if strand == "+" else (e, e + PROMOTER)
        L = lens[scaf]; ss = max(0, ss); ee = min(L, ee)
        if ss < ee:
            seg = arr[scaf][ss:ee]
            seg[(seg == INTERGENIC) | (seg == DOWNSTREAM_C)] = PROMOTER_C
    for scaf, s, e, strand in genes:
        paint(scaf, s, e, INTRON)
    for scaf, s, e in exons:
        paint(scaf, s, e, EXON)
    return arr

def load_bed(path):
    iv = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return iv

def distribute(label, bed_path, arr, lens, out_dir):
    iv = load_bed(bed_path)
    total = sum(len(v) for v in iv.values())
    assigned = np.zeros(len(CATS), dtype=np.int64)
    overlap = np.zeros(len(CATS), dtype=np.int64)
    for scaf, ivs in iv.items():
        a = arr.get(scaf)
        if a is None:
            continue
        for s, e in ivs:
            seg = a[s:e]
            if seg.size == 0:
                continue
            codes = np.unique(seg)
            assigned[codes.min()] += 1
            for c in codes:
                overlap[c] += 1
    bg = np.zeros(len(CATS), dtype=np.int64)
    for s, a in arr.items():
        u, c = np.unique(a, return_counts=True)
        for code, cnt in zip(u, c):
            bg[code] += cnt
    genome_bp = bg.sum()

    region_total = np.zeros(len(CATS), dtype=np.int64)
    region_hit = np.zeros(len(CATS), dtype=np.int64)
    for scaf, a in arr.items():
        mask = np.zeros(a.shape[0], dtype=bool)
        for s, e in iv.get(scaf, []):
            if e > mask.size:
                e = mask.size
            if s < mask.size:
                mask[s:e] = True
        change = np.empty(a.shape[0], dtype=bool)
        change[0] = True
        change[1:] = a[1:] != a[:-1]
        run_starts = np.flatnonzero(change)
        run_ends = np.empty_like(run_starts)
        run_ends[:-1] = run_starts[1:]
        run_ends[-1] = a.shape[0]
        run_codes = a[run_starts]
        csum = np.concatenate(([0], np.cumsum(mask)))
        run_has = (csum[run_ends] - csum[run_starts]) > 0
        for code in range(len(CATS)):
            sel = run_codes == code
            region_total[code] += int(sel.sum())
            region_hit[code] += int((sel & run_has).sum())

    os.makedirs(out_dir, exist_ok=True)
    t1 = os.path.join(out_dir, f"table1_{label}.tsv")
    with open(t1, "w") as o:
        o.write("Region\tN_structures\tFraction\tBackground_bp_fraction\tEnrichment\n")
        for i, cat in enumerate(CATS):
            frac = assigned[i] / total if total else 0
            bgf = bg[i] / genome_bp if genome_bp else 0
            enr = (frac / bgf) if bgf > 0 else float("nan")
            o.write(f"{cat}\t{assigned[i]}\t{frac:.4f}\t{bgf:.4f}\t{enr:.2f}\n")
        o.write(f"TOTAL\t{total}\t1.0000\t1.0000\t-\n")
    t2 = os.path.join(out_dir, f"table2_{label}.tsv")
    with open(t2, "w") as o:
        o.write("Region\tN_regions\tN_regions_with_structure\tFraction_regions_with_structure\n")
        for i, cat in enumerate(CATS):
            fr = region_hit[i] / region_total[i] if region_total[i] else 0
            o.write(f"{cat}\t{region_total[i]}\t{region_hit[i]}\t{fr:.4f}\n")
    return dict(label=label, total=total, assigned=assigned.copy(), overlap=overlap.copy(),
                bg=bg.copy(), genome_bp=genome_bp, region_total=region_total.copy(),
                region_hit=region_hit.copy(), t1=t1, t2=t2)

def main():
    if len(sys.argv) < 2:
        print("usage: genomic_distribution.py LABEL1=file1.bed [LABEL2=file2.bed ...]"); sys.exit(1)
    items = []
    for a in sys.argv[1:]:
        label, path = a.split("=", 1)
        items.append((label, path))
    lens = scaffold_lengths(GENOME)
    sys.stderr.write("Building genomic feature partition...\n")
    arr = build_category_arrays(lens)
    bg = np.zeros(len(CATS), dtype=np.int64)
    for s, a in arr.items():
        u, c = np.unique(a, return_counts=True)
        for code, cnt in zip(u, c):
            bg[code] += cnt
    gbp = bg.sum()
    sys.stderr.write("Genomic bp by category:\n")
    for i, cat in enumerate(CATS):
        sys.stderr.write(f"  {cat:12s} {bg[i]:>12,} ({bg[i]/gbp*100:5.2f}%)\n")
    for label, path in items:
        r = distribute(label, path, arr, lens, "results/distribution")
        print(f"\n=== {label}: {r['total']} structures ===")
        print("Region        N      Frac   BgFrac  Enrich")
        for i, cat in enumerate(CATS):
            frac = r["assigned"][i] / r["total"] if r["total"] else 0
            bgf = r["bg"][i] / r["genome_bp"]
            enr = frac / bgf if bgf > 0 else float("nan")
            print(f"  {cat:12s} {r['assigned'][i]:>5} {frac:7.3f} {bgf:7.3f} {enr:6.2f}")

if __name__ == "__main__":
    main()
