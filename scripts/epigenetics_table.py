import csv, re, sys
from collections import defaultdict

GFF = "data/annotation.gff"
DOMTBL = "results/hmmer/epi_domtbl.txt"
FAMILIES = "hmm/families.tsv"

def attrs(col):
    d = {}
    for kv in col.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d

fam_meta = {}
name_to_pfam = {}
with open(FAMILIES) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        fam_meta[row["pfam"]] = row

gene_info = {}
prot2locus = {}
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        ftype = f[2]
        a = attrs(f[8])
        if ftype == "gene":
            lt = a.get("locus_tag") or a.get("Name")
            if lt:
                gene_info[lt] = (f[0], int(f[3]), int(f[4]), f[6])
        elif ftype == "CDS":
            pid = a.get("protein_id")
            lt = a.get("locus_tag")
            if pid and lt:
                prot2locus[pid] = lt

prot_fam = defaultdict(dict)
fam_name = {}
with open(DOMTBL) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        c = line.split()
        if len(c) < 13:
            continue
        prot = c[0]
        qname = c[3]
        qacc = c[4]
        pfam = qacc.split(".")[0] if qacc.startswith("PF") else qacc
        full_eval = float(c[6])
        fam_name[pfam] = qname
        prev = prot_fam[prot].get(pfam)
        if prev is None or full_eval < prev:
            prot_fam[prot][pfam] = full_eval

fam_gene = defaultdict(dict)
for prot, fams in prot_fam.items():
    lt = prot2locus.get(prot)
    if lt is None:
        lt = prot
    for pfam, ev in fams.items():
        cur = fam_gene[pfam].get(lt)
        if cur is None or ev < cur[0]:
            fam_gene[pfam][lt] = (ev, prot)

import os
os.makedirs("results", exist_ok=True)
rows = []
for pfam, genes in fam_gene.items():
    meta = fam_meta.get(pfam, {})
    for lt, (ev, prot) in genes.items():
        gi = gene_info.get(lt)
        if gi:
            scaf, s, e, strand = gi
            coord = f"{scaf}:{s}-{e}({strand})"
        else:
            coord = "NA"
        rows.append({
            "family_pfam": pfam,
            "family_name": fam_name.get(pfam, meta.get("name", "")),
            "category": meta.get("category", ""),
            "description": meta.get("description", ""),
            "gene_locus_tag": lt,
            "coordinates": coord,
            "best_protein": prot,
            "best_evalue": f"{ev:.1e}",
        })

cat_order = {c: i for i, c in enumerate([
    "DNA methylation", "DNA methylation read", "DNA modification",
    "Histone fold", "Histone write (methyl)", "Histone write (acetyl)",
    "Histone erase (deacetyl)", "Histone erase (demethyl)", "Histone read",
    "Chromatin remodeling", "Cohesin/SMC", "Cohesin"])}
rows.sort(key=lambda r: (cat_order.get(r["category"], 99), r["family_pfam"], r["best_evalue"]))

with open("results/epigenetic_genes.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

with open("results/epigenetic_summary.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["family_pfam", "family_name", "category", "description", "n_genes"])
    for pfam, meta in fam_meta.items():
        n = len(fam_gene.get(pfam, {}))
        w.writerow([pfam, meta["name"], meta["category"], meta["description"], n])

n_genes_total = len({lt for genes in fam_gene.values() for lt in genes})
n_fam_hit = sum(1 for pfam in fam_meta if fam_gene.get(pfam))
print(f"Families searched: {len(fam_meta)} | with >=1 gene: {n_fam_hit}")
print(f"Distinct epigenetic genes identified: {n_genes_total}")
print(f"Detailed rows: {len(rows)} -> results/epigenetic_genes.tsv")
print(f"Summary       -> results/epigenetic_summary.tsv")
