#!/usr/bin/env python3
"""Map HMMER epigenetic-family hits to genes and build the family->gene table.

Inputs:
  data/annotation.gff           NCBI GFF3
  results/hmmer/epi_domtbl.txt  hmmsearch --domtblout (Pfam --cut_ga)
  hmm/families.tsv              pfam<TAB>name<TAB>category<TAB>description
Outputs:
  results/epigenetic_genes.tsv  one row per (family, gene)
  results/epigenetic_summary.tsv  one row per family with gene counts
"""
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

# ---- families metadata (hmm NAME -> (pfam, category, description)) ----
fam_meta = {}          # pfam acc -> row
name_to_pfam = {}      # hmm NAME -> pfam acc (filled from domtbl 'tname'/acc)
with open(FAMILIES) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        fam_meta[row["pfam"]] = row

# ---- parse GFF: gene coords by locus_tag; protein_id -> locus_tag ----
gene_info = {}             # locus_tag -> (scaffold, start, end, strand)
prot2locus = {}           # protein_id -> locus_tag
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

# ---- parse hmmsearch domtblout ----
# columns: tname tacc tlen qname qacc(pfam) qlen Eval(full) score bias ...
prot_fam = defaultdict(dict)   # protein -> {pfam: best_evalue}
fam_name = {}                  # pfam acc -> hmm name
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

# ---- build family -> gene rows (collapse isoforms to gene/locus_tag) ----
# best evalue per (pfam, locus_tag)
fam_gene = defaultdict(dict)   # pfam -> {locus_tag: (best_eval, protein_id)}
for prot, fams in prot_fam.items():
    lt = prot2locus.get(prot)
    if lt is None:
        # fallback: strip version? keep protein as pseudo-gene
        lt = prot
    for pfam, ev in fams.items():
        cur = fam_gene[pfam].get(lt)
        if cur is None or ev < cur[0]:
            fam_gene[pfam][lt] = (ev, prot)

# ---- write detailed table ----
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

# sort by category then family then evalue
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

# ---- summary per family ----
with open("results/epigenetic_summary.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["family_pfam", "family_name", "category", "description", "n_genes"])
    # include all curated families, even with 0 hits
    for pfam, meta in fam_meta.items():
        n = len(fam_gene.get(pfam, {}))
        w.writerow([pfam, meta["name"], meta["category"], meta["description"], n])

n_genes_total = len({lt for genes in fam_gene.values() for lt in genes})
n_fam_hit = sum(1 for pfam in fam_meta if fam_gene.get(pfam))
print(f"Families searched: {len(fam_meta)} | with >=1 gene: {n_fam_hit}")
print(f"Distinct epigenetic genes identified: {n_genes_total}")
print(f"Detailed rows: {len(rows)} -> results/epigenetic_genes.tsv")
print(f"Summary       -> results/epigenetic_summary.tsv")
