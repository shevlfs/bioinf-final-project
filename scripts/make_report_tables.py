#!/usr/bin/env python3
import csv, os

os.makedirs("results/markdown", exist_ok=True)

rows = list(csv.DictReader(open("results/epigenetic_summary.tsv"), delimiter="\t"))
best = {}
for r in csv.DictReader(open("results/epigenetic_genes.tsv"), delimiter="\t"):
    pf = r["family_pfam"]
    ev = float(r["best_evalue"])
    if pf not in best or ev < best[pf][0]:
        best[pf] = (ev, r["gene_locus_tag"], r["coordinates"])

with open("results/markdown/epigenetic_families.md", "w") as o:
    o.write("| Pfam | Family | Category | n genes | Example gene | Coordinates | E-value |\n")
    o.write("|------|--------|----------|--------:|--------------|-------------|---------|\n")
    for r in rows:
        pf = r["family_pfam"]
        b = best.get(pf)
        if b:
            o.write(f"| {pf} | {r['family_name']} | {r['category']} | {r['n_genes']} | "
                    f"{b[1]} | {b[2]} | {b[0]:.0e} |\n")
        else:
            o.write(f"| {pf} | {r['family_name']} | {r['category']} | {r['n_genes']} | "
                    f"— | — | — |\n")

with open("results/markdown/epigenetic_genes_full.md", "w") as o:
    o.write("| Checked family | Gene (locus_tag) | Gene coordinates |\n")
    o.write("|----------------|------------------|------------------|\n")
    for r in csv.DictReader(open("results/epigenetic_genes.tsv"), delimiter="\t"):
        o.write(f"| {r['family_name']} ({r['family_pfam']}) | {r['gene_locus_tag']} | {r['coordinates']} |\n")

STRUCTS = [("G4", "Quadruplexes"), ("Zhunt", "Z-Hunt")]
if os.path.exists("results/distribution/table1_ZDNABERT.tsv"):
    STRUCTS.append(("ZDNABERT", "Z-DNABERT"))

def load_t1(label):
    return {r["Region"]: r for r in csv.DictReader(
        open(f"results/distribution/table1_{label}.tsv"), delimiter="\t")}

t1 = {lab: load_t1(lab) for lab, _ in STRUCTS}
regions = ["Exon", "Intron", "Promoter", "Downstream", "Intergenic", "TOTAL"]
with open("results/markdown/distribution_table1.md", "w") as o:
    head = "| Region |"
    sep = "|--------|"
    for _, disp in STRUCTS:
        head += f" {disp} (N) | {disp} frac | {disp} enrich |"
        sep += "-----:|-----:|-----:|"
    head += " Background (bp frac) |"; sep += "-----:|"
    o.write(head + "\n" + sep + "\n")
    for reg in regions:
        row = f"| {reg} |"
        for lab, _ in STRUCTS:
            r = t1[lab][reg]
            row += f" {r['N_structures']} | {r['Fraction']} | {r['Enrichment']} |"
        bg = t1[STRUCTS[0][0]][reg]["Background_bp_fraction"] if reg != "TOTAL" else "1.0000"
        row += f" {bg} |"
        o.write(row + "\n")

def load_t2(label):
    return {r["Region"]: r for r in csv.DictReader(
        open(f"results/distribution/table2_{label}.tsv"), delimiter="\t")}

t2 = {lab: load_t2(lab) for lab, _ in STRUCTS}
with open("results/markdown/distribution_table2.md", "w") as o:
    head = "| Region | N regions |"
    sep = "|--------|----------:|"
    for _, disp in STRUCTS:
        head += f" Regions w/ {disp} | Frac ({disp}) |"
        sep += "-----:|-----:|"
    o.write(head + "\n" + sep + "\n")
    for reg in ["Exon", "Intron", "Promoter", "Downstream", "Intergenic"]:
        nreg = t2[STRUCTS[0][0]][reg]["N_regions"]
        row = f"| {reg} | {nreg} |"
        for lab, _ in STRUCTS:
            r = t2[lab][reg]
            row += f" {r['N_regions_with_structure']} | {r['Fraction_regions_with_structure']} |"
        o.write(row + "\n")

print("Markdown tables written to results/markdown/")
for f in sorted(os.listdir("results/markdown")):
    print("  ", f)
