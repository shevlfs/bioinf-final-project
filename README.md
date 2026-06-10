# *Caenorhabditis auriculariae* — genetic, epigenetic and secondary-structure annotation

**Course:** HSE minor, "Early evolution of invertebrates" project (2026)
**Taxon (group):** genus *Caenorhabditis* (Nematoda; рус. *ценорабдитис*)
**Individual organism:** *Caenorhabditis auriculariae* (strain NKZ352)
**Assembly:** `GCA_904845305.2` (CAUJ3.4), GenBank, chromosome-scale scaffolds

This repository contains my individual part of the project: download of the genome
and annotation, identification of epigenetic-machinery gene families, and genome-wide
annotation of two non-B DNA secondary structures — **G-quadruplexes** and **Z-DNA** —
followed by their distribution across genomic features.

---

## 1. The organism

*Caenorhabditis auriculariae* is a free-living, bacterivorous nematode. It was
morphologically described in 1999 and re-isolated by Kanzaki et al. (2021) from a
*Platydema* sp. darkling beetle (Coleoptera: Tenebrionidae) in Japan, where it lives
in a **necromenic / phoretic association** with the beetle and feeds on bacteria in
decaying wood and fungal (mushroom) substrates. The species epithet *auriculariae*
refers to the jelly fungus *Auricularia*, the substrate of the original description.
The laboratory strain **NKZ352** is maintained on standard nematode growth medium
(NGM) at ~20–25 °C in Taisei Kikuchi's laboratory (University of Miyazaki).

**Why it matters for this project.** Molecular phylogenetics (rRNA + 269 single-copy
genes) places *C. auriculariae* together with *C. sonorae* and *C. monodelphis* at the
**most basal clade of the genus *Caenorhabditis***. It is therefore the natural
"early-branching" outgroup of the genus, and a good organism for asking *when* an
epigenetic mark or secondary structure first appeared during nematode evolution.

| Property | Value |
|---|---|
| NCBI accession | GCA_904845305.2 (assembly CAUJ3.4) |
| Assembly level | Chromosome-scale scaffolds (Hi-C / 3D-DNA) |
| Total length | 111,858,439 bp (~112 Mb) |
| Scaffolds | 6 → **Chr1–Chr5 + ChrX** (genus karyotype n = 6) |
| Contigs | 35 |
| Scaffold N50 | 18,788,797 bp |
| Contig N50 | 14,535,289 bp |
| Gaps (N) | low; 35 contigs in 6 scaffolds, i.e. ~29 internal gaps |
| GC content | 38 % |
| Protein-coding genes | 15,814 (18,753 protein isoforms) |
| Annotated exons | 389,135 |
| Sequencing | Oxford Nanopore + Illumina HiSeq |
| Submitter | University of Miyazaki (Kikuchi lab) |
| Genome paper | Kanzaki et al. 2021, *Sci. Rep.* 11:6720 (PMID 33762598) |

PubMed hits for `"Caenorhabditis auriculariae"`: a handful (the 2021 genome paper and
phylogenetic notes); the species itself is little-studied, but the genus *Caenorhabditis*
(via the model *C. elegans*) is one of the best-characterised in all of biology, which
provides the epigenetic background below.

### Epigenetic background of the genus (from *C. elegans*)
- **DNA methylation:** *Caenorhabditis* famously **lacks canonical 5-methyl-cytosine
  DNA methylation** — no DNMT1/DNMT3 cytosine methyltransferases. Only the conserved
  **DNMT2/TRDMT1** (a tRNA methyltransferase) is retained. N6-methyl-adenine (6mA) has
  been reported but is contentious.
- **Histone modifications are extensive:** H3K4me (SET-2/COMPASS), H3K9me (MET-2/SET-25),
  H3K27me (PRC2 = MES-2/MES-3/MES-6), H3K36me (MES-4/MET-1), plus HATs, HDACs and
  sirtuins. These drive germline immortality, dosage compensation and heterochromatin.
- **3D genome:** *Caenorhabditis* has cohesin/condensin (SMC complexes) but **no CTCF**
  ortholog — TADs are organised differently from vertebrates.

---

## 2. Methods & parameters

All steps are pure-Python / standard tools; scripts are in [`scripts/`](scripts/).

| Step | Tool / method | Key parameters |
|---|---|---|
| Genome, annotation, proteome | NCBI Datasets REST API v2 | accession GCA_904845305.2; FASTA + GFF3 + protein FASTA |
| Epigenetic gene families | **HMMER 3.4** `hmmsearch` vs 26 curated Pfam HMMs | Pfam gathering thresholds (`--cut_ga`) |
| G-quadruplexes | regex on **both strands**, case-insensitive | pattern `(G{3,5}[ATGC]{1,7}){3,}G{3,5}` (and C-pattern for the − strand) |
| Z-DNA (thermodynamic) | **Z-Hunt-II** (canonical C tool), cross-validated vs a NumPy re-implementation | window/maxdin = 12, mindin = 6 dinucleotides; **z-score > 400** |
| Z-DNA (deep learning) | **Z-DNABERT** (DNABERT-6, HG-kouzine model) on Apple-Silicon **GPU (MPS)** | 6-mers, 512-token windows; **p > 0.5**, region length > 10 |
| Feature distribution | per-base priority partition in NumPy | Exon > Intron > Promoter (1000 bp up of TSS) > Downstream (200 bp) > Intergenic |

**Notes on the tools**

- **HMMER / Pfam.** The protein FASTA (18,753 sequences) was searched against 26 HMM
  profiles fetched from InterPro covering the three required classes — **DNA
  (de)methylation, histone-fold / histone-like proteins, and histone-modifying enzymes**
  — plus chromatin remodelers and the cohesin complex (for the group part). Hits were
  collapsed to genes via the GFF `locus_tag`.
- **Z-Hunt.** The Z-DNA calls (`results/zhunt.bed`) were produced with the **canonical
  Z-Hunt-II C tool** (Ho et al., *EMBO J.* 1986; `tools/zhunt-reference.c`, run via
  `scripts/run_zhunt_c.py`, parallelised over genome chunks on 14 cores with `windowsize
  12 / min 6 / max 12`). To make the result trustworthy it was **cross-validated against
  an independent NumPy re-implementation** of the same algorithm (`scripts/zhunt.py`),
  which agrees with the C tool to **mean |Δz| = 0.07** over 200 kb and **83,225 vs 83,213
  Z-DNA regions genome-wide (0.01 % difference)** — the only differences are a handful of
  borderline (z ≈ 400) windows where the two tools resolve exact anti/syn-energy ties
  differently. `scripts/zhunt.py --selftest` also reproduces a scalar transliteration of
  the C code to **0 difference**. Z-Hunt is a property of the **duplex** (strand-symmetric),
  so one pass annotates both strands; the genomic cut-off **z-score > 400** was applied and
  adjacent passing positions merged.
- **Z-DNABERT.** As a second, independent Z-DNA caller, the deep-learning model Z-DNABERT
  (Beknazarov et al. 2023; DNABERT-6 fine-tuned on human Kouzine Z-flipons,
  `mitiau/Z-DNABERT`) was run genome-wide on the **Apple-Silicon GPU via PyTorch MPS**
  (`scripts/run_zdnabert.py`): sequence → overlapping 6-mers → 512-token windows →
  per-token Z-DNA probability → threshold p > 0.5 → regions > 10 bp. It is far more
  conservative than Z-Hunt (13,375 vs 83,225 sites) but **92.6 % of its regions overlap a
  Z-Hunt site** — two orthogonal methods (physics-based and ML-based) independently flag
  the same loci, which strongly validates the Z-DNA annotation.

---

## 3. Results — epigenetic gene families

**24 of 26 families** searched have at least one gene; **196 distinct genes** were
identified. Full per-gene table: [`results/epigenetic_genes.tsv`](results/epigenetic_genes.tsv)
(required `family → gene → coordinates` format also in
[`results/markdown/epigenetic_genes_full.md`](results/markdown/epigenetic_genes_full.md)).
Per-family summary: [`results/epigenetic_summary.tsv`](results/epigenetic_summary.tsv).

| Pfam | Family | Category | n genes | Example gene | Coordinates | E-value |
|------|--------|----------|--------:|--------------|-------------|---------|
| PF00145 | DNA_methylase | DNA methylation | 1 | CAUJ4_LOCUS8523 | Chr4:7384415-7386554(−) | 9e-44 |
| PF12047 | DNMT1-RFD | DNA methylation | **0** | — | — | — |
| PF02008 | zf-CXXC | DNA methylation | 1 | CAUJ4_LOCUS6682 | Chr3:11566140-11568680(−) | 1e-15 |
| PF01429 | MBD | DNA methylation read | 2 | CAUJ4_LOCUS7549 | Chr3:18689295-18696356(−) | 2e-19 |
| PF12851 | Tet_JBP | DNA modification | **0** | — | — | — |
| PF00125 | Histone | Histone fold | 29 | CAUJ4_LOCUS5957 | Chr3:7811358-7811871(+) | 5e-39 |
| PF00538 | Linker_histone | Histone fold | 2 | CAUJ4_LOCUS7428 | Chr3:17871622-17872536(+) | 5e-31 |
| PF00856 | SET | Histone write (methyl) | 25 | CAUJ4_LOCUS3706 | Chr2:8972070-8975227(+) | 9e-27 |
| PF05033 | Pre-SET | Histone write (methyl) | 2 | CAUJ4_LOCUS7449 | Chr3:18023212-18035018(−) | 1e-16 |
| PF01853 | MOZ_SAS | Histone write (acetyl) | 6 | CAUJ4_LOCUS12038 | Chr5:10589978-10594008(−) | 2e-82 |
| PF00583 | Acetyltransf_1 | Histone write (acetyl) | 30 | CAUJ4_LOCUS12077 | Chr5:10752817-10755962(−) | 2e-21 |
| PF08214 | HAT_KAT11 | Histone write (acetyl) | 1 | CAUJ4_LOCUS6454 | Chr3:10468719-10482099(+) | 3e-64 |
| PF00850 | Hist_deacetyl | Histone erase (deacetyl) | 9 | CAUJ4_LOCUS776 | Chr1:7539332-7542321(+) | 3e-90 |
| PF02146 | SIR2 | Histone erase (deacetyl) | 4 | CAUJ4_LOCUS9578 | Chr4:12851386-12856837(+) | 1e-65 |
| PF02373 | JmjC | Histone erase (demethyl) | 7 | CAUJ4_LOCUS7669 | Chr4:560945-580136(−) | 3e-41 |
| PF02375 | JmjN | Histone erase (demethyl) | 2 | CAUJ4_LOCUS7669 | Chr4:560945-580136(−) | 6e-15 |
| PF01593 | Amino_oxidase | Histone erase (demethyl) | 10 | CAUJ4_LOCUS614 | Chr1:6547484-6555393(+) | 7e-70 |
| PF00628 | PHD | Histone read | 21 | CAUJ4_LOCUS4874 | Chr2:17612649-17629295(+) | 1e-19 |
| PF00439 | Bromodomain | Histone read | 15 | CAUJ4_LOCUS930 | Chr1:8609816-8631023(−) | 9e-102 |
| PF00385 | Chromo | Histone read | 13 | CAUJ4_LOCUS13687 | Chr6/X:4826609-4832140(+) | 2e-38 |
| PF00855 | PWWP | Histone read | 1 | CAUJ4_LOCUS1552 | Chr1:13130063-13135482(+) | 3e-16 |
| PF00567 | Tudor | Histone read | 11 | CAUJ4_LOCUS3266 | Chr2:6738331-6742652(−) | 3e-25 |
| PF01426 | BAH | Histone read | 3 | CAUJ4_LOCUS930 | Chr1:8609816-8631023(−) | 8e-34 |
| PF00176 | SNF2_N | Chromatin remodeling | 16 | CAUJ4_LOCUS7272 | Chr3:16077121-16096379(+) | 1e-107 |
| PF02463 | SMC_N | Cohesin / SMC | 8 | CAUJ4_LOCUS5632 | Chr3:5445890-5463321(+) | 1e-70 |
| PF04825 | Rad21_Rec8 | Cohesin | 3 | CAUJ4_LOCUS13968 | Chr6/X:6763387-6779857(−) | 1e-40 |

(Coordinates use the GenBank scaffold ids `CAJGYM0200000XX.1`; here abbreviated Chr1–Chr5/ChrX.)

### 3.1 Key finding — DNA methylation machinery is reduced
- Exactly **one** C-5 cytosine-methyltransferase gene (PF00145), a **330-aa** protein —
  the size of **DNMT2/TRDMT1** (a tRNA methyltransferase), *not* the ~1600-aa DNMT1.
- **No DNMT1** (the DNMT1-RFD targeting domain PF12047 is absent) and **no DNMT3**.
- **No TET/JBP** 5mC-oxidase (PF12851 absent).
- 1 CXXC zinc-finger and 2 MBD (methyl-CpG-reader) genes are present, but the **writers**
  of 5mC are missing.

➡ *C. auriculariae*, the basal-most *Caenorhabditis*, **already lacks canonical 5mC DNA
methylation**, exactly like the derived *C. elegans*. The loss of DNMT1/DNMT3 therefore
predates the radiation of the genus.

### 3.2 Group part — CTCF & cohesin
- **Cohesin is present:** SMC subunits (PF02463, 8 genes incl. SMC1/SMC3-type) and the
  kleisin **RAD21/REC8** (PF04825, 3 genes). The complex that extrudes chromatin loops
  exists.
- **CTCF is absent:** no CTCF-type multi-zinc-finger insulator was recovered, consistent
  with the well-established absence of CTCF in nematodes. Loop anchoring in this lineage
  does not use a CTCF/cohesin system.

### 3.3 Histone machinery is complete
All histone-modifying classes are richly represented — **writers** (SET methyl-, MYST &
GNAT acetyl-transferases, 25 + 36 genes), **erasers** (HDAC, sirtuin, Jumonji & LSD1
demethylases), **readers** (PHD, bromo, chromo, Tudor, PWWP, BAH) and **histone-fold**
genes (29) plus linker histone H1 (2). The chromatin-remodeling ATPase family SNF2 has
16 members.

---

## 4. Results — non-B DNA secondary structures

| Structure | Method | Count | Strands |
|---|---|---:|---|
| G-quadruplex (PQS) | regex `(G{3,5}[ATGC]{1,7}){3,}G{3,5}` | **4,918** | 2,652 (+) / 2,266 (−) |
| Z-DNA (Z-Hunt) | Z-Hunt-II, canonical C (z-score > 400) | **83,225** | duplex (strand-symmetric) |
| Z-DNA (Z-DNABERT) | DNABERT-6, MPS GPU (p > 0.5) | **13,375** | duplex |

BED files: [`results/quadruplexes.bed`](results/quadruplexes.bed),
[`results/zhunt.bed`](results/zhunt.bed), [`results/zdnabert.bed`](results/zdnabert.bed).
Both Z-DNA callers spread evenly over all six chromosomes; **92.6 % of Z-DNABERT regions
overlap a Z-Hunt site**. G4 median length 22 bp; Z-Hunt z-score median 619 (max ~1×10⁷);
Z-DNABERT probability median 0.91.

### Table 1 — distribution of structures across features (structure-centric, vs background)
"Fraction" = share of all structures of that type; "Background" = share of genomic **bp**
in that feature; "enrich" = fraction / background ( >1 enriched, <1 depleted).

Columns per method: N structures, fraction of that method's calls, and enrichment
(fraction ÷ background bp-fraction; **>1 enriched, <1 depleted**).

| Region | Quadr. (N) | Quadr. frac | Quadr. enr | Z-Hunt (N) | Z-Hunt frac | Z-Hunt enr | Z-DNABERT (N) | ZDB frac | ZDB enr | Background |
|--------|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Exon | 519 | 0.106 | **0.56** | 18,871 | 0.227 | **1.21** | 3,980 | 0.298 | **1.58** | 0.188 |
| Intron | 2,946 | 0.599 | **1.29** | 35,329 | 0.424 | 0.92 | 5,313 | 0.397 | 0.86 | 0.463 |
| Promoter (1 kb up TSS) | 324 | 0.066 | **0.64** | 9,957 | 0.120 | **1.15** | 1,399 | 0.105 | 1.01 | 0.104 |
| Downstream (200 bp) | 40 | 0.008 | 0.43 | 942 | 0.011 | 0.60 | 105 | 0.008 | 0.42 | 0.019 |
| Intergenic | 1,089 | 0.221 | 0.98 | 18,126 | 0.218 | 0.96 | 2,578 | 0.193 | 0.85 | 0.227 |
| **Total** | **4,918** | 1.000 | — | **83,225** | 1.000 | — | **13,375** | 1.000 | — | 1.000 |

### Table 2 — fraction of feature regions that contain a structure (region-centric)

| Region | N regions | w/ quadr. | frac | w/ Z-Hunt | frac | w/ Z-DNABERT | frac |
|--------|----------:|----------:|-----:|----------:|-----:|-------------:|-----:|
| Exon | 169,476 | 172 | 0.0010 | 7,497 | 0.0442 | 3,874 | 0.0229 |
| Intron | 154,089 | 2,125 | 0.0138 | 9,439 | 0.0613 | 4,400 | 0.0286 |
| Promoter | 13,497 | 280 | 0.0207 | 3,155 | **0.2338** | 1,248 | 0.0925 |
| Downstream | 10,797 | 40 | 0.0037 | 372 | 0.0345 | 109 | 0.0101 |
| Intergenic | 8,305 | 650 | 0.0783 | 3,078 | **0.3706** | 1,612 | 0.1941 |

### Interpretation (comparison with background)
- **G-quadruplexes are enriched in introns (1.29×)** and **depleted in exons (0.56×) and
  promoters (0.64×)**. Unlike mammals — where G4s mark promoters/CpG islands — in this
  GC-poor (38 %) nematode genome PQS avoid coding/regulatory sequence and accumulate in
  introns and intergenic DNA.
- **Z-DNA is enriched in exons and promoters** by *both* callers — Z-Hunt (exon 1.21×,
  promoter 1.15×) and Z-DNABERT (exon **1.58×**, promoter 1.01×) — and depleted downstream
  (0.60/0.42×) and in introns. Almost a quarter of all promoters (23 %) carry a Z-Hunt
  site. This reflects Z-DNA's preference for GC / purine-pyrimidine alternations in
  transcriptionally active 5′ and coding regions — a pattern conserved with vertebrates.
  The conservative, experimentally-trained Z-DNABERT sharpens the exon signal, and its
  92.6 % overlap with Z-Hunt shows the two methods converge on the same biology.

---

## 5. Repository layout

```
data/        genome.fna, proteins.faa, annotation.gff, cds.fna  (from GCA_904845305.2)
hmm/         families.tsv + 26 Pfam HMM profiles + pressed DB
results/
  quadruplexes.bed                 G4 predictions (BED6, both strands)
  zhunt.bed                        Z-DNA predictions, Z-Hunt z-score > 400 (BED6)
  zdnabert.bed                     Z-DNA predictions, Z-DNABERT p > 0.5 (BED6)
  epigenetic_genes.tsv             family → gene → coordinates (per gene)
  epigenetic_summary.tsv           per-family gene counts (26 families)
  hmmer/                           raw hmmsearch domtbl/tbl output
  distribution/                    Table 1 & Table 2 per structure (TSV)
  markdown/                        the tables above, as markdown
scripts/
  epigenetics_table.py             map HMMER hits → genes
  find_quadruplexes.py             G4 regex, both strands
  zhunt.py                         NumPy Z-Hunt-II port (+ --selftest), validation
  run_zhunt_c.py                   run canonical Z-Hunt C tool, multicore (production)
  run_zdnabert.py                  Z-DNABERT (DNABERT-6) inference on GPU/MPS
  genomic_distribution.py          feature partition + Tables 1 & 2
  make_report_tables.py            markdown table generator
tools/       zhunt-reference.c     reference Z-Hunt C source (production Z-DNA tool)
models/      zdnabert/6-new-12w-0/ Z-DNABERT DNABERT-6 weights (gitignored; see §6)
```

## 6. Reproduce

```bash
# 0. genome.fna and annotation.gff are committed gzipped (GitHub 100 MB limit) — unpack
gunzip -k data/genome.fna.gz data/annotation.gff.gz

# 1. data (NCBI Datasets REST API) — or regenerate the package from the accession
curl -L "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_904845305.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF" -o data/caur_package.zip

# 2. epigenetic gene families (HMMER 3.4)
hmmsearch --cut_ga --cpu 6 --domtblout results/hmmer/epi_domtbl.txt \
          hmm/epigenetic_families.hmm data/proteins.faa
python3 scripts/epigenetics_table.py

# 3. G-quadruplexes
python3 scripts/find_quadruplexes.py

# 4. Z-DNA: canonical Z-Hunt C tool (cross-checked by the NumPy port)
gcc tools/zhunt-reference.c -O2 -o tools/zhunt_c -lm
python3 scripts/zhunt.py --selftest          # NumPy port == scalar reference (Δ=0)
python3 scripts/run_zhunt_c.py --mindin 6 --maxdin 12 --threshold 400 --procs 14 \
                         --genome data/genome.fna --out results/zhunt.bed

# 5. Z-DNABERT (deep-learning Z-DNA) on Apple-Silicon GPU (MPS); weights from Google Drive
#    (gdown the 5 files of the 'HG kouzine' model into models/zdnabert/6-new-12w-0/)
.venv/bin/python scripts/run_zdnabert.py --genome data/genome.fna \
        --out results/zdnabert.bed --batch 64 --fp16

# 6. distribution tables (all three structures)
python3 scripts/genomic_distribution.py G4=results/quadruplexes.bed \
        Zhunt=results/zhunt.bed ZDNABERT=results/zdnabert.bed
python3 scripts/make_report_tables.py
```

## 7. Conclusions (for the "early evolution" question)

- **DNA 5mC methylation:** *C. auriculariae* — the **earliest-branching** member of the
  genus — **already lacks the DNMT1/DNMT3 writers and TET erasers**, keeping only the
  ancestral DNMT2/TRDMT1. So the loss of cytosine-DNA-methylation machinery is **ancient
  within nematodes**, present before the *Caenorhabditis* genus diversified; this mark is
  *not* a derived feature of *C. elegans* alone. By contrast, vertebrates and many other
  invertebrates retain DNMT1/DNMT3 — the methylation toolkit is the ancestral eukaryotic
  state that was *lost* in this lineage.
- **Histone modifications** (writers, erasers, readers of methyl/acetyl marks) and the
  histone-fold proteins are **fully present** in the basal nematode — these epigenetic
  layers are ancient and conserved across the eukaryotes/invertebrates.
- **3D-genome proteins:** **cohesin is present but CTCF is absent**, confirming that
  CTCF-based insulation is a vertebrate/bilaterian innovation not used by nematodes.
- **Secondary structures:** both G-quadruplex- and Z-DNA-forming potential are abundant
  in this basal nematode (4,918 PQS; 83,225 Z-Hunt and 13,375 Z-DNABERT Z-DNA sites, with
  92.6 % agreement between the two independent Z-DNA methods), with Z-DNA enriched at
  promoters/exons — i.e. these non-B structures are deep, conserved features of genome
  organisation rather than recent additions.

## References
- Kanzaki N. *et al.* (2021) *Additional description and genome analyses of Caenorhabditis
  auriculariae representing the basal lineage of genus Caenorhabditis.* **Sci. Rep.**
  11:6720. PMID 33762598.
- Ho P.S. *et al.* (1986) *A computer-aided thermodynamic approach for predicting the
  formation of Z-DNA in naturally occurring sequences.* **EMBO J.** 5:2737–2744. (Z-Hunt)
- Beknazarov N. *et al.* (2023) *Biological roles for Z-DNA and Z-RNA revealed by deep
  learning.* **Life Sci. Alliance** 6:e202301962. (Z-DNABERT; `github.com/mitiau/Z-DNABERT`)
- Ji Y. *et al.* (2021) *DNABERT: pre-trained BERT for DNA-language in genome.*
  **Bioinformatics** 37:2112. (base model)
- Mistry J. *et al.* (2021) *Pfam: The protein families database in 2021.* **NAR** 49:D412.
- Eddy S.R. (2011) *Accelerated profile HMM searches.* **PLoS Comput. Biol.** 7:e1002195. (HMMER)
