#!/usr/bin/env python3
import re, sys

GENOME = "data/genome.fna"
OUT = "results/quadruplexes.bed"

plus_re = re.compile(r"(?:G{3,5}[ACGT]{1,7}){3,}G{3,5}", re.IGNORECASE)
minus_re = re.compile(r"(?:C{3,5}[ACGT]{1,7}){3,}C{3,5}", re.IGNORECASE)

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)

n_plus = n_minus = 0
with open(OUT, "w") as out:
    for name, seq in read_fasta(GENOME):
        for m in plus_re.finditer(seq):
            out.write(f"{name}\t{m.start()}\t{m.end()}\tG4_plus\t0\t+\n")
            n_plus += 1
        for m in minus_re.finditer(seq):
            out.write(f"{name}\t{m.start()}\t{m.end()}\tG4_minus\t0\t-\n")
            n_minus += 1
        sys.stderr.write(f"  {name}: len={len(seq):,}\n")

print(f"G4 on + strand: {n_plus}")
print(f"G4 on - strand: {n_minus}")
print(f"Total PQS: {n_plus + n_minus} -> {OUT}")
