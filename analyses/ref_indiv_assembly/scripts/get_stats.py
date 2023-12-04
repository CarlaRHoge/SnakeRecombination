#!/usr/bin/env python

import sys


def safediv(n,d):
    return n/d if d!= 0 else "NaN"

# Process fasta
for line in sys.stdin:
    line = line.strip()
    if line.startswith(">"):
        seq = ""
    else:
        seq += line

# Report stats
rmsk = sum([n.islower() for n in seq])
seq = seq.upper()
l = float(len(seq))
#m = seq.count("MM")
c = seq.count("C") #+ m
g = seq.count("G") #+ m
gc = c+g
cpg = sum([seq.count(n) for n in ["CG","MM"]])
ns = seq.count("N")

print(safediv(gc,l),
      l,
      safediv(rmsk,l),
      safediv(ns, l), 
      safediv(cpg,l),
      safediv(safediv(cpg,l),safediv(c,l)*safediv(g,l)), 
      sep="\t")



