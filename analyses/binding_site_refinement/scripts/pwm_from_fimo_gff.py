#/usr/bin/env python

import sys

seqs = []
compl = {"A":"T","C":"G","G":"C","T":"A"}

for line in sys.stdin:
    if not line.startswith("#"):
        fields = line.split()
        strand = fields[6]
        dna = fields[-1].split("=")[-1].replace(";","").upper()
        seqs.append(dna)

nucs = ["A", "C", "G", "T"]
freqs = {}
for i,n in enumerate(nucs):
    freqs[n] = []
    for i2,a in enumerate(seqs[0]):
        sel_a = [s[i2] for s in seqs]
        f = sel_a.count(n)/len(sel_a)
        freqs[n].append("{:.3f}".format(f))
        
for n in nucs:
    print("{}:   {}".format(n, "   ".join(map(str,freqs[n]))))


    