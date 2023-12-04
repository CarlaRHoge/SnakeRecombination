#!/usr/bin/env python

import sys

header = []
motif = []
names = []

for i,line in enumerate(sys.stdin):
    line = line.strip()
    fields = line.split()
    if i<=28:
        header.append(line)
    else:
        if len(fields)>0:
            motif.append(line)
        else:
            motif_name = motif[0].split()[1]
            names.append(motif_name)
            with open(f"streme/{motif_name}.txt","w") as fh:
                for h in header:
                    fh.write(h + "\n")
                for m in motif:
                    fh.write(m + "\n")
            motif = []
            
with open("streme/motif_names.txt", "w") as fh:
    for n in names:
        fh.write(n + "\n")


