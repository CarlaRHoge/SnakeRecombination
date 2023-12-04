#!/usr/bin/env python

#from optparse import OptionParser
import sys
import numpy as np

# usage = "usage: %prog [options] arg1 arg2"
# parser = OptionParser(usage=usage)
# parser.add_option("-r", "--ref", action="store", type="str", dest="ref",help="Reference species")
# (options, args) = parser.parse_args()

haplotypes = {}
substrate = {}
diffs = {}
nucleotides = ["A","C","T","G"]

for line in sys.stdin:
    if line.startswith("#"):
        continue
    elif line.startswith("a"):
        if len(haplotypes)>0:
            for sp,dna in haplotypes.items():
                for n1,n2 in zip(haplotypes[ref], haplotypes[sp]):
                    if n1 in nucleotides and n2 in nucleotides:
                        if sp not in substrate:
                            substrate[sp] = 0
                            diffs[sp] = 0
                        substrate[sp] += 1
                        if n1!=n2:
                            diffs[sp] += 1
            
    elif line.startswith("s"):
        fields = line.split()
        sp = fields[1].split(".")[0]
        dna = fields[-1]
        if len(haplotypes)==0:
            ref = sp
        haplotypes[sp] = dna.upper()

for sp in diffs:
    print(ref,
          sp,
          diffs[sp], 
          substrate[sp],
          diffs[sp]/substrate[sp] if substrate[sp]>0 else np.nan, 
          sep="\t")
        