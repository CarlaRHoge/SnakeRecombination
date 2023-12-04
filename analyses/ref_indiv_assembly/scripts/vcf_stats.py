#!/usr/bin/env python

import sys
from optparse import OptionParser
import itertools
import numpy as np

def mean(l):
    if len(l)==0:
        return "nan"
    return "{:.4f}".format(sum(l)/float(len(l)))

def safedivision(n,d):
    return n/float(d) if d!=0 else np.nan

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--info", action="store", type="str", dest="info",help="Comma separated list of INFO fields")
(options, args) = parser.parse_args()

info_params = options.info.split(",")
info_cumul = {p:[] for p in info_params}
ads = []
dps = []
afs = []
miss = []
polys = 0

for line in sys.stdin:

    line = line.strip()
    fields = line.split()
    polys += 1

    # INFO field params
    info = {x.split("=")[0]:float(x.split("=")[1]) for x in fields[7].split(";")}
    for p in info_params:
        if p in info:
            info_cumul[p].append(info[p])

    # Genotype field params
    genotypes = [g.split(":") for g in fields[9:]]
    guide = fields[8].split(":")

    # Jump if genotype is malformated for some sample
    if any(len(g)<3 for g in genotypes):
        continue

    # Get variant allele frequency
    gts = [sum(map(int, gt[0].replace("|","/").split("/"))) if "." not in gt[0] else "." for gt in genotypes]
    clean_gts = [g for g in gts if g!="."]
    af = safedivision(sum(clean_gts),(len(clean_gts)*2))
    afs.append(af if af<=0.5 else 1-af)
    
    # Missingness
    miss.append(safedivision(gts.count("."),len(gts)))
            
    for param in ["AD","DP"]:

        i = guide.index(param)

        # AD
        if param=="AD":
            het_ad = [list(map(int,gt[i].split(","))) for gt in genotypes if gt[0]=="0/1" or gt=="0|1" or gt=="1|0"]
            balance = ([safedivision(h[0],sum(h)) for h in het_ad])
            ads += balance

        if param=="DP":

            depths = [int(gt[i]) for gt in genotypes if gt[i]!="."]
            dps += depths

# Nucleotide diversity
# haps = {}
# diffs = []
# if len(haplotypes)>0:
#     for i,n in enumerate(haplotypes[0]):
#         haps[i] = [h[i] for h in haplotypes]
#     for pair in itertools.combinations(haps.values(),2):
#         d = sum([abs(n1-n2) for n1,n2 in zip(pair[0],pair[1]) if "."!=n1 and "."!=n2])
#         diffs.append(d)

# Missingness
# concat_haps = sum(list(haps.values()), [])
# if len(concat_haps)>0:
#     missingness = "{:.4f}".format(concat_haps.count(".")/len(concat_haps))
# else:
#     missingness = "nan"

# Print information
stats = ""
for param in info_params:
    stats += mean(info_cumul[param]) + "\t"
stats += mean(afs) + "\t"
stats += mean(ads) + "\t"
stats += mean(dps) + "\t"
stats += mean(miss) + "\t"
stats += str(polys)

print(stats)
                


    

    

    


