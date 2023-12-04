#/usr/bin/env python

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from optparse import OptionParser
import itertools

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-o", "--pdf", action="store", type="str", dest="pdf")
parser.add_option("-f", "--file_list", action="store", type="str", dest="file_list")
parser.add_option("-l", "--label_list", action="store", type="str", dest="label_list")
(options, args) = parser.parse_args()

files = options.file_list.split(",")
labels = options.label_list.split(",")
dfs = []

nucs = ["A","C","T","G"]

for fh,label in zip(files, labels):
    with open(fh, "r") as ofh:
        rows = []
        for line in ofh:
            if line.startswith(">"):
                continue
            fields = line.split()
            rows.append(list(map(float,fields[1:])))
    df = pd.DataFrame(rows).transpose()
    df.columns = nucs
    df["pwm"] = label
    dfs.append(df)
    
mdf = pd.concat(dfs)

labels = ['original', 'hot', 'hotFarCpGi']

f,axs0 = plt.subplots(1,3, figsize=[9,3])
axs = axs0.reshape(-1)

i = 0
for comb in list(itertools.combinations(labels,2)):
    l1,l2 = comb
    if l1!=l2:
        d1 = mdf[mdf.pwm==l1].reset_index(drop=True)
        d2 = mdf[mdf.pwm==l2].reset_index(drop=True)
        rdf_l = []
        for n in nucs:
            rdf = pd.DataFrame([d1[n],d2[n]]).transpose()
            rdf.columns = [l1, l2]
            rdf["nuc"] = n
            rdf_l.append(rdf)
            sns.regplot(data=rdf, 
                        x=l1, 
                        y=l2,
                        scatter_kws = {"s":8},
                        ax=axs[i], 
                        label=n)
        mrdf = pd.concat(rdf_l)
    axs[i].axline((0, 0), 
                  slope=1, 
                  color='lightgray', 
                  linestyle="dotted")
    axs[i].set_xlim([0,1])
    axs[i].set_ylim([0,1])
    if i==2:
        axs[i].legend()
    i += 1
        
plt.tight_layout()
sns.despine()
plt.savefig(options.pdf)