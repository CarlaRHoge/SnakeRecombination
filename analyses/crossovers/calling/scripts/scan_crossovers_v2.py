#!/usr/bin/env python

import pandas as pd
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
import sys

def read_hapi_genotypes(file):
    df = pd.read_csv(file)
    df = df[df["Unnamed: 0"].str.contains("Super")]
    df["start"] = df.apply(lambda r: int(r["Unnamed: 0"].split("_")[-1]), axis=1)
    #df["end"] = [df.iloc[i+1]["start"] if i<df.index[-1] else np.nan for i,r in df.iterrows()]
    ndf = df[[c for c in df.columns if "-" in c or c in ["start"]]]
    ndf = ndf.replace("0", np.nan).replace("A", 1).replace("a", np.nan).replace("B", 2).replace("b", np.nan)
    return df, ndf

def scan_phase_changes(hapi, parents):
    
    df_complete = hapi[parents + ["start"]]
    df_informative = df_complete[df_complete.apply(lambda x: "A" in x.values or "B" in x.values, axis=1)]
    df_informative = df_informative.replace("A",2).replace("B",1).replace("0",np.nan)
    df_nomiss = df_informative.dropna()
    
    indexes = df_nomiss.index
    ph = {chrom:[] for chrom in parents}

    for i in range(0, len(indexes)-1):

        i1 = indexes[i]
        i2 = indexes[i+1]
        window = df_nomiss.loc[i1:i2]

        phase_changes = [set(window[chrom].values)=={1,2} for chrom in parents]
        if len(set(phase_changes))==1:
            continue
        
        # Distance between informative sites is >200kb, skip
        #if (window.start.max() - window.start.min()) > 200e3:
        #    continue

        minority_bool = True if phase_changes.count(True)<=2 else False
        minority_parents = [p for p,b in zip(parents, phase_changes) if b==minority_bool]
        for chrom in minority_parents:
            ph[chrom].append([i,i1,i2])
    
    #clean_crossovers(crossovers, df_complete)
    return ph, df_informative

def get_crossovers(phase_changes, df_polish, indiv):

    crossovers = []
    
    # Group changes of phase in a row, if changes are not odd
    block = []
    groups = []
    min_info_markers = 10
    
    for i,change in enumerate(phase_changes):
        df_loc, start, end = change
        if i==0:
            prev = df_loc
            if prev >= min_info_markers:
                block.append(change)
        else:
            if df_loc-prev <= min_info_markers:
                block.append(change)
            else:
                groups.append(block)
                block = [change]
            prev = df_loc

    groups.append(block)
    if len(phase_changes)==2:
        groups = []

    for block in groups:
        # Crossovers are only odd changes of phase
        if len(block)%2!=0:
            
            # Try to narrow down the crossover, if only 1 clear phase change
            if len(block)==1:
                itr, start, end = block[0]
                refine = df_polish.loc[start:end][[indiv,"start"]].dropna()
                refine = refine#.reset_index() 
                for absi,(i,r) in enumerate(refine.iterrows()):
                    if absi==0:
                        phase = r[indiv]
                    else:
                        if r[indiv]!=phase:
                            cross_start = refine.iloc[absi-1]["start"]
                            cross_end   = refine.iloc[absi]["start"]
                            break
                        phase = r[indiv]
                # If the template is staying in phase, but the other changing
                if len(set(refine[indiv].dropna().values))==1:
                    cross_start = refine.start.min()
                    cross_end = refine.start.max()
                    
                crossovers.append([cross_start, cross_end, start, end])
            
            # If multiple (but odd) changes of chase, get narrower range possible
            else:
                start = block[0][1]
                end   = block[-1][-1]
                cross_start = df_polish.loc[start]["start"]
                cross_end   = df_polish.loc[end]["start"]
                crossovers.append([cross_start, cross_end, start, end])
    
    index_list = list(df_informative.index)
    maxl = len(index_list)
    info_markers_limit = 30
    
    clean_crossovers = []
    for c in crossovers:
        istart = index_list.index(c[2])
        iend = index_list.index(c[3])
        if istart<=info_markers_limit or iend>=(maxl-info_markers_limit):
            continue
        clean_crossovers.append(c)
            
    
    return clean_crossovers


usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-s", "--scaf", action="store", type="str", dest="scaf",help="Scaf string")
parser.add_option("-f", "--family", action="store", type="str", dest="family",help="Family string")
parser.add_option("-c", "--complex", action="store", type="str", dest="famcomplex",help="Family complex string")
(options, args) = parser.parse_args()

scaf = options.scaf
fam = options.family
famcomplex = options.famcomplex
hapi_dir = "/moto/palab/projects/recombination_snakes/shared_corn_rattle/crossovers/calling"
df,ndf = read_hapi_genotypes("{}/hapi_output/{}.{}/{}.01.csv".format(hapi_dir, scaf, fam, famcomplex))

maternal = [c for c in ndf.columns if "M" in c]
paternal = [c for c in ndf.columns if "P" in c]
crossovers = {}

# Detect changes of phase and potential crossovers
for parents in [maternal, paternal]:
    ph_change, df_informative = scan_phase_changes(df, parents)
    for p in parents:
        c = get_crossovers(ph_change[p], df_informative, p)
        crossovers[p] = c

bed_l = []
for pchrom in crossovers:
    for cross in crossovers[pchrom]:
        bed_l.append([scaf, int(cross[0]), int(cross[1]), pchrom, pchrom.split("-")[0]])
bed_df = pd.DataFrame(bed_l)

# If no crossovers in scaffold
if len(bed_df)==0:
    f,ax = plt.subplots()
    subd = ndf
    sns.heatmap(subd[[c for c in ndf if "-" in c]].transpose(), ax=ax, cmap="viridis")

    xticklabel = [[i,l] for i,l in enumerate(subd.start) if i%int(len(subd)/15)==0]
    ax.set_xticks([x[0] for x in xticklabel])
    ax.set_xticklabels(["{:,}".format(x[1]) for x in xticklabel])
    plt.tight_layout()
    plt.savefig("hapi_pdfs/{}.{}.jpg".format(scaf, fam), dpi=200)
    sys.exit()
    
bed_df.columns = ["chrom", "start", "end", "chrom", "sex"]

# Clean up crossovers next to each other in one parent
discard_db = pd.DataFrame()

for sex,df in bed_df.groupby("sex"):
    bed = BedTool.from_dataframe(df)
    bed = bed.sort()
    bed_merged = bed.merge(d=5000, c=4, o="collapse")
    df_merged = bed_merged.to_dataframe()
    if len(df_merged)==0:
        continue
    df_to_discard = df_merged[df_merged.apply(lambda x: len(x["name"].split(",")) > 1, axis=1)]
    to_discard = BedTool.from_dataframe(df_to_discard)
    df_discarded = bed.intersect(to_discard).to_dataframe()
    df_kept = bed.intersect(to_discard, v=True).to_dataframe()
    discard_db = pd.concat([discard_db, df_discarded])

discard_ids = {}
if len(discard_db)>0:
    discard_ids = {c:[[r.start, r.end] for i,r in df.iterrows()] for c,df in discard_db.groupby("name")}

# Plotting and output
f,ax = plt.subplots()

start = 1
end = 10e6

subd = ndf.loc[start:end]
sns.heatmap(subd[[c for c in ndf if "-" in c]].transpose(), ax=ax, cmap="viridis")

xticklabel = [[i,l] for i,l in enumerate(subd.start) if i%int(len(subd)/15)==0]
ax.set_xticks([x[0] for x in xticklabel])
ax.set_xticklabels(["{:,}".format(x[1]) for x in xticklabel])

for chrom in crossovers:
    y = list(subd.columns).index(chrom) + 0.5
    for cross in crossovers[chrom]:
        start_pos,end_pos = cross[0],cross[1]
        loc = cross[-2] + ((cross[-1]-cross[-2])/2)
        if [start_pos, end_pos] in sum(list(discard_ids.values()),[]):
            ax.scatter(loc-start, y, 
                       color="red", 
                       edgecolor="black", 
                       linewidth=0.2)
        else:
            print("{}\t{}\t{}\t{}".format(scaf, int(cross[0]), int(cross[1]), chrom))
            ax.scatter(loc-start, y, 
                       color="lightgreen" if end_pos-start_pos<20e3 else "orange",
                       edgecolor="black", 
                       linewidth=0.2)

plt.tight_layout()
plt.savefig("hapi_pdfs/{}.{}.jpg".format(scaf, fam), dpi=200)

