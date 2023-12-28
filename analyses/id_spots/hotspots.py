import sys
from glob import glob
import os
import subprocess
import pandas as pd
import _pickle as pickle
import numpy as np
import random
import statistics
import csv
import math

#Input is the output of LDHelmet or Pyrho for 1 scaffold, after filtering as described in the supplementary materials. 

file = "LDHelmet_output_filtered_PerScaffold.bed" #File name
chrom = "scaffold" #Scaffold ID
size = 1000 #Size of bin for map

kb_1_rates={}
kb_1_rates[chrom] = {}
estimates=pd.read_csv(file, sep="	", header=None)
#For each row in the output
for row in range(0, len(estimates.iloc[:,0])): 
    #Find the start and ending base
    start= estimates.iloc[row, 0]
    end = estimates.iloc[row, 1]
    #Iterate over those bases 
    for base, index in zip(range(start,end), range(0, len(range(start,end)))): 
        #Binning by size, add the estimated rate at that base and keep track of the number of bases in the bin.
        if math.floor(base/size) in kb_1_rates[chrom]: 
            kb_1_rates[chrom][math.floor(base/size)][0] += estimates.iloc[row, 2]
            kb_1_rates[chrom][math.floor(base/size)][1] += 1
        else: 
            kb_1_rates[chrom][math.floor(base/size)] = [estimates.iloc[row, 2], 1]

#To identify hotspots, scan the genome in 2kb windows overlapping by 1kb. The below code makes these bins out of the 1kb bins above. 

kb_2_rates={}
for chrom in kb_1_rates.keys():
    kb_2_rates[chrom]={}
    for key in sorted(kb_1_rates[chrom].keys()):
        #If it's the lowest bin, just add it to the next highest bin
        if key == min(kb_1_rates[chrom].keys()): 
            kb_2_rates[chrom][key+1] = [kb_1_rates[chrom][key][0], kb_1_rates[chrom][key][1]]
        #If it's the highest bin, add it to the last bin.
        elif key == max(kb_1_rates[chrom].keys()): 
            kb_2_rates[chrom][key][0] += kb_1_rates[chrom][key][0]
            kb_2_rates[chrom][key][1] += kb_1_rates[chrom][key][1]
        #Otherwise, add it to this bin and the next one.     
        else: 
            if key in kb_2_rates[chrom]:
                kb_2_rates[chrom][key+1] = [kb_1_rates[chrom][key][0], kb_1_rates[chrom][key][1]]
                kb_2_rates[chrom][key][0] += kb_1_rates[chrom][key][0]
                kb_2_rates[chrom][key][1] += kb_1_rates[chrom][key][1]
            else:
                kb_2_rates[chrom][key+1] = [kb_1_rates[chrom][key][0], kb_1_rates[chrom][key][1]]
                kb_2_rates[chrom][key] = [kb_1_rates[chrom][key][0], kb_1_rates[chrom][key][1]]

#Store the background recombination rate for each 2kb window in kb_2_rates
bgrates={}
total = 0
bgnum = 21 #Set how wide the background region is ((bgnum-1)x1000 = region). Here, it's 20kb on either side. 
bgrates[chrom]={}
#If there are fewer than 10 bins, skip this scaffold.
if len(kb_2_rates[chrom]) < 10: continue
minKey = min(kb_2_rates[chrom].keys())
maxKey = max(kb_2_rates[chrom].keys())
for key in kb_2_rates[chrom]: 
    #If there's less than 1/4 of the bases in this 2kb bin after filtering, skip it.
    if kb_2_rates[chrom][key][1] < 500: continue
    bg = 0
    #count keeps track of the number of missing 1kb bins.
    count = 0
    nope = 0
    for newkey in range(key-bgnum,key+bgnum+2): 
        if newkey not in kb_1_rates[chrom].keys(): 
            count += 1
    # If more than 4 background bins are completely missing, skip this bin.
    if count > 4: continue
    count = 0
    #For each 1kb region in the background region
    for newkey in range(key-bgnum,key+bgnum+2):
        #Except the central region +/- 1kb
        if newkey in range(key-1, key+2): continue
        #Store the sum of the recombination rates per bp and the total number of of bases. 
        if (newkey in kb_1_rates[chrom] and newkey in kb_2_rates[chrom]):
            bg += kb_1_rates[chrom][newkey][0]  
            count += kb_1_rates[chrom][newkey][1]
    #If missing more than half of the bases in the background, skip it. 
    if count > 20000: 
        #save the average recombination rate in the background and the recombination rate in the focal region. 
        bgrates[chrom][key] = [bg/count, kb_2_rates[chrom][key][0]/kb_2_rates[chrom][key][1]]
        total += 1


#Find hotspots. 
combdis = {}
for heat in [5]: #Set heat values.
    with open(hotspots.bed, "w") as f:
        f.write("")
    size=1000 #Set size of bin
    #For each scaffold
    for chrom in bgrates.keys(): 
        hotspotloc = []
        strength = []
        relstrength=[]
        distance = []
        maxstrengths = []
        maxrelstrengths = []
        #If the central 2kb has a heat >5x the background heat, save the hotspot location, relative heat, and recombination rate.
        for element in bgrates[chrom]: 
            location = element*size
            if bgrates[chrom][element][1] > heat*bgrates[chrom][element][0]: 
                hotspotloc.append(location)
                relstrength.append(bgrates[chrom][element][1]/bgrates[chrom][element][0])
                strength.append(bgrates[chrom][element][1])
        cleaned_hotspots = []
        donotadd = []
        #For each hotspot identified in the scaffold
        for key, i in zip(hotspotloc, range(0,len(hotspotloc))):
            #If the hotspot hasn't already been identified as a duplicate, add it. 
            if i not in donotadd:
                hotspotlocs = [key]
                strengths = [strength[i]]
                relstrengths = [relstrength[i]]
                maxloc = hotspotloc[i] 
                #For hotspots that are further along the chromosome than the central one
                for j in range(i+1,len(hotspotloc)): 
                    #If the distance is less than 5000
                    #if they are within 5kb of the focal hotspot, keep track of their stats so that they can be merged with the focal hotspot.
                    if abs(maxloc-hotspotloc[j])<5000: 
                        if abs(maxloc-hotspotloc[j]) in combdis: 
                            combdis[abs(maxloc-hotspotloc[j])] += 1
                        else: combdis[abs(maxloc-hotspotloc[j])] = 1
                        #If the location is already in the dictionary, append another location
                        maxloc=hotspotloc[j]
                        hotspotlocs.append(hotspotloc[j])
                        strengths.append(strength[j])
                        relstrengths.append(relstrength[j])
                        donotadd.append(j)
                #Merge hotspots (and stats) within 5kb of eachother.
                total=sum(hotspotlocs)
                count=len(hotspotlocs)
                cleaned_hotspots.append(total/count)
                #Also save the distance between min and max called hotspots
                distance.append(max(hotspotlocs)-min(hotspotlocs))
                maxstrengths.append(max(strengths))
                maxrelstrengths.append(max(relstrengths))
        #Output all hotspots after merging.
        with open(hotspots.bed), "a") as f:
            for i in range(0, len(cleaned_hotspots)):
                f.write(chrom.split("_bp")[0] + "\t" + str(int(cleaned_hotspots[i])-size-int(distance[i]/2)) + "\t" + str(int(cleaned_hotspots[i])+size+int(distance[i]/2)+1) + "\t" + str(distance[i]+2*size) + "\t" + str(maxstrengths[i])+ "\t" + str(maxrelstrengths[i])+"\n")