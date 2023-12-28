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

recmap = "recombination_map.dict" #In 1kb bins, generated as in "hotspots.py"
allhots = pd.read_csv("hotspots.bed", sep="\t", header=None)
hotnucs = pd.read_csv(os.path.join("hotspots.nuc"), sep="\t") #output of bedtools.nuc for hotspot regions
chrom = "scaffold" #Scaffold ID
chromhots = allhots[allhots.iloc[:, 0] == chrom]
# Find mean rate per chromosome
count = 0
total = 0
for chromo in kb_1_rates.keys():
    for key in kb_1_rates[chromo]: 
        count+= kb_1_rates[chromo][key][0]
        total+= kb_1_rates[chromo][key][1]
halfmean = count / total / 2
#For every hotspot width, identify potential coldspots.
for wi in chromhots.iloc[:, 3].unique(): 
    outfile = "potential_coldspots.bed"
    with open(outfile, "w") as out: 
        out.write("#CHROM\tCOLDSTART\tCOLDEND\tHOTSTART\tHOTEND\tWIDTH\n")
    hots = chromhots[chromhots.iloc[:, 3] == wi]
    width = int(wi/1000)
    #Load GC content for given width and chromosome
    nucfile = "Potential_coldspots_nohots_" + str(wi) + ".nuc" #Output of "bedtools makewindows -g Genome.genome -w " + str(wi) + " -s 1000 | bedtools intersect -v -a stdin -b hotpots.bed | bedtools nuc -fi genome.fa -bed stdin > Potential_coldspots_noHotpots" + str(wi) + ".nuc"
    nucs = pd.read_csv(nucfile, sep="\t")
    nucchrom = nucs[nucs.iloc[:, 0] == chrom]
    hotnucchrom = hotnucs[hotnucs.iloc[:, 0] == chrom]
    # for each hotspot on the chrom of a given width
    for hotrow in range(len(hots.iloc[:, 0])): 
        coldcount = 0
        #Grab the GC content
        hotspotnuc = hotnucchrom[hotnucchrom.iloc[:, 1] == hots.iloc[hotrow, 1]].iloc[0, 7]
        # Find start and end of hotspot
        startkey = math.floor(hots.iloc[hotrow, 1]/1000)
        endkey = math.ceil(hots.iloc[hotrow, 2]/1000)
        #find background rates for 5mb surrounding a given hotspot
        winSize_rates={}
        for chrom in kb_1_rates.keys():
            #For each potential coldspot location (within 5mb of hotspot), find the rec rate.
            for key in range(startkey-5000, endkey+5000): 
                count = 0
                for i in range(key - math.floor(width/2), key + math.ceil(width/2)): 
                    count += 1
                    if i in kb_1_rates[chrom]: 
                        if key not in winSize_rates: winSize_rates[key] = [kb_1_rates[chrom][i][0], kb_1_rates[chrom][i][1]]
                        else: 
                            winSize_rates[key][0] += kb_1_rates[chrom][i][0]
                            winSize_rates[key][1] += kb_1_rates[chrom][i][1]
            with open(outfile, "a") as out: 
                for key in winSize_rates: 
                    #For each potential coldspot, continue if doesn't have at least half of the bases
                    if winSize_rates[key][1] < width/2: continue
                    # If the rate is less than 1/2 of the scaffold average
                    if winSize_rates[key][0] / winSize_rates[key][1] < halfmean:
                        # Find the GC content 
                        startcold = (key-math.floor(width/2))*1000
                        endcold = (key+math.ceil(width/2))*1000
                        nucspot = nucchrom[nucchrom.iloc[:, 1] == startcold]
                        if len(nucspot.iloc[:, 0]) != 1: continue
                        # If it's more than 1% different, continue
                        if abs(hotspotnuc - nucspot.iloc[0, 4]) > .01: continue
                        coldcount += 1
                        # otherwise, write it down. 
                        out.write(chrom + "\t" + str(startcold)+ "\t" +str(endcold)+ "\t"+ str(round(winSize_rates[key][0] / winSize_rates[key][1], 7)) + "\t"+ str(hots.iloc[hotrow, 1]) + "\t" + str(hots.iloc[hotrow, 2]-1) + "\t" + str(hots.iloc[hotrow, 3]) + "\n")
