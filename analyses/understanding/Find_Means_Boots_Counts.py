#!/usr/bin/env python

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

#The folder containing bed input_files resulting from "bedtools closest" where bed input_file "a" is the recombination map and bed input_file "b" is the features of interest
fold = sys.argv[1]
#Bin width of the recombination map - 100bp 
binwidth = int(sys.argv[2])
#Maximum distance to consider (i.e. only looking 20kb from feature of interest)
dismax = int(sys.argv[3])
#Returns a list of len num of bootstrapped means from values.
def bootstrapList(values, num):
    results = []
    for i in range(0, num): 
        currset = np.random.choice(values,size=len(values), replace=True)
        currmean = sum(currset)/len(values)
        results.append(currmean)
    return results

feature=fold.split("/")[-1]
print("Working!")
print(feature)
#Will output the mean for each distance, the number of data points in each distance bin, and the bootstraps for each distance. 
meansout = fold + "_all_Means_dict.txt"
bootout = meansout.split("_Means")[0] + "_Boots_dict.txt"
countout = meansout.split("_Means")[0] + "_Counts_dict.txt"
distofeat={}
#For each bed input_file
for input_file in glob(os.path.join(fold, "*.bed")):
    print(input_file)
    # Remove it if there isn't anything in it. 
    if os.stat(input_file).st_size == 0: 
        subprocess.call(["rm", input_file])
        continue
    #Read in the bed input_file
    rates = pd.read_csv(input_file, header=None, sep="	")
    #Drop rows with NAs 
    rates.dropna(inplace=True)
    #For each datapoint in the map
    for row in range(0, len(rates.iloc[:, 0])): 
        # Find the distance (binned by the bin width) to the nearest feature (in the last column)
        dis = int(round(rates.iloc[row, len(rates.iloc[0, :])-1]/binwidth) * binwidth)
        #If it's larger than the largest distance considered, keep moving
        if dis > dismax: continue
        # Save the recombination rate
        rate = rates.iloc[row, 3]
        # Add the rec rate to a nested dictionary, where the first key is the distance to the feature, and the second is the rec rate. Keep track of how many times one rec rate is seen at a given distance to feature. 
        if dis in distofeat:
            if rate in distofeat[dis]: 
                distofeat[dis][rate] += 1
            else: distofeat[dis][rate] = 1
        else: 
            distofeat[dis] = {}
            distofeat[dis][rate] = 1
mean = {}
boots = {}
counts = {}
#For each distance in the dictionary
for dis in distofeat: 
    counts[dis] = 0
    listrates = []
    for rate in distofeat[dis]: 
        #Keep track of how many data points in a given distance
        counts[dis] += distofeat[dis][rate]
        #Make a list of all recombination rates at a given distance
        for i in range(distofeat[dis][rate]): 
            listrates.append(rate)
    #Find the mean of the list of rec rates
    mean[dis] = sum(listrates)/len(listrates)
    #Bootstrap from that list
    boots[dis] = bootstrapList(listrates, 100)
with open(meansout, "wb") as f: 
    pickle.dump(mean, f)
with open(countout, "wb") as f: 
    pickle.dump(counts, f)
with open(bootout, "wb") as f: 
    pickle.dump(boots, f)
