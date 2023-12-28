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
import matplotlib.pyplot as plt
from scipy import stats

#Identify sections of the position weight matricies of lengths corresponding to 5 to 13 zinc fingers, which are shared between alleles. 
pieces = {}
counter = {}
maxes = {}
#For each length 5 to 13 zinc fingers
for k in range(5, 13): 
    pieces[k] = {}
    counter[k] = {}
    maxes[k] = {}
    count = 0
    for file in glob("*_PWM.txt"): 
        PWM = pd.read_csv(file, sep="   ", skiprows=1, header=None, index_col=0)
        #If the motif is shorter than the length of the necessary binding site, skip.
        if len(PWM.iloc[0, :]) < 3*k: continue
        #Moving over the position weight matrix (by 3s, as every 3 rows = 1 ZF)
        for i in range(int(len(PWM.iloc[0, :])/3) - k+1): 
            start = i * 3
            end = (i + k) * 3
            currPWM = PWM.copy().iloc[:, start:end]
            #Renumber the columns starting at 0.
            currPWM.columns=range(len(currPWM.iloc[0, :]))
            found = False
            #If this segment is already in pieces, add 1 observation
            for df in pieces[k]: 
                if pieces[k][df].equals(currPWM): 
                    counter[k][df] += 1
                    found = True
            #Otherwise, add the observation to pieces (observed once)
            if found == False: 
                pieces[k][count] = currPWM
                counter[k][count] = 1
                count += 1
    #Find the maximum number of observations
    maxcounter = 0
    for i in counter[k]: 
        if counter[k][i] > maxcounter: maxcounter = counter[k][i]
    count = 0
    #Find the PWMs associated with the max number of observations of that size.
    for i in counter[k]: 
        if counter[k][i] == maxcounter: 
            maxes[k][str(count) + "_" + str(maxcounter)] = pieces[k][i]
            count += 1

# Convert from PWM representation to letter representation (one letter per ZF binding site)
max_lets = {}
for k in maxes: 
    max_lets[k] ={}
    for j in maxes[k]: 
        PWM = maxes[k][j]
        name = ""
        for i in range(int(len(PWM.iloc[0, :])/3)): 
            start = i * 3
            end = (i + 1) * 3
            currpwm = PWM.copy().iloc[:, start:end]
            currpwm.columns = range(1, 4)
            for z in ZFgrp: 
                if ZFgrp[z].equals(currpwm): 
                    name += z
        max_lets[k][j] = name
for k in max_lets: 
    print("Kmer ", k)
    for i in max_lets[k]: 
        print(max_lets[k][i], " Count ", i.split("_")[-1])

#Remove PWM pieces which are entirely contained within larger pieces and found the same number of times.
max_lets = {}
for k in maxes: 
    max_lets[k] ={}
    for j in maxes[k]: 
        PWM = maxes[k][j]
        name = ""
        for i in range(int(len(PWM.iloc[0, :])/3)): 
            start = i * 3
            end = (i + 1) * 3
            currpwm = PWM.copy().iloc[:, start:end]
            currpwm.columns = range(1, 4)
            for z in ZFgrp: 
                if ZFgrp[z].equals(currpwm): 
                    name += z
        max_lets[k][j] = name
for k in sorted(max_lets, reverse=True): 
    for seq in max_lets[k]: 
        for j in range(k-1, 5, -1): 
            for seq2 in list(max_lets[j]): 
                if seq2 not in maxes[j]: continue
                if max_lets[j][seq2] in max_lets[k][seq]: 
                    count1 = int(seq2.split("_")[-1])
                    count2 = int(seq.split("_")[-1])
                    if count1 == count2: 
                        maxes[j].pop(seq2)
                        max_lets[j].pop(seq2)
#Output as PWM files.
for k in maxes: 
    count = 0
    for i in maxes[k]: 
        outfile = os.path.join(folder, "PRDM9_PWMs", "KmerZFs_" + str(k) + "_" + str(count) + "_PWM.txt")
        with open(outfile, "w") as f: 
            f.write("k" + str(k) + "_count" + i.split("_")[-1] + "_" + str(count) + "\n")
        count += 1
        maxes[k][i].to_csv(outfile, sep="\t", mode="a", header=False)