
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


ZFgrp = {}
number = 1
#Assign the binding site associated with each ZF a letter.
#For each allele's position weight matrix
for file in glob("/moto/palab/users/crh2152/projects/CS_recomb/PRDM9_PWMs/*_PWM.txt"): 
    PWM = pd.read_csv(file, sep="   ", skiprows=1, header=None, index_col=0)
    #Identify each set of 3 columns, corresponding to the binding site of 1 ZF
    for i in range(int(len(PWM.iloc[0, :])/3)): 
        start = i * 3
        end = (i + 1) * 3
        currpwm = PWM.copy().iloc[:, start:end]
        currpwm.columns = range(1, 4)
        found = False
        #If the binding site has been identified before, continue
        for j in ZFgrp: 
            if ZFgrp[j].equals(currpwm): 
                found = True
                continue
        #Otherwise, assign it the next letter.
        if found == False: 
            let = chr(ord('@')+number)
            ZFgrp[let] = currpwm
            number += 1

#For each allele, output the letters associated with its binding array.
#For each position weight matrix
for file in sorted(glob("/moto/palab/users/crh2152/projects/CS_recomb/PRDM9_PWMs/*_PWM.txt")):
    PWM = pd.read_csv(file, sep="   ", skiprows=1, header=None, index_col=0)
    name = ""
    #For each zinc finger, find the letter associated with that binding site, and add that letter to name.
    for i in range(int(len(PWM.iloc[0, :])/3)): 
        start = i * 3
        end = (i + 1) * 3
        currpwm = PWM.copy().iloc[:, start:end]
        currpwm.columns = range(1, 4)
        for j in ZFgrp: 
            if ZFgrp[j].equals(currpwm): 
                name += j
    print(file.split("/")[-1].split("_PWM")[0], name)

#The letters and (rough) similarity between the position weight matricies that they represent were used to organize and color figure 1C. 