from glob import glob
import os
import subprocess
import pandas as pd
import _pickle as pickle
import numpy as np
import random
import statistics
import csv
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import pylab 
import math
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
from scipy import stats
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn import metrics
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, precision_recall_curve, average_precision_score
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from statsmodels.regression import linear_model
from sklearn.preprocessing import PolynomialFeatures
%matplotlib inline

vcf_folder = "/moto/palab/users/crh2152/projects/CS_vcfs/"
rec_folder = "/moto/palab/users/crh2152/projects/CS_recomb/"
folder = "/moto/palab/users/crh2152/projects/CS_recomb/Model"
script_folder = "/moto/palab/users/crh2152/Scripts/"
reference = os.path.join(vcf_folder, "NewGenome", "Pantherophis_guttatus_masked_01092021.fa")
log_folder = os.path.join(rec_folder, "Log")

#Model on whole genome
Test = pd.read_csv("All_wLog_wDis_thinned.bed", sep="\t")
X_CP = Test.copy()[["PRDM9", "CpG", "TSS", "CpG:TSS", "P9:CpGOrTSS","GC", "LogRate_1mb", "DistoTelo"]]
y = Test["LogRate"]
X_CP = sm.add_constant(X_CP) 
est = sm.OLS(y, X_CP, missing="drop").fit() 
est.summary()

#Model for only scaffolds on macrochromosomes
Test = pd.read_csv("All_wLog_large_wDis_thinned.bed", sep="\t")
X_CP = Test.copy()[["PRDM9", "CpG", "TSS", "CpG:TSS", "P9:CpGOrTSS","GC", "LogRate_1mb", "DistoTelo"]]
y = Test["LogRate"]
X_CP = sm.add_constant(X_CP) 
l_est = sm.OLS(y, X_CP, missing="drop").fit() 
l_est.summary()

#Model for only scaffolds on microchromosomes
Test = pd.read_csv(os.path.join(folder, "All_wLog_small_wDis_thinned.bed"), sep="\t")
X_CP = Test.copy()[["PRDM9", "CpG", "TSS", "CpG:TSS", "P9:CpGOrTSS","GC", "LogRate_1mb", "DistoTelo"]]
y = Test["LogRate"]
X_CP = sm.add_constant(X_CP) 
s_est = sm.OLS(y, X_CP, missing="drop").fit() 
s_est.summary()