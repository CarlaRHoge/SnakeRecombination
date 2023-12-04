#!/usr/bin/env python

pref = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"

with open("sras.txt", "r") as fh:
    for line in fh:
        sra = line.strip()
        print(pref, 
              sra[:6], 
              "00{}".format(sra[-1]), 
              sra, 
              "{}_1.fastq.gz".format(sra), 
              sep="/")

        print(pref, 
              sra[:6], 
              "00{}".format(sra[-1]), 
              sra, 
              "{}_2.fastq.gz".format(sra), 
              sep="/")
    