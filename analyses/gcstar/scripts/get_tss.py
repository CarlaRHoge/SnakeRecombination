#!/usr/bin/env python

import sys

for line in sys.stdin:
    line = line.strip()
    fields = line.split()
    #scaffold-mi4    maker   gene    59160   72969   .       +
    strand = fields[6]
    if strand=="+":
        start = int(fields[3])-2000 
        start = 1 if start < 1 else start
        print(fields[0], start , fields[3], sep="\t")
    if strand=="-":
        print(fields[0], fields[4], int(fields[4])+2000, sep="\t")
    