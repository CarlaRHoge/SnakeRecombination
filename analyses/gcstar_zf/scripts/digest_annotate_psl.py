#!/usr/bin/env python

import sys
from optparse import OptionParser
import numpy as np

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-b", "--branch", action="store", type="str", dest="branch", help="Species")
parser.add_option("-s", "--spot", action="store", type="str", dest="spot", help="Spot kind") 
(options, args) = parser.parse_args()

spots = {}

# Read lifted coordinates from stdin, save the number of base-pairs that map to
# a spot of the same kind in the other species
for line in sys.stdin:
    line = line.strip()
    fields = line.split()
    start = int(fields[1])
    end = int(fields[2])
    spot_id = fields[3]
    frac = float(fields[4])
    dist = int(fields[5])
    if spot_id not in spots:
        spots[spot_id] = {"frac":[], "dist":[]}
    spots[spot_id]["frac"].append((end-start)*frac)
    spots[spot_id]["dist"].append(dist if dist>=0 else np.nan)

# Read the original coordinate bed file, annotate the number calculated above
with open("misc_files/{}.{}.bed".format(options.branch, options.spot), "r") as fh:
    for line in fh:
        line = line.strip()
        fields = line.split()
        spot_id = fields[4]
        mapped_bps = np.nan if spot_id not in spots else sum(spots[spot_id]["frac"])
        distance = np.nan if spot_id not in spots else min(spots[spot_id]["dist"])
        print(line, mapped_bps, distance, sep="\t")
    
   
    