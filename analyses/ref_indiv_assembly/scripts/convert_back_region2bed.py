#!/usr/bin/env python

import sys

for line in sys.stdin:
    line = line.strip()
    fields = line.split()
    region = fields[1]
    scaf = region.split(":")[0]
    start = int(region.split(":")[1].split("-")[0])
    end   = int(region.split(":")[1].split("-")[1])
    print(scaf, 
          int(start-5e3), 
          int(end+5e3), 
          "\t".join(fields[2:]),
          fields[0],
          sep="\t")

