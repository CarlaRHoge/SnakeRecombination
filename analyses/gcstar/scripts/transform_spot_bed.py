#!/usr/bin/env python

import sys

for line in sys.stdin:
    line = line.strip()
    fields = line.split()
    scaf = fields[0]
    slop_start = int(fields[1])
    slop_end   = int(fields[2])
    og_str = fields[3]
    start = int(og_str.split(":")[1].split("-")[0])
    end   = int(og_str.split(":")[1].split("-")[1])
    print(og_str,start,end)

