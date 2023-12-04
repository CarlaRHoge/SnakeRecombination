#!/usr/bin/env python

import sys

#['Super_scaffold_2', '2716000', '2718000', '0.0826659750000002', '5.157762315983274', '6151000', '6153000', '0.0085208', '2000']

for line in sys.stdin:
    line = line.strip()
    fields = line.split('\t')
    hot_scaf = fields[0]
    hot_start = int(fields[1])
    hot_end = int(fields[2])
    hot_rate = float(fields[3])
    hot_heat = float(fields[4])

    cold_start = int(fields[5])
    cold_end = int(fields[6])
    cold_rate = float(fields[7])
    l = int(fields[8])

    hot_id = "{}:{}-{}".format(hot_scaf, hot_start, hot_end)
    cold_id = "{}:{}-{}".format(hot_scaf, cold_start, cold_end)
    
    print(hot_id, hot_rate/hot_heat, "hotspot", sep='\t')
    print(cold_id, cold_rate, "coldspot", sep='\t')
    