#!/usr/bin/env python

import sys

seen = []
gcs = []
target = float(sys.argv[1])

for line in sys.stdin:
    fields = line.split()
    gc = float(fields[-1])/100
    gcs.append(abs(gc - target))
    seen.append(fields)
    
idx = gcs.index(min(gcs))
sel = seen[idx]
print(fields[0].replace(":","\t").replace("-", "\t"), 
      "{:.3f}".format(min(gcs)),
      sep="\t")


    