#!/usr/bin/env python

import sys

stats = {}

for line in sys.stdin:
    line = line.strip()
    fields = line.split()
    kind = fields[-1]
    if kind not in stats:
        stats[kind] = [0,0]
    stats[kind][0] += 1
    stats[kind][1] += int(fields[2]) - int(fields[1])

for k in stats:
    print("{}\t{}\t{}".format(stats[k][0], stats[k][1], k))

    