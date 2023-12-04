#!/usr/bin/env python

import sys

for line in sys.stdin:
    fields = line.split()
    if not line.startswith("#"):
        geno_flds = fields[9:]
        geno = [g for g in geno_flds]
        
        