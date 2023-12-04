#!/usr/bin/env python

import sys
import pandas as pd

wd = "/moto/palab/projects/recombination_snakes/shared_corn_rattle/gcstar_zf"
misc_file = "{}/misc_files/zf2zf.txt".format(wd)
r2r = pd.read_csv(misc_file, sep="\t", header=None)
r2r.columns = ["og_scaf", "og_l", "t_scaf", "t_l", "offset"]
d_scaf = r2r.set_index("og_scaf")["t_scaf"].to_dict()
d_offs = r2r.set_index("t_scaf")["offset"].to_dict()

for line in sys.stdin:
    line = line.strip()
    # If in corn snake assembly
    if line.startswith("Super_scaffold"):
        print(line)
    # If crotalus
    else:
        fields = line.split()
        og_scaf = fields[0]
        start   = int(fields[1])
        end     = int(fields[2])
        t_scaf  = d_scaf[og_scaf]
        offset  = d_offs[t_scaf]
        print(t_scaf, start-offset, end-offset, "\t".join(fields[3:]), sep="\t")
    