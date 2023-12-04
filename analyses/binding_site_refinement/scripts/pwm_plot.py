#/usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt
from optparse import OptionParser

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-o", "--pdf", action="store", type="str", dest="pdf")
(options, args) = parser.parse_args()

rows = []
for line in sys.stdin:
    if line.startswith(">"):
        continue
    fields = line.split()
    rows.append(list(map(float,fields[1:])))

df = pd.DataFrame(rows).transpose()
df.plot(kind='bar', stacked=True)

# labels for x & y axis
plt.xlabel('Position')
plt.ylabel('Nucleotide frequency')
plt.legend(labels=["A","C","T","G"], loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig(options.pdf)