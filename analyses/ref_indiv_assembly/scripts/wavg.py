import sys
import numpy as np

if not sys.stdin.isatty():
    rs = []
    ws = []
    for line in sys.stdin:
        line = line.strip()
        fields = line.split()
        rs.append(float(fields[-1]))  
        ws.append(float(fields[3])-float(fields[2]))
    if len(ws)==0:
        print(np.nan)
    else:
        print(np.average(rs, weights=ws))