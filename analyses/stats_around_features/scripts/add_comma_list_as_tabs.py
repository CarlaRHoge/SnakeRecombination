#!/usr/bin/env python

import sys
from optparse import OptionParser

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-l", "--list", action="store", type="str", dest="additions", 
                  help="Comma separated list of strings to add tab separated") 
(options, args) = parser.parse_args()

addons = "\t".join(options.additions.split(","))

for line in sys.stdin:
    line = line.strip()
    print(line, addons, sep="\t")
