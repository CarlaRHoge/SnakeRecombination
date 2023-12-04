#!/usr/bin/env python

from optparse import OptionParser
import itertools
import gzip
import sys

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-v", "--vcf", action="store", type="str", dest="vcf",help="VCF file")
parser.add_option("-o", "--out", action="store", type="str", dest="outputsuff",help="Output suffix")
parser.add_option("-f", "--family", action="store", type="str", dest="family",help="Family ped file")
(options, args) = parser.parse_args()

# Read families
with open(options.family, "r") as fh:
    family, parents = {}, {}
    i2f = {}
    i2sex = {}
    for line in fh:
        if line.startswith("#"):
            continue
        line = line.strip()
        fields = line.split()
        family_id, kid, father, mother, sex, phenotype = fields
        i2sex[father] = "1"
        i2sex[mother] = "2"
        i2sex[kid] = "0"
        for indiv in [kid, father, mother]:
            i2f[indiv] = family_id
        if family_id not in family:
            family[family_id] = {}
        family[family_id][kid] = [father, mother]
        parents[kid] = [father, mother]

# Read VCF
#with open(sys.stdin, "r") as vcf_file:
with open(options.vcf, 'r') if options.vcf is not "-" else sys.stdin as vcf_file:

    for i,line in enumerate(vcf_file):
        line = line.strip()#.decode("utf-8").strip()
        fields = line.split()
        
        # Header & open hapi files
        if line.startswith("#"): 
            if line.startswith("#CHROM"):
                samples = fields[9:]
                s2i = {s:i for i,s in enumerate(samples)}
                header_size = i
                flat_recomb = 1/1e6
                
                #Open output files
                map_file = open("{}.map.txt".format(options.outputsuff), "w")
                marker_file = open("{}.marker.txt".format(options.outputsuff), "w")
                indiv2family = {s:k for s in samples for k,v in family.items() if s in v}

                # Write first lines in geno file
                geno_list = []
                geno_list.append([i2f[s] for s in samples])
                geno_list.append([s for i,s in enumerate(samples)])
                geno_list.append([parents[s][0] if s in indiv2family else "x" for s in samples])
                geno_list.append([parents[s][1] if s in indiv2family else "x" for s in samples])
                geno_list.append([i2sex[s] for i,s in enumerate(samples)])
                geno_list.append(["0" for s in samples])

            continue

        # Genotypes & fill hapi lists
        else:
            jump = False
            chrom, pos = fields[0], int(fields[1])
            
            genotypes = fields[9:]
            sgenotypes = [g.split(":")[0].replace("|","/") if "./.:" not in g else "0" for g in genotypes]
            sgenotypes = ["0/0" if g=="0" else "/".join([str(n+1) for n in map(int,g.split("/"))]) if g!="." else "0/0" for g in sgenotypes]
            
            # All hets?
            if all(g==sgenotypes[0] for g in sgenotypes):
                continue

            # Filter sites not consistent with Mendelian segregation
            for kid in parents:
                ikid = s2i[kid]
                ifather = s2i[parents[kid][0]]
                imother = s2i[parents[kid][1]]
                
                k = tuple(sgenotypes[ikid].split("/"))
                # If kid missing jump
                if "0" in k:
                    continue

                p1 = sgenotypes[ifather].split("/")
                p2 = sgenotypes[imother].split("/")
                # If any parent is missing
                if "0" in p1 or "0" in p2:
                    jump = True
                    break
                # If both parents homozygous or both heterozygous, or no heterozygous in any
                if sorted(p1)==sorted(p2) or not any([sorted(pg)==["1","2"] for pg in [p1,p2]]):
                    jump = True
                    break
                # Check segregation
                if k not in list(itertools.product(p1,p2)) + list(itertools.product(p2,p1)):
                    jump = True
                    break
                
            if jump:
                continue

            geno_list.append(sgenotypes)
            snp_id = "{}_{}".format(chrom, pos)
            map_file.write("1\t{}\t{}\n".format(snp_id, pos*flat_recomb))
            marker_file.write("{}\t{}\n".format("M", snp_id))
            
# Output genotype file
geno_file = open("{}.geno.txt".format(options.outputsuff), "w")
tranposed_list = list(map(list, zip(*geno_list)))
divider = 6
for t in tranposed_list:
    geno_file.write("\t".join(t[:divider]) + "\t" + " ".join(t[divider:]) + "\n")

