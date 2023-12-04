#!/usr/bin/env/python

import sys
from optparse import OptionParser
from collections import Counter
from Bio import Seq
from Bio.SeqUtils import IUPACData
from Bio.SeqUtils import seq3
import re

def split_sequences(sequences):
    '''Splits sequences in dictionaries into chunks without any N'''

    result = {species: [[]] for species in sequences.keys()}
    sequence_length = len(next(iter(sequences.values())))

    for i in range(sequence_length):
        if any(sequence[i] == 'N' for sequence in sequences.values()):
            for species in sequences.keys():
                result[species].append([])
            continue
        for species, sequence in sequences.items():
            result[species][-1].append(sequence[i])

    for species in sequences.keys():
        result[species] = ["".join(subseq) for subseq in result[species] if subseq]

    return result


def read_maf(maf_file, species, ref):
    '''Returns dictionary of fasta sequences for each species in the MAF'''
    
    with open(maf_file,"r") if maf_file is not "-" else sys.stdin as maf:

        skipNext = False
        compare = [sp for sp in species if sp!=ref]
        results = {sp:{"loses":0, "gains":0} for sp in compare}
        
        for line in maf:
            line = line.strip()
            fields = line.split()    
            
            # Skip the single empty line after the header and before sequence blocks
            if skipNext:
                skipNext = False
                continue
            
            # Process the header, take the name of 'species' in the maf file
            if line.startswith("#"):
                if line.startswith("# hal"):
                    sequences = {sp:[] for sp in species}
                    tmp_sequences = {sp:"" for sp in species}
                    branches_seen = []
                    skipNext = True
                    
            # Process alignment lines, add to dictionary of sequences
            if line.startswith("s"):
                if fields[0]=="s":
                    branch = fields[1].split(".")[0]
                    if branch not in species: continue
                    dna = fields[-1].upper()
                    dna = dna.replace("-","N")
                    if branch==ref:
                        chrom = fields[1].split(".")[1]
                        chunkid = f"{chrom}:{fields[2]}"
                        ref_chunk = fields[-1]
                    tmp_sequences[branch] += dna.upper()
                    branches_seen.append(branch)
                        
            # Only keep track of sequences if species were represented only once in previous block
            if line=="":
                branch_count = Counter(branches_seen)
                if len(branch_count)==len(species) and all(v==1 for k,v in branch_count.items()):
                    split_seqs = split_sequences(tmp_sequences)
                    n_seqs = len(split_seqs[ref])

                    if n_seqs==0: continue
                    for i,_ in enumerate(split_seqs[species[0]]):
                        for sp in species:
                            print(f">{sp}|{chunkid}|{i}")
                            print(split_seqs[sp][i])
                                
                tmp_sequences = {sp:"" for sp in species}
                branches_seen = []

    return results

def main():
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file_list", action="store", type="str", dest="sp_list", help="List of species")
    parser.add_option("-m", "--maf", action="store", type="str", dest="maf", help="Maf file, - if stdin")
    parser.add_option("-r", "--ref", action="store", type="str", dest="ref", help="Reference genome")
    (options, args) = parser.parse_args()
    
    # Which species to include
    with open(options.sp_list,"r") as species_file:
        species = [line.strip() for line in species_file]
        
    # Read alignment and score gains and loses
    read_maf(options.maf, species, options.ref)

if __name__ == "__main__":
    main()
