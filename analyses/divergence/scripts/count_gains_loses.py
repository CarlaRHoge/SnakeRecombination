#!/usr/bin/env/python

import sys
from optparse import OptionParser
from collections import Counter
from Bio import Seq
from Bio.SeqUtils import IUPACData
from Bio.SeqUtils import seq3
import re

# def find_motif(iupac_motif, sequence):
#     # Convert IUPAC motif to regular expression
#     regex_motif = "".join(IUPACData.ambiguous_dna_values[i] for i in iupac_motif.upper())
#     motif = Seq.Seq(regex_motif)
    
#     # Convert string sequence to BioPython sequence object
#     sequence = Seq.Seq(sequence)
    
#     # Find motif in sequence
#     motif_indices = []
#     motif_len = len(motif)
#     for i in range(len(sequence) - motif_len + 1):
#         if sequence[i:i+motif_len] == motif:
#             motif_indices.append(i)
    
#     return motif_indices

def find_motif(iupac_motif, sequence):
    
    regex_motif = "".join(i for i in iupac_motif.upper())
    regex_motif = regex_motif.replace('N', '.')
    sequence = Seq.Seq(sequence)
    
    # Find motif in sequence using regular expressions
    motif_indices = [m.start() for m in re.finditer('(?={})'.format(regex_motif), str(sequence))]
    
    return motif_indices


def reverse_complement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}
    complement_seq = "".join(complement_dict[base] for base in seq.upper())
    return complement_seq[::-1]

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

def read_maf(maf_file, species, ref, target):
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

                    for i in range(n_seqs):
                        fhits = {sp:find_motif(target, split_seqs[sp][i]) for sp in species}
                        rhits = {sp:find_motif(target, reverse_complement(split_seqs[sp][i])) for sp in species}
                        hits = {sp:fhits[sp]+[int(x+1e8) for x in rhits[sp]] for sp in species}
                        if all(v==hits[ref] for k,v in hits.items()): continue
                        if hits[compare[0]]==hits[compare[1]]: continue

                        # Gains
                        for c1 in compare:
                            c2 = [c for c in compare if c!=c1][0]
                            for h in hits[c1]:
                                if h not in hits[ref] and h not in hits[c2]:
                                    results[c1]["gains"] += 1
                        # Loses
                        for h in hits[ref]:
                            for c1 in compare:
                                c2 = [c for c in compare if c!=c1][0]
                                if h not in hits[c1] and h in hits[c2] and h in hits[ref]:
                                    results[c1]["loses"] += 1
#                        print(hits)
#                        print(results)
                                
                tmp_sequences = {sp:"" for sp in species}
                branches_seen = []

    return results

def main():
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file_list", action="store", type="str", dest="sp_list", help="List of species")
    parser.add_option("-m", "--maf", action="store", type="str", dest="maf", help="Maf file, - if stdin")
    parser.add_option("-t", "--target", action="store", type="str", dest="target", help="Target DNA motif")
    parser.add_option("-r", "--ref", action="store", type="str", dest="ref", help="Reference genome")
    (options, args) = parser.parse_args()
    
    # Which species to include
    with open(options.sp_list,"r") as species_file:
        species = [line.strip() for line in species_file]
        
    # Read alignment and score gains and loses
    results = read_maf(options.maf, species, options.ref, options.target)
    for sp,balance in results.items():
        for kind,instances in balance.items():
            print(f"{sp}\t{kind}\t{instances}")

if __name__ == "__main__":
    main()
