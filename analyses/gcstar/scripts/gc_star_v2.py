#!/usr/bin/env/python

import sys
from optparse import OptionParser
import numpy as np
import pandas as pd
from collections import Counter

def read_nwk_tree(nwk_string):
    '''Returns list of species in newick tree'''
    
    clean_nwk_string = nwk_string.replace("(","").replace(")","").replace(";","").replace(",","")
    species = clean_nwk_string.split(":1")
    return species
    
def print_sequences(seq_dict, start, end):
    '''Prints each sequences in the dictionary in one line'''
    
    for sp,seq in seq_dict.items():
        pad = 20 - len(sp)
        sys.stdout.write("{}\t{}\n".format(sp+" "*pad, seq[start:end]))
    
def pad_string(string, maxl):
    '''Adds extra white space to string until it is maxl long'''
    
    return string + " "*(maxl-len(string))

def read_maf(maf_file, species, ref):
    '''Returns dictionary of fasta sequences for each species in the MAF'''
    
    with open(maf_file,"r") if maf_file is not "-" else sys.stdin as maf:
        
        for line in maf:
            line = line.strip()
            fields = line.split()    
            
            # Process the header, take the name of 'species' in the maf file
            if line.startswith("#"):
                if line.startswith("# hal"):
                    sequences = {sp:"" for sp in species}
                    tmp_sequences = {sp:"" for sp in species}
                    branches_seen = []
                    
            # Process alignment lines, add to dictionary of sequences
            if line.startswith("s"):
                if fields[0]=="s":
                    branch = fields[1].split(".")[0]
                    if branch not in species:
                        continue
                    dna = fields[-1].upper()
                    dna = dna.replace("-","N")
                    if branch==ref:
                        ref_chunk = fields[-1]
                    tmp_sequences[branch] += dna.upper()
                    branches_seen.append(branch)
                        
            # Only keep track of sequences if species were represented only once in previous block
            if line=="":
                branch_count = Counter(branches_seen)
                for branch in species:
                    if branch not in branches_seen:
                        sequences[branch] += "N"*len(ref_chunk)
                    elif branch_count[branch]>1:
                        sequences[branch] += "N"*len(ref_chunk)
                    else:
                        sequences[branch] += tmp_sequences[branch]

                tmp_sequences = {sp:"" for sp in species}
                branches_seen = []

    return sequences

def count_derived_alleles(sequences, anc_sp):
    '''Counts number of derived alleles in branch'''
    
    derived_counts = {sp:0 for sp in sequences}
    for i,anc in enumerate(sequences[anc_sp]):
        for sp in sequences:
            if sequences[sp][i]!=anc:
                derived_counts[sp] += 1
    
    for sp,counts in derived_counts.items():
        sys.stderr.write("{}\t{}\n".format(pad_string(sp, 25), 
                                           counts))

def safe_division(n,d):
    '''Returns NaN if denominator is 0'''
    
    if d==0:
        return np.nan
    else:
        return n/float(d)
    
def compute_gcstar(sequences, anc):
    '''Computes GC* given an ancestral sequence and sliding window size'''
    
    descendant_sp = [sp for sp in sequences if sp!=anc]
    derived_changes = {sp:{"AT": 0,"GC": 0} for sp in descendant_sp}
    ancestral_states = {sp:{"AT": 0,"GC": 0} for sp in descendant_sp}
    nucs = ["A","T","G","C"]
    at = 0
    gc = 0
    analysed_bps = 0

    for i,anc in enumerate(sequences[anc]):
        
        sp2alleles = {sp:sequences[sp][i].upper() for sp in descendant_sp}
        obs_nonN = [a for a in sp2alleles.values() if a in nucs]
        
        # Only if all species are observed
        if len(obs_nonN)==len(descendant_sp):
            n_seen = len(set(obs_nonN))
            counts = Counter(obs_nonN)
            # If monomorphic
            if n_seen==1:
                kind_anc = "AT" if anc in "AT" else "GC"
                for sp in descendant_sp:
                    ancestral_states[sp][kind_anc] += 1
            # If only one derived mutation
            if n_seen==2: #and any(c==1 for n,c in counts.items()):
                #der = [n for n,c in counts.items() if c==1][0]
                der = [n for n,c in counts.items() if n!=anc][0]
                kind_anc = "AT" if anc in "AT" else "GC"
                kind_der = "AT" if der in "AT" else "GC"
                for sp,a in sp2alleles.items():
                    if a==der:
                        derived_changes[sp][kind_der] += 1
                        ancestral_states[sp][kind_anc] += 1
                    elif a==anc:
                        ancestral_states[sp][kind_anc] += 1
    
    # Compute GC*
    gcstar_sp = {}
    total_derived = {}
    total_analysed = {}
    for sp in derived_changes:
        gcstar = safe_division(safe_division(derived_changes[sp]["GC"], ancestral_states[sp]["AT"]),
                               safe_division(derived_changes[sp]["GC"], ancestral_states[sp]["AT"]) + safe_division(derived_changes[sp]["AT"], ancestral_states[sp]["GC"]))
        gcstar_sp[sp] = gcstar
        total_derived[sp] = sum(derived_changes[sp].values())
        total_analysed[sp] = sum(ancestral_states[sp].values())
    #sys.stderr.write("{} analysed base-pairs".format(analysed_bps))
    return gcstar_sp, total_derived, total_analysed

def compute_gcstar_parsimony(sequences):
    '''Computes GC* given an ancestral sequence and sliding window size'''
    
    descendant_sp = [sp for sp in sequences if sp!="None"]
    derived_changes = {sp:{"AT": 0,"GC": 0} for sp in descendant_sp}
    ancestral_states = {sp:{"AT": 0,"GC": 0} for sp in descendant_sp}
    nucs = ["A","T","G","C"]
    at = 0
    gc = 0
    analysed_bps = 0

    for i,anyn in enumerate(sequences[descendant_sp[0]]):
        
        sp2alleles = {sp:sequences[sp][i].upper() for sp in descendant_sp}
        obs_nonN = [a for a in sp2alleles.values() if a in nucs]
        
        # Only if all species are observed
        if len(obs_nonN)==len(descendant_sp):
            n_seen = len(set(obs_nonN))
            counts = Counter(obs_nonN)
            # If monomorphic
            if n_seen==1:
                anc = obs_nonN[0]
                kind_anc = "AT" if anc in "AT" else "GC"
                for sp in descendant_sp:
                    ancestral_states[sp][kind_anc] += 1
            # If only one derived mutation
            if n_seen==2 and any(c==1 for n,c in counts.items()):
                anc = [n for n,c in counts.items() if c>1][0]
                der = [n for n,c in counts.items() if c==1][0]
                kind_anc = "AT" if anc in "AT" else "GC"
                kind_der = "AT" if der in "AT" else "GC"
                for sp,a in sp2alleles.items():
                    if a==der:
                        derived_changes[sp][kind_der] += 1
                        ancestral_states[sp][kind_anc] += 1
                    elif a==anc:
                        ancestral_states[sp][kind_anc] += 1
    
    # Compute GC*
    gcstar_sp = {}
    total_derived = {}
    total_analysed = {}
    for sp in derived_changes:
        gcstar = safe_division(safe_division(derived_changes[sp]["GC"], ancestral_states[sp]["AT"]),
                               safe_division(derived_changes[sp]["GC"], ancestral_states[sp]["AT"]) + safe_division(derived_changes[sp]["AT"], ancestral_states[sp]["GC"]))
        gcstar_sp[sp] = gcstar
        total_derived[sp] = sum(derived_changes[sp].values())
        total_analysed[sp] = sum(ancestral_states[sp].values())
    #sys.stderr.write("{} analysed base-pairs".format(analysed_bps))
    return gcstar_sp, total_derived, total_analysed
        
def replace_CpGs_gaps_by_Ns(sequences):
    '''Substitutes CpGs in any sequence by Ns'''
    
    return {sp:sequences[sp].upper().replace("CG","NN").replace("-","N") for sp in sequences}


def dict2df(dictionary, cols):
    '''Returns dataframe with keys as fist col and values as second'''
    
    df = pd.DataFrame([dictionary.keys(), dictionary.values()]).transpose()
    df.columns = cols
    return df

def extract_vars(fields):
    
    scaf = fields[0].decode('UTF-8')
    start = int(fields[1])
    end = int(fields[2])
    hotspot_width = int(fields[3])
    rrate = float(fields[4])
    heat = float(fields[5])
    cpgi_distance = float(fields[6])

    return scaf, start, end, hotspot_width, rrate, heat, cpgi_distance

def main():
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--scaf", action="store", type="str", dest="scaf", help="Scaffold name")
    parser.add_option("-s", "--start", action="store", type="int", dest="start", help="Start coordinate")
    parser.add_option("-e", "--end", action="store", type="int", dest="end", help="End coordinate")
    parser.add_option("-r", "--ref", action="store", type="str", dest="ref", help="Reference genome")
    parser.add_option("-d", "--dist", action="store", type="str", dest="distances", help="List of distances")
    parser.add_option("-f", "--file_list", action="store", type="str", dest="sp_list", help="List of species")
    parser.add_option("-m", "--maf", action="store", type="str", dest="maf", help="Maf file, - if stdin")
    parser.add_option("-l", "--lift", action="store", type="str", dest="lift", help="Lift base-pairs")
    parser.add_option("-a", "--anc", action="store", type="str", dest="anc", help="Ancestor branch", default="None")
    parser.add_option("-n", "--nature", action="store", type="str", dest="nature", help="Hotspots or coldspot?")
    parser.add_option("-p", "--distance_to_other_sp", action="store", type="str", dest="dist_other_sp", help="Distance to other species spot")
    (options, args) = parser.parse_args()
    
    # Which species to include
    with open(options.sp_list,"r") as species_file:
        species = [line.strip() for line in species_file]

    if options.anc not in species:
        species.append(options.anc)
        
    # Read alignment
    sequences = read_maf(options.maf, species, options.ref)
    masked_sequences = sequences
    #masked_sequences = replace_CpGs_gaps_by_Ns(sequences)
    
    # Compute GC star for region
    results = []
    slide = 500
    maxl = len(masked_sequences[species[0]])-1

    for slide_start in range(0, maxl, 100):
        window_sequences = {sp:masked_sequences[sp][slide_start:slide_start+slide-1] for sp in masked_sequences}
        if options.anc!="None":
            gc_star,total_derived,total_analysed = compute_gcstar(window_sequences, options.anc)
        else:
            gc_star,total_derived,total_analysed = compute_gcstar_parsimony(window_sequences)
            
        df_star = dict2df(gc_star, ["species", "gcstar"]).set_index("species").join(dict2df(total_derived, ["species", "total_derived"]).set_index("species")).join(dict2df(total_analysed, ["species", "bps"]).set_index("species")).reset_index()
        df_star["position"] = slide_start-(maxl/2)
        df_star["distances"] = options.distances
        df_star["hotspot"] = "{}:{}-{}".format(options.scaf, options.start, options.end)
        df_star["nature"] = options.nature
        df_star["lift_bp"] = options.lift
        df_star["dist_bp"] = options.dist_other_sp
        results.append(df_star)

        results_df = pd.concat(results)

    print(results_df.to_string(index=False, header=False))

if __name__ == "__main__":
    main()
