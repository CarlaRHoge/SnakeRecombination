#!/usr/bin/env python


def parse_motifs(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    motifs = []
    motif = None
    new_motif = ''
    processing_matrix = False
    for line in lines:
        if line.startswith('MOTIF'):
            motif = line.split()[1].split("-")[-1]
            processing_matrix = False
            new_motif = ''
        elif line.startswith('letter-probability matrix'):
            processing_matrix = True
        elif processing_matrix:
            probs = list(map(float, line.split()))
            max_prob = max(probs)
            if max_prob < 0.8:
                new_motif += 'N'
            else:
                max_prob_index = probs.index(max_prob)
                if max_prob_index == 0:
                    new_motif += 'A'
                elif max_prob_index == 1:
                    new_motif += 'C'
                elif max_prob_index == 2:
                    new_motif += 'G'
                elif max_prob_index == 3:
                    new_motif += 'T'
        if motif and new_motif and len(new_motif) == len(motif):
            motifs.append(new_motif)
            motif = None
            new_motif = ''
            processing_matrix = False

    return motifs

filename = "/moto/palab/users/crh2152/projects/CS_recomb/LDHelmet/Hotspots/Streme_NoPromoters/streme_p05.txt"
#filename = "/moto/palab/users/crh2152/projects/CS_recomb/LDHelmet/Hotspots/Streme_NoPromoters/streme_p10.txt"
motifs = parse_motifs(filename)
for motif in motifs:
    print(motif.rstrip("N").lstrip("N"))

