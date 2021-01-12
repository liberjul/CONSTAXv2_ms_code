def match_score(seq_1, seq_2):
    '''
    Function to match two strings and return number of substitutions.

    Assumes seq_1 can be degenerate, while seq_2 is not.
    '''
    # Dictionary of the possible non-mismatch letters for a primer base.
    degen_bases = {"A" : "ARWMDHVN", "C" : "CYSMBHVN", "G" : "GRSKBDVN", "T" : "TYWKBDHN",
                    "R" : "AGRDVN", "Y" : "CTYBHN", "S" : "GCSBVN", "W" : "ATWDHN",
                    "K" : "GTKBDN", "M" : "ACMHVN", "B" : "CGTBSYKN", "D" : "AGTDRWKN",
                    "H" : "ACTHMWYN", "V" : "ACGVMRSN", "N" : "ACTGRYSWKMBDHVN"}
    seq_1 = seq_1.upper()
    seq_2 = seq_2.upper()

    mismatch = 0
    for i in range(len(seq_1)):
        if seq_2[i] not in degen_bases[seq_1[i]]: # Degen and mismatch
            mismatch += 1
        else:
            continue
    return mismatch

def match_positions(primer, template):
    '''
    For a given primer seq, find the number of matches at the beginning position.
    '''
    scores = np.zeros(len(template) - len(primer))
    for i in range(len(template) - len(primer)):
        scores[i] = match_score(primer, template[i:i + len(primer)])
    return scores

def reverse_complement(seq):
    '''
    Returns a reverse complemented sequence string.
    '''
    comp_dict = {"A" : "T", "T" : "A", "G" : "C", "C" : "G", "N" : "N",
                "R" : "Y", "Y" : "R", "S" : "S", "W" : "W", "K" : "M",
                "M" : "K", "B" : "V", "V" : "B", "D" : "H", "H" : "D"}
    rev_seq = seq.upper()[::-1]
    rc_seq = ""
    for i in rev_seq:
        rc_seq = F"{rc_seq}{comp_dict[i]}"
    return rc_seq

def best_match(primer_forward, primer_reverse, template):
    '''
    Returns the index and score of the best match position.
    [forward index, forward score, reverse index, reverse score]
    '''
    m_f = match_positions(primer_forward, template)
    m_r = match_positions(primer_reverse, reverse_complement(template))[::-1]
    return [np.argmin(m_f), np.min(m_f), np.argmin(m_r), np.min(m_r)]

def pcr_prod(primer_forward, primer_reverse, template, max_mis = 3, include_primers = False):
    f_ind, f_score, r_ind, r_score = best_match(primer_forward, primer_reverse, template)
    if f_score > max_mis:
        print(F"Forward primer does not match below mismatch score. Score = {f_score}")
    if r_score > max_mis:
        print(F"Reverse primer does not match below mismatch score. Score = {r_score}")
    if f_score > max_mis or r_score > max_mis:
        return ""
    else:
        amplicon = template[f_ind:r_ind + len(primer_reverse)]
        if include_primers:
            return amplicon
        else:
            return amplicon[len(primer_forward) : -len(primer_reverse)]

def pcr_fastas(in_fasta, out_fasta, primer_forward, primer_reverse, max_mis = 3, include_primers = False, verbose = False):
    missed = 0
    with open(in_fasta, "r") as ifile, open(out_fasta, "w") as ofile:
        line = ifile.readline()
        while line != "":
            if line[0] == ">":
                header = line
                line = ifile.readline()
            else:
                template = ""
                while line != "" and line[0] != ">":
                    template += line.strip()
                    line = ifile.readline()
                template = template.replace("U", "T")
                amp = pcr_prod(primer_forward, primer_reverse, template, max_mis = max_mis, include_primers = include_primers)
                if len(amp) > 0:
                    ofile.write(F"{header}{amp}\n")
                else:
                    if verbose:
                        print(header.strip())
                    missed += 1
        print(F"{missed} records did not amplify")

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input FASTA file")
parser.add_argument("-o", "--output", type=str, help="output FASTA file")
parser.add_argument("-f", "--forward", type=str, help="forward primer sequence")
parser.add_argument("-r", "--reverse", type=str, help="reverse primer sequence")
parser.add_argument("-m", "--max_mis", type=int, help="maximum number of mismatches to include", default=3)
parser.add_argument("-p", "--include_primers", type=str, help="include primers?", default="False")
parser.add_argument("-v", "--verbose", type=str, help="print mismatches?", default="False")

args = parser.parse_args()

pcr_fastas(args.input, args.output, args.forward, args.reverse,
    max_mis=args.max_mis,
    include_primers=bool(args.include_primers),
    verbose=bool(args.verbose))
