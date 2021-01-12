import numpy as np
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--in_file", type=str, default="/mnt/research/common-data/Bio/UserDownloads/CONSTAX/DB/sh_general_release_fungi_35077_RepS_04.02.2020.fasta",help="Reference FASTA to sample")
parser.add_argument("-q", "--quer_len", type=int, help="Number of sequences to sample for query")
parser.add_argument("-r", "--ref_len", type=int, help="Number of sequences to sample for reference")
parser.add_argument("-k", "--k_iter", type=int, help="Number of separate sampling runs")
parser.add_argument("--qname", type=str, default="query_iter", help="Query file name")
parser.add_argument("--rname", type=str, default="ref_iter" help="Reference file name")
args = parser.parse_args()

seq_dict = {}

with open(args.in_file, "r") as ifile:
    line = ifile.readline()
    count = 0
    while line != "":
        header = line[1:].strip()
        count += 1
        line = ifile.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq += line.strip()
            line = ifile.readline()
        if header in seq_dict.keys():
            if len(seq) > len(seq_dict[header]): # Use the longer sequence
                seq_dict[header] = seq
        else:
            seq_dict[header] = seq

keys = np.array(list(seq_dict.keys()))
print("Number of seqs", len(keys))

for k in range(args.k_iter):
    print(F"{k} iteration...")
    np.random.shuffle(keys)
    export_keys = keys[:args.quer_len]
    print("Query:", len(export_keys))
    ref_keys = keys[args.quer_len:(args.quer_len + args.ref_len)]
    print("Reference:", len(ref_keys))
    buffer = ""
    for key in export_keys:
        buffer += F">{key}\n{seq_dict[key]}\n"
    with open(F"{args.qname}{k}.fasta", "w") as ofile:
        ofile.write(buffer)
    buffer = ""
    for key in ref_keys:
        buffer += F">{key}\n{seq_dict[key]}\n"
    with open(F"{args.rname}{k}.fasta", "w") as ofile:
        ofile.write(buffer)
