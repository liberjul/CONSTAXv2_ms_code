import argparse
import pandas as pd
import numpy.random as npr
import numpy as np

def partition(df, level):
    df = df.astype('str')
    quer_ids, ref_ids = {}, {}
    uniq_tax = np.unique(df.iloc[:,level])
    print(len(uniq_tax))
    for taxon in uniq_tax:
        sub_df = df[df.iloc[:,level] == taxon]
        sub_level_taxa = sub_df.iloc[:,level + 1]
        sub_level_unique_1, counts = np.unique(sub_level_taxa, return_counts=True)
        sub_level_unique = sub_level_unique_1[sub_level_unique_1 != "-"]
        if len(sub_level_unique) > 1:
            counts = counts[sub_level_unique_1 != "-"] # Counts of named taxa
            quer_taxa = list(sub_level_unique[counts == 2]) # Query taxa included all with a count of two
            sub_level_unique = sub_level_unique[counts != 2]
            npr.shuffle(sub_level_unique)

            spl = int(len(sub_level_unique)/2) - len(quer_taxa)
            quer_taxa.extend(sub_level_unique[:spl])
            ref_taxa = sub_level_unique[spl:]
            for q in quer_taxa:
                for i in sub_df[sub_df.iloc[:,level + 1] == q][0]:
                    quer_ids[i] = None
            for r in ref_taxa:
                for i in sub_df[sub_df.iloc[:,level + 1] == r][0]:
                    ref_ids[i] = None
    return quer_ids, ref_ids


parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tax", type=str, help="header taxonomy file")
parser.add_argument("-d", "--db", type=str, help="database file input")
parser.add_argument("-f", "--format", type=str, help="database file format")
parser.add_argument("-q", "--query", type=str, help="query file output")
parser.add_argument("-r", "--ref", type=str, help="reference file output")
args = parser.parse_args()

tax = pd.read_table(args.tax, sep=";", header=None)

fam_quer, fam_ref = partition(tax, 5)
gen_quer, gen_ref = partition(tax, 6)

rec_dict = {}
with open(args.db, "r") as ifile:
    line = ifile.readline()
    while line != "":
        header = line
        if args.format == "SILVA":
            key = line.split('.')[0]
        else:
            key = F">{line.split('|')[1]}"
        line = ifile.readline()
        seq = ""
        while line != "" and line[0] != ">":
            seq += line.strip()
            line = ifile.readline()
        rec_dict[key] = F"{header}{seq}\n"
with open(args.query + "_fam.fasta", "w") as file:
    for key in fam_quer.keys():
        file.write(rec_dict[key])
with open(args.ref + "_fam.fasta", "w") as file:
    for key in fam_ref.keys():
        file.write(rec_dict[key])
with open(args.query + "_gen.fasta", "w") as file:
    for key in gen_quer.keys():
        file.write(rec_dict[key])
with open(args.ref + "_gen.fasta", "w") as file:
    for key in gen_ref.keys():
        file.write(rec_dict[key])
