import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

tax_dict = {"0" : ""}
taxid_dict = {}
for k in range(1):
    for r in ["gen"]:
        with open(F"{args.dir}/sh_general_release_fungi_35077_RepS_04.02.2020__RDP_taxonomy_trained.txt", "r") as ifile:
            with open(F"./unite_test/data/unite_tax.txt", "w") as ofile:
                line = ifile.readline()
                line = ifile.readline()
                while line != "":
                    id, taxa, parent, rank, rank_name = line.strip().split("*")
                    par_string = tax_dict[parent]
                    new_string = F"{par_string}{taxa};"
                    tax_dict[id] = new_string
                    taxid_dict["Root;" + new_string.strip(";")] = id
                    ofile.write(F"{new_string}\t{id}\t{rank_name.lower()}\t138\n")
                    line = ifile.readline()
        with open(F"{args.dir}/sh_general_release_fungi_35077_RepS_04.02.2020__RDP_trained.fasta", "r") as ifile:
            with open(F"./unite_test/seqid2taxid.map", "w") as ofile:
                line = ifile.readline()
                while line != "":
                    if line[0] == ">":
                        seqid, tax_string = line[1:].strip().split("\t")
                        taxid = taxid_dict[tax_string]
                        ofile.write(F"{seqid}\t{taxid}\n")
                    line = ifile.readline()
