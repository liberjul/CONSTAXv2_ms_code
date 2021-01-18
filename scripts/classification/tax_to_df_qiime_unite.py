# Sort qiime output taxonomy file by RDP taxonomy file
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

for k in range(5):
    for r in ["gen", "fam"]:
        for reg in ["", "_itsx.ITS1", "_itsx.ITS2"]:
            print(k, r, reg)
            with open(F"{args.dir}/uni_parts_{k}_{r}_qiime_query{reg}_taxonomy.tsv", "r") as ifile:
                buffer = "OTU_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
                line = ifile.readline()
                line = ifile.readline()
                while line != "":
                    tax_string = line.split("\t")[1]
                    tax_list1 = tax_string.strip().split(";")
                    tax_list2 = [x[3:] for x in tax_list1]
                    i = 0
                    tax_list3 = [""]*7
                    taxa = tax_list2[i]
                    while "unclassified" not in taxa:
                        tax_list3[i] = taxa
                        i += 1
                        if i < len(tax_list2):
                            taxa = tax_list2[i]
                        else:
                            break
                    tax_out = '\t'.join(tax_list3)
                    buffer = F"{buffer}{line.split('|')[1]}\t{tax_out}\n"
                    line = ifile.readline()
            with open(F"{args.dir}/uni_parts_{k}_{r}_qiime_query{reg}_taxonomy_comb.txt", "w") as ofile:
                ofile.write(buffer)
