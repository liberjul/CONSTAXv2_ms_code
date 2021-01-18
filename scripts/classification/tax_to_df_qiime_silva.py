# Sort qiime output taxonomy file by RDP taxonomy file
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

for k in range(5):
    for r in ["gen", "fam"]:
        for reg in ["_full", "_v3-4", "_v4"]:
            print(k, r, reg)
            if reg == "_full":
                filename = F"{args.dir}/sil_parts_{k}_query_sub_1k_{r}_qiime_taxonomy.tsv"
            else:
                filename = F"{args.dir}/sil_parts_{k}_query_{r}{reg}_qiime_taxonomy.tsv"
            with open(filename, "r") as ifile:
                buffer = "OTU_ID\tRank_1\tRank_2\tRank_3\tRank_4\tRank_5\tRank_6\tRank_7\n"
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
                    buffer = F"{buffer}{line.split('.')[0]}\t{tax_out}\n"
                    line = ifile.readline()
            with open(F"{args.dir}/sil_parts_{k}_{r}_qiime_query{reg}_taxonomy_comb.txt", "w") as ofile:
                ofile.write(buffer)
