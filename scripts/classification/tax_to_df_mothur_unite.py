import argparse
# Sort mothur output taxonomy file by RDP taxonomy file
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

for k in range(5):
    for r in ["gen", "fam"]:
        for c in ["wang", "knn"]:
            for reg in ["", "_itsx.ITS1", "_itsx.ITS2"]:
                print(k, r, c, reg)
                mothur_d = {}
                with open(F"{args.dir}/unite_partition_{k}_query_{r}{reg}.uni_parts_{k}_{r}.{c}.taxonomy", "r") as mothur:
                    line = mothur.readline()
                    while line != "":
                         mothur_d[line.split("|")[1]] = line
                         line = mothur.readline()
                buffer = ""
                with open(F"{args.dir}/query_dbs/unite_partition_{k}_query_{r}{reg}__RDP_taxonomy.txt", "r") as rdp:
                    line = rdp.readline()
                    line = rdp.readline()
                    while line != "":
                        out_line = mothur_d[line.split('\t')[0]]
                        buffer = F"{buffer}{out_line}"
                        line = rdp.readline()
                with open(F"{args.dir}/sorted_unite_mothur_{c}_{k}_query_{r}{reg}.tax", "w") as ofile:
                    ofile.write(buffer)
                with open(F"{args.dir}/sorted_unite_mothur_{c}_{k}_query_{r}{reg}.tax", "r") as ifile:
                    buffer = "OTU_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
                    line = ifile.readline()
                    while line != "":
                        tax_string = line.split("\t")[-1]
                        tax_list1 = tax_string.strip(";\n").split(";")
                        tax_list2 = [x[3:] for x in tax_list1]
                        tax_list3 = [x.split("(")[0] for x in tax_list2]
                        i = 0
                        tax_list4 = [""]*7
                        taxa = tax_list3[i]
                        while "unclassified" not in taxa:
                            tax_list4[i] = taxa
                            i += 1
                            if i < len(tax_list3):
                                taxa = tax_list3[i]
                            else:
                                break
                        tax_out = '\t'.join(tax_list4)
                        buffer = F"{buffer}{line.split('|')[1]}\t{tax_out}\n"
                        line = ifile.readline()
                with open(F"{args.dir}/sorted_unite_mothur_{c}_{k}_query_{r}{reg}_comb.txt", "w") as ofile:
                    ofile.write(buffer)
