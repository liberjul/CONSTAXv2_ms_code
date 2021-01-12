import argparse
# Sort mothur output taxonomy file by RDP taxonomy file
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

for k in range(5):
    for r in ["gen", "fam"]:
        for c in ["wang", "knn"]:
            for reg in ["_v3-4", "_v4"]:
                print(k, r, c, reg)
                mothur_d = {}
                with open(F"{args.dir}/silva_partition_{k}_query_{r}{reg}.silva_partition_{k}_ref_{r}.{c}.taxonomy", "r") as mothur:
                    line = mothur.readline()
                    while line != "":
                         mothur_d[line.split("\t")[0].split(".")[0]] = line
                         line = mothur.readline()
                buffer = ""
                with open(F"{args.dir}/query_dbs/silva_partition_{k}_query_{r}{reg}__RDP_taxonomy.txt", "r") as rdp:
                    header = rdp.readline()
                    line = rdp.readline()
                    while line != "":
                        out_line = mothur_d[line.split('\t')[0]]
                        buffer = F"{buffer}{out_line}"
                        line = rdp.readline()
                with open(F"{args.dir}/sorted_silva_mothur_{c}_{k}_query_{r}{reg}.tax", "w") as ofile:
                    ofile.write(buffer)
                with open(F"{args.dir}/sorted_silva_mothur_{c}_{k}_query_{r}{reg}.tax", "r") as ifile:
                    buffer = header.replace("Seq_ID", "OTU_ID")
                    rank_count = len(header.split("\t")) - 1
                    line = ifile.readline()
                    while line != "":
                        ID = line.split('\t')[0].split(".")[0]
                        tax_string = line.split("\t")[-1]
                        tax_list1 = tax_string.strip(";\n").split(";")
                        tax_list2 = [x for x in tax_list1]
                        tax_list3 = [x.split("(")[0] for x in tax_list2]
                        i = 0
                        tax_list4 = [""]*rank_count
                        taxa = tax_list3[i]
                        while "unclassified" not in taxa and i < rank_count:
                            tax_list4[i] = taxa
                            i += 1
                            if i < len(tax_list3):
                                taxa = tax_list3[i]
                            else:
                                break
                        tax_out = '\t'.join(tax_list4)
                        buffer = F"{buffer}{ID}\t{tax_out}\n"
                        line = ifile.readline()
                with open(F"{args.dir}/sorted_silva_mothur_{c}_{k}_query_{r}{reg}_comb.txt", "w") as ofile:
                    ofile.write(buffer)
            print(k, r, c, "full")
            mothur_d = {}
            with open(F"{args.dir}/silva_partition_{k}_query_sub_1k_{r}.silva_partition_{k}_ref_{r}.{c}.taxonomy", "r") as mothur:
                line = mothur.readline()
                while line != "":
                     mothur_d[line.split("\t")[0].split(".")[0]] = line
                     line = mothur.readline()
            buffer = ""
            with open(F"{args.dir}/query_dbs/silva_partition_{k}_query_sub_1k_{r}__RDP_taxonomy.txt", "r") as rdp:
                header = rdp.readline()
                line = rdp.readline()
                while line != "":
                    out_line = mothur_d[line.split('\t')[0]]
                    buffer = F"{buffer}{out_line}"
                    line = rdp.readline()
            reg = "_full"
            with open(F"{args.dir}/sorted_silva_mothur_{c}_{k}_query_{r}{reg}.tax", "w") as ofile:
                ofile.write(buffer)
            with open(F"{args.dir}/sorted_silva_mothur_{c}_{k}_query_{r}{reg}.tax", "r") as ifile:
                buffer = header.replace("Seq_ID", "OTU_ID")
                rank_count = len(header.split("\t")) - 1
                line = ifile.readline()
                while line != "":
                    ID = line.split('\t')[0].split(".")[0]
                    tax_string = line.split("\t")[-1]
                    tax_list1 = tax_string.strip(";\n").split(";")
                    tax_list2 = [x for x in tax_list1]
                    tax_list3 = [x.split("(")[0] for x in tax_list2]
                    i = 0
                    tax_list4 = [""]*rank_count
                    taxa = tax_list3[i]
                    while "unclassified" not in taxa and i < rank_count:
                        tax_list4[i] = taxa
                        i += 1
                        if i < len(tax_list3):
                            taxa = tax_list3[i]
                        else:
                            break
                    tax_out = '\t'.join(tax_list4)
                    buffer = F"{buffer}{ID}\t{tax_out}\n"
                    line = ifile.readline()
            with open(F"{args.dir}/sorted_silva_mothur_{c}_{k}_query_{r}{reg}_comb.txt", "w") as ofile:
                ofile.write(buffer)
