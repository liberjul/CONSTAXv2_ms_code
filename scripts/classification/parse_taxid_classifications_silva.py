tax_from_id = {'0' : "\t\t\t\t\t\t"}
with open("./silva_test/data/tax_slv_ssu_138.txt", "r") as ifile:
    line = ifile.readline()
    while line != "":
        tax_string, taxid = line.split("\t")[:2]
        tax_list = tax_string.strip("uncultured;").split(";")
        if "Bacteria" == tax_list[0]:
            tax_tab_string = "\t".join(tax_list) + "\t"*(7-len(tax_list))
            tax_from_id[taxid] = tax_tab_string
        line = ifile.readline()

for k in range(5):
    for r in ["gen", "fam"]:
        for reg in ["v4", "v3-4", "full"]:
            if reg == "full":
                datafile = F"silva_kraken2_p{k}_query_sub_1k_{r}.txt"
            else:
                datafile = F"silva_kraken2_p{k}_query_{r}_{reg}.txt"
            with open(datafile, "r") as ifile:
                buffer = "Seq_ID\tRank_1\tRank_2\tRank_3\tRank_4\tRank_5\tRank_6\tRank_7\n"
                line = ifile.readline()
                while line != "":
                    spl = line.split("\t")
                    seqid = spl[1].split(".")[0]
                    class_id = spl[2].split("taxid ")[1].split(")")[0]
                    buffer = F"{buffer}{seqid}\t{tax_from_id[class_id]}\n"
                    line = ifile.readline()
            with open(F"sil_parts_{k}_{r}_kraken2_query_{reg}_taxonomy_comb.txt", "w") as ofile:
                ofile.write(buffer)
