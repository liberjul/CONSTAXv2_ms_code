tax_from_id = {'0' : "\t\t\t\t\t\t"}
with open("./unite_test/data/unite_tax.txt", "r") as ifile:
    line = ifile.readline()
    while line != "":
        tax_string, taxid = line.split("\t")[:2]
        tax_list = tax_string.strip(";").split(";")
        tax_tab_string = "\t".join(tax_list) + "\t"*(7-len(tax_list))
        tax_from_id[taxid] = tax_tab_string
        line = ifile.readline()

for k in range(5):
    for r in ["gen", "fam"]:
        for reg in ["its1", "its2", "full"]:
            if reg == "full":
                datafile = F"unite_kraken2_p{k}_query_{r}.txt"
            else:
                datafile = F"unite_kraken2_p{k}_query_{r}_{reg}.txt"
            with open(datafile, "r") as ifile:
                buffer = "Seq_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
                line = ifile.readline()
                while line != "":
                    spl = line.split("\t")
                    seqid = spl[1].split("|")[1]
                    class_id = spl[2].split("taxid ")[1].split(")")[0]
                    buffer = F"{buffer}{seqid}\t{tax_from_id[class_id]}\n"
                    line = ifile.readline()
            with open(F"uni_parts_{k}_{r}_kraken2_query_{reg}_taxonomy_comb.txt", "w") as ofile:
                ofile.write(buffer)
