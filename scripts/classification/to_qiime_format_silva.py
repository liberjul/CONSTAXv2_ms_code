import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input database file")
parser.add_argument("--otax", type=str)
parser.add_argument("--oalign", type=str)
args = parser.parse_args()
# Formatted SILVA 138 release files downloaded from https://github.com/qiime2/docs/blob/master/source/data-resources.rst

ID_dict = {}
align_buf = ""
tax_buf = ""
line_count = 0
with open(args.input, "r") as ifile, open(args.oalign, "w") as oalign, open(args.otax, "w") as otax:
    line = ifile.readline()
    while line != "":
        line_count += 1
        if line[0] == ">":
            spl_line = line.strip().replace(" Bacteria;", "|Bacteria;").split("|")
            ID, tax = spl_line
            align_buf = F"{align_buf}{ID}\n"
            tax_list = tax.replace(" ", "_").split(";")
            tax_list.extend(["unidentified"]*(7-len(tax_list)))
            tax_str = ""
            for i,x in enumerate("dpcofgs"):
                tax_str = F"{tax_str}{x}__{tax_list[i]};"
            tax_buf = F"{tax_buf}{ID[1:]}\t{tax_str[:-1]}\n"
        else:
            line = line.replace("U", "T").replace("u", "t")
            align_buf = F"{align_buf}{line}"
        if (line_count % 1000) == 0:
            print(line_count)
            oalign.write(align_buf)
            otax.write(tax_buf)
            align_buf = ""
            tax_buf = ""
        line = ifile.readline()
    oalign.write(align_buf)
    otax.write(tax_buf)
