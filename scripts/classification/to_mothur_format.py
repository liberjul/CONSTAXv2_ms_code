import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input database file")
parser.add_argument("--itax", type=str)
parser.add_argument("--otax", type=str)
parser.add_argument("--ialign", type=str)
parser.add_argument("--oalign", type=str)
args = parser.parse_args()
# >Entoloma_vindobonense|JX454802|SH1569086.08FU|refs|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Entolomataceae;g__Entoloma;s__Entoloma_vindobonense

ID_dict = {}
with open(args.input, "r") as ifile:
    line = ifile.readline()
    while line != "":
        spl_line = line.split("|")
        if len(spl_line) > 3:
            ID = F"{spl_line[2]}_{spl_line[1]}_{spl_line[3]}"
            ID_dict[ID] = None
        line = ifile.readline()
with open(args.itax, "r") as itax, open(args.otax, "w") as otax:
    line = itax.readline()
    while line != "":
        spl_line = line.split("\t")
        if spl_line[0] in ID_dict:
            otax.write(line)
        line = itax.readline()
with open(args.ialign, "r") as ialign, open(args.oalign, "w") as oalign:
        line = ialign.readline()
        while line != "":
            if line[1:].strip() in ID_dict:
                oalign.write(line)
                line = ialign.readline()
                oalign.write(line)
                line = ialign.readline()
            else:
                line = ialign.readline()
                line = ialign.readline()
