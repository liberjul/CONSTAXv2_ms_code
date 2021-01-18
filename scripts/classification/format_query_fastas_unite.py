import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

for k in range(5):
    for r in ["gen", "fam"]:
        for reg in ["", "_itsx.ITS1", "_itsx.ITS2"]:
            print(k, r, reg)
            buffer = ""
            line_count = 0
            with open(F"{args.dir}/unite_partition_{k}_query_{r}{reg}.fasta", "r") as ifile, open(F"{args.dir}/unite_partition_{k}_query_{r}{reg}_qm_fmt.fasta", "w") as ofile:
                line = ifile.readline()
                while line != "":
                    fmt_line = line.upper().split("\\")
                    if len(fmt_line) == 1:
                        buffer = F"{buffer}{fmt_line[0]}"
                    else:
                        buffer = F"{buffer}{fmt_line[0]}\n"
                    line_count += 1
                    if line_count == 1000:
                        ofile.write(buffer)
                        buffer = ""
                        line_count = 0
                    line = ifile.readline()
                if line_count > 0:
                    ofile.write(buffer)
