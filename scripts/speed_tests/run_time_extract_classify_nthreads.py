import sys

with open(sys.argv[1],"r") as ifile:
    buf = "Iteration,step,threads,algorithm,seq_count,time\n"
    iter = 0
    line = ifile.readline()
    while line != "":
        if "Started at" in line:
            start, threads, algorithm, q = line.strip().split("|")
            start = float(start.split("at ")[1])
        elif "Training completed at" in line:
            training = float(line.strip().split("at ")[1])
            buf = F"{buf}{iter},train,{threads},{algorithm},{q},{training - start}\n"
        elif "Finished at" in line:
            finish = float(line.strip().split("at ")[1])
            buf = F"{buf}{iter},classify,{threads},{algorithm},{q},{finish - training}\n"
            iter += 1
        line = ifile.readline()
with open(sys.argv[2], "w") as ofile:
    ofile.write(buf)
