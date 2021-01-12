import sys

with open(sys.argv[1],"r") as ifile:
    buf = "Iteration,step,algorithm,seq_count,time\n"
    iter = 0
    line = ifile.readline()
    while line != "":
        if "Started at" in line:
            start, tail = line.strip().split(" Iteration ")
            iter, tail = tail.split(" ")
            algorithm, seq_count = tail.split("|")
            start = float(start.split("at ")[1])
        elif "Training completed at" in line:
            training = float(line.strip().split("at ")[1])
            buf = F"{buf}{iter},train,{algorithm},{seq_count},{training - start}\n"
        elif "Finished at" in line:
            finish = float(line.strip().split("at ")[1])
            buf = F"{buf}{iter},classify,{algorithm},{seq_count},{finish - training}\n"
        line = ifile.readline()
with open(sys.argv[2], "w") as ofile:
    ofile.write(buf)
