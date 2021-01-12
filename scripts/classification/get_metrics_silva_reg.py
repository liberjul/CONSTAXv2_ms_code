import numpy as np
import pandas as pd
import argparse

def sensitivity(predict, target):
    '''
    True positives / known labels count
    '''
    tp = np.sum(predict == target)
    print("True positives: ", tp)
    count = len(target)
    print("Known labels: ", count)
    return tp / count

def misclass(predict, target):
    '''
    Misclassified false positives (known but wrong) / known labels count
    '''
    fp_mis = np.sum(predict != target)
    print("False positive MC: ", fp_mis)
    count = len(target)
    return fp_mis / count
def overclass(predict):
    '''
    Overclassified false positive (unknown but assigned) / novel label count
    '''
    fp_over = np.sum(predict != "-")
    return fp_over/len(predict)

def errors_per_query(pok, pou, target):
    '''
    (MC false positives + OC false positives) / (known labels + novel labels)
    '''
    fp_mis = np.sum(pok != target)
    fp_over = np.sum(pou != "-")
    return (fp_mis + fp_over)/(len(target) + len(pou))

def score_df(pred_df, known_df, level_known):
    known_tax = known_df.iloc[:, 1:level_known + 1]
    flat_k_tax = known_tax.to_numpy().flatten()

    pred_tax = pred_df.iloc[:, :level_known]
    flat_pok_tax = pred_tax.to_numpy(dtype="str").flatten()
    flat_pok_tax[np.core.defchararray.find(flat_pok_tax,"ncertae_sedis")!=-1] = "-"
    flat_pok_tax = [x.split("_")[0] for x in flat_pok_tax]
    flat_pok_tax = np.array(flat_pok_tax)

    print(flat_k_tax)
    print("Shape known: ", known_tax.shape)
    print(flat_pok_tax)
    print("Shape predicted: ", pred_tax.shape)
    unk_pred_tax = pred_df.iloc[:, level_known:]
    flat_pou_tax = unk_pred_tax.to_numpy(dtype="str").flatten()
    flat_pou_tax[np.core.defchararray.find(flat_pou_tax,"ncertae_sedis")!=-1] = "-"

    sens = sensitivity(flat_pok_tax, flat_k_tax)
    mc = misclass(flat_pok_tax, flat_k_tax)
    oc = overclass(flat_pou_tax)
    epq = errors_per_query(flat_pok_tax, flat_pou_tax, flat_k_tax)

    return [sens, mc, oc, epq]

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

level_dict = {"fam" : 5, "gen" : 6}
region_dict = {"Full" : "", "V4" : "_v4", "V3-4" : "_v3-4"}
buffer = "database,partition_level,region,k_iteration,classifier,sensitivity,MC,OC,EPQ\n"
for i in ["fam", "gen"]:
    for j in range(5):
        for d in ["silva"]: #, "silva"]:
            for region in region_dict.keys():
                print(F"{args.dir}/out_p{j}_{i}")
                if region == "Full":
                    pred_file = F"{args.dir}/out_p{j}_{i}_1k_cons/combined_taxonomy.txt"
                else:
                    pred_file = F"{args.dir}/out_p{j}_{i}{region_dict[region]}_cons/combined_taxonomy.txt"
                comb_tax = pd.read_table(pred_file)
                rdp = comb_tax.iloc[:,1::4]
                blast = comb_tax.iloc[:,2::4]
                sintax = comb_tax.iloc[:,3::4]
                cons = comb_tax.iloc[:,4::4]
                rdp, blast, sintax, cons = [x.fillna("-") for x in [rdp, blast, sintax, cons]]

                name_dict = {}

                with open(pred_file, "r") as ifile:
                    line = ifile.readline()
                    line = ifile.readline()
                    while line != "":
                        name = line.split(".")[0]
                        name_dict[name] = None
                        line = ifile.readline()
                if region == "Full":
                    tax_file = F"{args.dir}/query_dbs/silva_partition_{j}_query_sub_1k_{i}__RDP_taxonomy.txt"
                else:
                    tax_file = F"{args.dir}/query_dbs/silva_partition_{j}_query_{i}{region_dict[region]}__RDP_taxonomy.txt"
                with open(tax_file, "r") as ifile:
                    with open(F"{args.dir}/filtered_tax.txt", "w") as ofile:
                        line = ifile.readline()
                        ofile.write(line)
                        line = ifile.readline()
                        while line != "":
                            name = line.split("\t")[0]
                            if name not in name_dict:
                                print(name)
                            else:
                                ofile.write(line)
                            line = ifile.readline()

                known_df = pd.read_table(F"{args.dir}/filtered_tax.txt")
                for pred_df in [rdp, blast, sintax, cons]:
                    metrics = score_df(pred_df, known_df, level_known=level_dict[i])
                    metrics = [str(x) for x in metrics]
                    classifier = pred_df.columns[0].split("_")[-1]
                    buffer = F"{buffer}{d},{i},{region},{j},{classifier},{','.join(metrics)}\n"
with open(F"{args.dir}/part_cv_metrics_silva_reg_conservative.csv", "w") as ofile:
    ofile.write(buffer)
buffer = "database,partition_level,region,k_iteration,classifier,sensitivity,MC,OC,EPQ\n"
for i in ["fam", "gen"]:
    for j in range(5):
        for d in ["silva"]: #, "silva"]:
            for region in region_dict.keys():
                print(F"{args.dir}/out_p{j}_{i}")
                if region == "Full":
                    pred_file = F"{args.dir}/out_p{j}_{i}_1k/combined_taxonomy.txt"
                else:
                    pred_file = F"{args.dir}/out_p{j}_{i}{region_dict[region]}/combined_taxonomy.txt"
                comb_tax = pd.read_table(pred_file)
                rdp = comb_tax.iloc[:,1::4]
                blast = comb_tax.iloc[:,2::4]
                sintax = comb_tax.iloc[:,3::4]
                cons = comb_tax.iloc[:,4::4]
                rdp, blast, sintax, cons = [x.fillna("-") for x in [rdp, blast, sintax, cons]]

                name_dict = {}

                with open(pred_file, "r") as ifile:
                    line = ifile.readline()
                    line = ifile.readline()
                    while line != "":
                        name = line.split(".")[0]
                        name_dict[name] = None
                        line = ifile.readline()
                if region == "Full":
                    tax_file = F"{args.dir}/query_dbs/silva_partition_{j}_query_sub_1k_{i}__RDP_taxonomy.txt"
                else:
                    tax_file = F"{args.dir}/query_dbs/silva_partition_{j}_query_{i}{region_dict[region]}__RDP_taxonomy.txt"
                with open(tax_file, "r") as ifile:
                    with open(F"{args.dir}/filtered_tax.txt", "w") as ofile:
                        line = ifile.readline()
                        ofile.write(line)
                        line = ifile.readline()
                        while line != "":
                            name = line.split("\t")[0]
                            if name not in name_dict:
                                print(name)
                            else:
                                ofile.write(line)
                            line = ifile.readline()

                known_df = pd.read_table(F"{args.dir}/filtered_tax.txt")
                for pred_df in [rdp, blast, sintax, cons]:
                    metrics = score_df(pred_df, known_df, level_known=level_dict[i])
                    metrics = [str(x) for x in metrics]
                    classifier = pred_df.columns[0].split("_")[-1]
                    buffer = F"{buffer}{d},{i},{region},{j},{classifier},{','.join(metrics)}\n"
with open(F"{args.dir}/part_cv_metrics_silva_reg.csv", "w") as ofile:
    ofile.write(buffer)
