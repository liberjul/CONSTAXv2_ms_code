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
    n_known = len(flat_k_tax)
    n_novel = len(flat_pou_tax)
    n = n_known + n_novel

    return [sens, mc, oc, epq, n_known, n_novel, n]

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

level_dict = {"fam" : 5, "gen" : 6}
buffer = "database,partition_level,region,k_iteration,classifier,sensitivity,MC,OC,EPQ,N_known,N_novel,N\n"
region_dict = {"Full" : "_full" , "V3-4" : "_v3-4", "V4" : "_v4"}
for k in range(5):
    for r in ["gen", "fam"]:
        for region in region_dict.keys():
            print(k, r, region)
            # comb_tax = pd.read_table(F"uni_parts_{k}_{r}_qiime_query{region_dict[region][1]}_taxonomy_comb.txt")
            # cons = comb_tax.iloc[:,1:]
            # cons = cons.fillna("-")

            name_dict = {}
            if region == "Full":
                filename_known = F"{args.dir}/query_dbs/silva_partition_{k}_query_sub_1k_{r}__RDP_taxonomy.txt"
            else:
                filename_known = F"{args.dir}/query_dbs/silva_partition_{k}_query_{r}{region_dict[region]}__RDP_taxonomy.txt"
            with open(filename_known, "r") as ifile:
                line = ifile.readline()
                line = ifile.readline()
                while line != "":
                    name = line.split("\t")[0]
                    name_dict[name] = None
                    line = ifile.readline()
            with open(F"sil_parts_{k}_{r}_qiime_query{region_dict[region]}_taxonomy_comb.txt", "r") as ifile:
                with open("filtered_tax.txt", "w") as ofile:
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
            comb_tax = pd.read_table("filtered_tax.txt")
            cons = comb_tax.iloc[:,1:]
            cons = cons.fillna("-")

            # known_df = pd.read_table("filtered_tax.txt")
            known_df = pd.read_table(filename_known)
            metrics = score_df(cons, known_df, level_known=level_dict[r])
            metrics = [str(x) for x in metrics]
            buffer = F"{buffer}silva,{r},{region},{k},qiime2-Naive-Bayes,{','.join(metrics)}\n"
with open("part_cv_metrics_silva_qiime.csv", "w") as ofile:
    ofile.write(buffer)
