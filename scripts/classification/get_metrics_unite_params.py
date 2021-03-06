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

    return [sens, mc, oc, epq]

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", type=str, help="data directory")
args = parser.parse_args()

level_dict = {"fam" : 5, "gen" : 6}
buffer = "database,partition_level,k_iteration,conf,max_hits,classifier,sensitivity,MC,OC,EPQ\n"
for i in ["fam", "gen"]:
    for j in range(5):
        for d in ["unite"]: #, "silva"]:
            for conf in [0.2, 0.4, 0.6, 0.8, 1.0]:
                for mhits in [1, 3, 5, 10, 20]:
                    print(F"{args.dir}/out_p{j}_{i}")
                    try:
                        print(F"{args.dir}/out_param/combined_taxonomy_{i}_p{j}_c{conf}_m{mhits}.txt")
                        comb_tax = pd.read_table(F"{args.dir}/out_param/combined_taxonomy_{i}_p{j}_c{conf}_m{mhits}.txt")
                        rdp = comb_tax.iloc[:,1::4]
                        blast = comb_tax.iloc[:,2::4]
                        sintax = comb_tax.iloc[:,3::4]
                        cons = comb_tax.iloc[:,4::4]
                        rdp, blast, sintax, cons = [x.fillna("-") for x in [rdp, blast, sintax, cons]]

                        known_df = pd.read_table(F"{args.dir}/query_dbs/unite_partition_{j}_query_sub_1k_{i}__RDP_taxonomy.txt")
                        for pred_df in [rdp, blast, sintax, cons]:
                            metrics = score_df(pred_df, known_df, level_known=level_dict[i])
                            metrics = [str(x) for x in metrics]
                            classifier = pred_df.columns[0].split("_")[1]
                            buffer = F"{buffer}{d},{i},{j},{conf},{mhits},{classifier},{','.join(metrics)}\n"
                    except:
                        pass
with open(F"{args.dir}/part_cv_metrics_unite_params.csv", "w") as ofile:
    ofile.write(buffer)
