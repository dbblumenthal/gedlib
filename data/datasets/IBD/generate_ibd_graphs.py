#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2020 by David B. Blumenthal                              #
#                                                                          #
#   This file is part of GEDLIB.                                           #
#                                                                          #
#   GEDLIB is free software: you can redistribute it and/or modify it      #
#   under the terms of the GNU Lesser General Public License as published  #
#   by the Free Software Foundation, either version 3 of the License, or   #
#   (at your option) any later version.                                    #
#                                                                          #
#   GEDLIB is distributed in the hope that it will be useful,              #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           #
#   GNU Lesser General Public License for more details.                    #
#                                                                          #
#   You should have received a copy of the GNU Lesser General Public       #
#   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. #
#                                                                          #
#//////////////////////////////////////////////////////////////////////////#

import pandas as pd
import numpy as np
import argparse
import os.path

def write_gxl_files(samples, data, nlogratios, otus, means, stds, contained_in_num_samples, cutoff, min_samples):
    for sample in samples:
        gxl_file_name = os.path.join("data", "{}.gxl".format(sample))
        gxl_file = open(gxl_file_name, "w")
        gxl_file.write("<?xml version=\"1.0\"?>\n")
        gxl_file.write("<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n")
        gxl_file.write("<gxl>\n")
        gxl_file.write("<graph id=\"{}.gxl\" edgeids=\"false\" edgemode=\"undirected\">\n".format(sample))
        num_nodes = 0
        node_ids_to_otus = {}
        for otu in otus:
            if data.loc[sample, otu] > 0:
                node_ids_to_otus[num_nodes] = otu
                gxl_file.write("\t<node id=\"_{}\">\n".format(num_nodes))
                gxl_file.write("\t\t<attr name=\"OTU\"><int>{}</int></attr>\n".format(otu.split('_')[1]))
                gxl_file.write("\t</node>\n")
                num_nodes += 1
        for node_id_1 in range(num_nodes):
            otu_1 = node_ids_to_otus[node_id_1]
            for node_id_2 in range(node_id_1+1, num_nodes):
                otu_2 = node_ids_to_otus[node_id_2]
                nlogratio = nlogratios[sample].loc[otu_1, otu_2]
                if pd.isna(nlogratio):
                    raise Exception("normalized log-ratio is NA")
                zscore = cutoff + 1
                if contained_in_num_samples.loc[otu_1, otu_2] >= min_samples:
                    zscore = (nlogratio - means.loc[otu_1, otu_2]) / stds.loc[otu_1,otu_2]
                if np.abs(zscore) > cutoff:
                    gxl_file.write("\t<edge from=\"_{}\" to=\"_{}\">\n".format(node_id_1, node_id_2))
                    gxl_file.write("\t\t<attr name=\"nlogratio\"><float>{}</float></attr>\n".format(nlogratio))
                    gxl_file.write("\t</edge>\n")
        gxl_file.write("</graph>\n")
        gxl_file.write("</gxl>\n")
        gxl_file.close()

def compute_log_ratios():
    print("Computing log-ratios.")
    data = pd.read_csv("otu_abundances.csv", index_col=0)
    samples = data.index
    otus = data.columns
    nlogratios = {sample : pd.DataFrame(index=otus, columns=otus, dtype=float) for sample in samples}
    for sample in samples:
        for otu_id_1 in range(len(otus)):
            otu_1 = otus[otu_id_1]
            if data.loc[sample,otu_1] > 0:
                for otu_id_2 in range(otu_id_1+1, len(otus)):
                    otu_2 = otus[otu_id_2]
                    if data.loc[sample,otu_2] > 0:
                        nlogratios[sample].loc[otu_1,otu_2] = np.log(data.loc[sample,otu_1]/data.loc[sample,otu_2])
    return samples, otus, data, nlogratios

def normalize_log_ratios(samples, otus, nlogratios):
    print("Normalizing log-ratios.")                    
    global_max = 0
    global_min = 0
    for sample in samples:
        values = [value for value in nlogratios[sample].values.ravel().tolist() if not pd.isna(value)]
        global_max = max(global_max, np.max(values))
        global_min = min(global_min, np.min(values))
    for sample in samples:
        for otu_id_1 in range(len(otus)):
            otu_1 = otus[otu_id_1]
            for otu_id_2 in range(otu_id_1+1, len(otus)):
                otu_2 = otus[otu_id_2]
                if not pd.isna(nlogratios[sample].loc[otu_1,otu_2]):
                    nlogratios[sample].loc[otu_1,otu_2] = (nlogratios[sample].loc[otu_1,otu_2] - global_min) / (global_max - global_min)
    return nlogratios

def aggregate(samples, otus, nlogratios):
    print("Computing means and standard deviations.")               
    means = pd.DataFrame(index=otus, columns=otus, dtype=float)
    stds = pd.DataFrame(index=otus, columns=otus, dtype=float)
    contained_in_num_samples = pd.DataFrame(index=otus, columns=otus, dtype=int)
    for otu_id_1 in range(len(otus)):
        otu_1 = otus[otu_id_1]
        for otu_id_2 in range(otu_id_1+1, len(otus)):
            otu_2 = otus[otu_id_2]
            ratios = [nlogratios[sample].loc[otu_1, otu_2] for sample in samples if not pd.isna(nlogratios[sample].loc[otu_1, otu_2])]
            contained_in_num_samples.loc[otu_1,otu_2] = len(ratios) 
            if len(ratios) > 0:
                means.loc[otu_1,otu_2] = np.mean(ratios)    
                stds.loc[otu_1,otu_2] = np.std(ratios)
    return means, stds, contained_in_num_samples

def run_script():
    parser = argparse.ArgumentParser(description="Generate IBD graphs.")
    parser.add_argument("--cutoff", help="only edges with absolute Z-scores greater than cutoff are included in the graphs", type=float, default=2.0)
    parser.add_argument("--min-samples", help="apply cutoff only if the edge is present in at least min-samples samples", type=int, default=5)
    args = parser.parse_args()
    samples, otus, data, nlogratios = compute_log_ratios()
    nlogratios = normalize_log_ratios(samples, otus, nlogratios)
    means, stds, contained_in_num_samples = aggregate(samples, otus, nlogratios)
    write_gxl_files(samples, data, nlogratios, otus, means, stds, contained_in_num_samples, args.cutoff, args.min_samples)
    write_gxl(samples, data, nlogratios, otus, means, stds, contained_in_num_samples, args.cutoff, args.min_samples)

if __name__ == "__main__":
    run_script()