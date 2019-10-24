#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2019 by David B. Blumenthal                              #
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

##
# @file process_results.py
# @brief Processes the results of median tests.
#
# @details 
# Usage: 
# ```sh
# $ python process_results.py <RESULT_FILE> <OUTPUT_PREFIX>
# ```
#
# For more information, execute `$ python process_results.py --help`.

'''Processes the results of median tests.'''

import pandas as pd
import numpy as np
import argparse

def aggregate_results(csv_file):
    rs = pd.read_csv(csv_file)
    rs.itrs = rs.itrs.apply(lambda l: [int(x) for x in l.split(";")])
    grouped = rs.groupby(["percent", "init_type", "num_inits", "algo"], as_index=False)
    concat_lists = lambda ll: [x for l in ll for x in l if x > 0] 
    as_list = lambda ll: [x for x in ll]
    ratio_greater = lambda l, t: float(len([x for x in l if x > t])) / float(len(l))
    agg_rs = grouped.agg({"time": np.mean, "time_init": np.mean, "time_converged": np.mean, "sod": np.mean, "sod_init": np.mean, "sod_converged": np.mean, "itrs": concat_lists, "state": as_list})
    agg_rs["itrs_max"] = agg_rs.itrs.apply(lambda l: max(l) if len(l) > 0 else 0)
    agg_rs["itrs_min"] = agg_rs.itrs.apply(lambda l: min(l) if len(l) > 0 else None)
    agg_rs["itrs_mean"] = agg_rs.itrs.apply(lambda l: np.mean(l) if len(l) > 0 else None)
    agg_rs["itrs_lower_quartile"] = agg_rs.itrs.apply(lambda l: np.quantile(l, .25) if len(l) > 0 else None)
    agg_rs["itrs_median"] = agg_rs.itrs.apply(lambda l: np.quantile(l, .5)if len(l) > 0 else None)
    agg_rs["itrs_upper_quartile"] = agg_rs.itrs.apply(lambda l: np.quantile(l, .75)if len(l) > 0 else None)
    agg_rs["ratio_initialized"] = agg_rs.state.apply(lambda l: ratio_greater(l, 0) if len(l) > 0 else None)
    agg_rs["ratio_converged"] = agg_rs.state.apply(lambda l: ratio_greater(l, 1) if len(l) > 0 else None)
    agg_rs["ratio_refined"] = agg_rs.state.apply(lambda l: ratio_greater(l, 2) if len(l) > 0 else None)
    return agg_rs

def compute_rs_vs_percent(agg_rs):
    row_rs = pd.DataFrame([10 * i for i in range(1, 11)], columns = ["percent"])
    algos = ["BRANCH_FAST", "REFINE", "IPFP"]
    init_types = ["MAX", "MIN", "MEAN", "MEDOID", "RANDOM"]
    nums_inits = [2 ** i for i in range(6)]
    cols = ['time', 'time_init', 'time_converged', 'sod', 'sod_init', 'sod_converged', 'itrs', 'state', 'itrs_max', 'itrs_min', 'itrs_mean', 'itrs_lower_quartile', 'itrs_median', 'itrs_upper_quartile', 'ratio_initialized', 'ratio_converged', 'ratio_refined']
    for algo in algos:
        for init_type in init_types:
            for num_inits in nums_inits:
                if num_inits > 1 and init_type != "RANDOM":
                    continue
                prefix = algo + "_" + init_type + "_I" + str(num_inits) + "_"
                for col in cols:
                    row_rs[prefix + col] = agg_rs[(agg_rs.algo == algo) & (agg_rs.init_type == init_type) & (agg_rs.num_inits == num_inits)][col].reset_index()[col]
    return row_rs

def compute_rs_vs_num_inits(agg_rs):
    row_rs = pd.DataFrame([2 ** i for i in range(6)], columns = ["num_inits"])
    algos = ["BRANCH_FAST", "REFINE", "IPFP"]
    percents = [10 * i for i in range(1, 11)]
    cols = ['percent', 'time', 'time_init', 'time_converged', 'sod', 'sod_init', 'sod_converged', 'itrs', 'state', 'itrs_max', 'itrs_min', 'itrs_mean', 'itrs_lower_quartile', 'itrs_median', 'itrs_upper_quartile', 'ratio_initialized', 'ratio_converged', 'ratio_refined']
    for algo in algos:
        for percent in percents:
            prefix = algo + "_RANDOM_P" + str(percent) + "_"
            for col in cols:
                row_rs[prefix + col] = agg_rs[(agg_rs.algo == algo) & (agg_rs.init_type == "RANDOM") & (agg_rs.percent == percent)][col].reset_index()[col]
    return row_rs

def create_init_plots(rs_vs_num_inits_letter, rs_vs_num_inits_aids, rs_vs_num_inits_muta):
    tikz = open("../img/effect_num_inits.tex", "w")
    tikz.write("%!TEX root = ../main.tex\n")
    tikz.write("\\begin{tikzpicture}\n")
    tikz.write("\\begin{groupplot}[group style={group name=groupplot, group size=3 by 2, horizontal sep=1.1cm, vertical sep=1.5cm},\n")
    tikz.write("width=.36\linewidth,\n")
    tikz.write("height=.25\linewidth,\n")
    tikz.write("legend columns=3,\n")
    tikz.write("legend cell align=left,\n")
    tikz.write("legend style={align=left, draw=none, column sep=.5ex},\n")
    tikz.write("/tikz/font=\small]\n")
    create_lineplots_inits_vs_sod(rs_vs_num_inits_letter, tikz, "\\letter (\\SI{90}{\\percent})", True)
    create_lineplots_inits_vs_sod(rs_vs_num_inits_aids, tikz, "\\aids (\\SI{90}{\\percent})")
    create_lineplots_inits_vs_sod(rs_vs_num_inits_muta, tikz, "\\muta (\\SI{90}{\\percent})")
    create_barplots_inits_vs_state(rs_vs_num_inits_letter, tikz, "\\letter (\\SI{90}{\\percent})")
    create_barplots_inits_vs_state(rs_vs_num_inits_aids, tikz, "\\aids (\\SI{90}{\\percent})")
    create_barplots_inits_vs_state(rs_vs_num_inits_muta, tikz, "\\muta (\\SI{90}{\\percent})")
    tikz.write("\end{groupplot}\n")
    tikz.write("\\node at ($(groupplot c1r1.north) !.5! (groupplot c3r1.north)$) [inner sep=0pt,anchor=south, yshift=5ex] {\\pgfplotslegendfromname{grouplegend}};\n")
    tikz.write("\end{tikzpicture}")
    tikz.close()
    
def create_config_plots(rs_vs_percent_letter, rs_vs_percent_aids, rs_vs_percent_muta):
    tikz = open("../img/effect_configs.tex", "w")
    tikz.write("%!TEX root = ../main.tex\n")
    tikz.write("\\begin{tikzpicture}\n")
    tikz.write("\\begin{groupplot}[group style={group name=groupplot, group size=3 by 2, horizontal sep=1.1cm, vertical sep=1.5cm},\n")
    tikz.write("width=.36\linewidth,\n")
    tikz.write("height=.25\linewidth,\n")
    tikz.write("legend columns=3,\n")
    tikz.write("legend cell align=left,\n")
    tikz.write("legend style={align=left, draw=none, column sep=.5ex},\n")
    tikz.write("/tikz/font=\small]\n")
    create_paretoplots_sod_vs_time(rs_vs_percent_letter, tikz, "\\letter (\\SI{90}{\\percent})", True)
    create_paretoplots_sod_vs_time(rs_vs_percent_aids, tikz, "\\aids (\\SI{90}{\\percent})")
    create_paretoplots_sod_vs_time(rs_vs_percent_muta, tikz, "\\muta (\\SI{90}{\\percent})")
    create_lineplots_max_itrs_vs_percent(rs_vs_percent_letter, tikz, "\\letter")
    create_lineplots_max_itrs_vs_percent(rs_vs_percent_aids, tikz, "\\aids")
    create_lineplots_max_itrs_vs_percent(rs_vs_percent_muta, tikz, "\\muta")
    tikz.write("\end{groupplot}\n")
    tikz.write("\\node at ($(groupplot c1r1.north) !.5! (groupplot c3r1.north)$) [inner sep=0pt,anchor=south, yshift=5ex] {\\pgfplotslegendfromname{grouplegend}};\n")
    tikz.write("\end{tikzpicture}")
    tikz.close()

def create_boxplots_inits_vs_itrs(rs_vs_num_inits, tikz, title):
    tikz.write("\\nextgroupplot[\n")
    tikz.write("title=" + title + ",\n")
    tikz.write("boxplot/draw direction=y,\n")
    tikz.write("xtick={2,6,10,14,18,22},\n")
    tikz.write("xticklabels={1,2,4,8,16,32},\n")
    tikz.write("xlabel={initial solutions},\n")
    tikz.write("ylabel={iterations},\n")
    tikz.write("ylabel style={yshift=-.1cm},\n")
    tikz.write("xlabel style={yshift=.1cm},\n")
    tikz.write("cycle list={},\n")
    tikz.write("]\n")
    xshift = {"BRANCH_FAST" : -1.1, "REFINE" : 0, "IPFP" : 1.1}
    color = {"BRANCH_FAST" : "BRANCHFAST", "REFINE" : "REFINE", "IPFP" : "IPFP"}
    index = {2 ** i : i for i in range(6)}
    pos = {2 ** i : 4 * i + 2 for i in range(6)}
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        col_prefix = algo + "_RANDOM_P90_itrs_";
        for num_inits in [2 ** i for i in range(6)]:
            tikz.write("\\addplot+[boxplot prepared={\n")
            tikz.write("draw position=" + str(pos[num_inits] + xshift[algo]) + ",\n")
            tikz.write("lower whisker=" + str(rs_vs_num_inits.at[index[num_inits], col_prefix + "min"]) + ",\n")
            tikz.write("lower quartile=" + str(rs_vs_num_inits.at[index[num_inits], col_prefix + "lower_quartile"]) + ",\n")
            tikz.write("median=" + str(rs_vs_num_inits.at[index[num_inits], col_prefix + "median"]) + ",\n")
            tikz.write("upper quartile=" + str(rs_vs_num_inits.at[index[num_inits], col_prefix + "upper_quartile"]) + ",\n")
            tikz.write("upper whisker=" + str(rs_vs_num_inits.at[index[num_inits], col_prefix + "max"]) + "\n")
            tikz.write("},\n")
            tikz.write("fill=" + color[algo] + "!65,\n")
            tikz.write("draw=" + color[algo] + ",\n")
            tikz.write("] coordinates {};\n")
    
    
def create_paretoplots_sod_vs_time(rs_vs_percent, tikz, title, legend = False):
    index = 8
    marks = {"MAX" : "mark=o", "MIN" : "mark=square", "MEAN" : "mark=triangle", "MEDOID" : "mark=diamond", "RANDOM" : "mark=pentagon", "MAXconv" : "mark=+", "MINconv" : "mark=x", "MEANconv" : "mark=-", "MEDOIDconv" : "mark=|", "RANDOMconv" : "mark=star"}
    colors = {"BRANCH_FAST" : "draw=BRANCHFAST", "REFINE" : "draw=REFINE", "IPFP" : "draw=IPFP"} 
    num_inits = {"MAX" : 1, "MIN" : 1, "MEAN" : 1, "MEDOID" : 1, "RANDOMBRANCH_FAST" : 8, "RANDOMREFINE" : 2, "RANDOMIPFP" : 2}
    macros = {"BRANCH_FAST" : "BRANCHFAST", "REFINE" : "REFINE", "IPFP" : "IPFP"}
    tikz.write("\\nextgroupplot[\n")
    tikz.write("title=" + title + ",\n")
    tikz.write("xlabel={mean runtime in \\si{\\second}},\n")
    tikz.write("ylabel={mean \SOD at termination},\n")
    tikz.write("ylabel style={yshift=-.1cm},\n")
    tikz.write("xlabel style={yshift=.1cm},\n")
    tikz.write("xmin=0.2,\n")
    tikz.write("xmax=1500,\n")
    tikz.write("xmode=log,\n")
    if legend:
        tikz.write("legend to name=grouplegend,\n")
    tikz.write("]\n")
    if legend:
        for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
            for init_type in ["MAX", "MIN", "MEAN", "MEDOID", "RANDOM"]:
                style = "only marks, mark size=3pt, thick, " + marks[init_type] + ", " + colors[algo]
                tikz.write("\\addlegendimage{" + style + "};\n")
                tikz.write("\\addlegendentry{\\" + macros[algo] + " (" + init_type.lower() + " init., w/ tightening)};\n")
                style = "only marks, mark size=3pt, " + marks[init_type + "conv"] + ", " + colors[algo]
                tikz.write("\\addlegendimage{" + style + "};\n")
                tikz.write("\\addlegendentry{\\" + macros[algo] + " (" + init_type.lower() + " init., w/o tightening)};\n")
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        for init_type in ["MAX", "MIN", "MEAN", "MEDOID", "RANDOM"]:
            num_inits_str = ""
            if init_type == "RANDOM":
                num_inits_str = str(num_inits[init_type + algo])
            else:
                num_inits_str = str(num_inits[init_type])
            col_prefix = algo + "_" + init_type + "_I" + num_inits_str + "_"
            sod = rs_vs_percent.at[index, col_prefix + "sod"]
            time = rs_vs_percent.at[index, col_prefix + "time"]
            style = "only marks, mark size=3pt, thick, " + marks[init_type] + ", " + colors[algo]
            tikz.write("\\addplot[" + style + "] coordinates {(" + str(time) + "," + str(sod) + ")};\n")
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        for init_type in ["MAX", "MIN", "MEAN", "MEDOID", "RANDOM"]:
            num_inits_str = ""
            if init_type == "RANDOM":
                num_inits_str = str(num_inits[init_type + algo])
            else:
                num_inits_str = str(num_inits[init_type])
            col_prefix = algo + "_" + init_type + "_I" + num_inits_str + "_"
            style = "only marks, mark size=3pt, thick,  " + marks[init_type + "conv"] + ", " + colors[algo]
            sod = rs_vs_percent.at[index, col_prefix + "sod_converged"]
            time = rs_vs_percent.at[index, col_prefix + "time_converged"]
            tikz.write("\\addplot[" + style + "] coordinates {(" + str(time) + "," + str(sod) + ")};\n")
            
def create_lineplots_max_itrs_vs_percent(rs_vs_percent, tikz, title, legend = False):
    marks = {"MAX" : "mark=o", "MIN" : "mark=square", "MEAN" : "mark=triangle", "MEDOID" : "mark=diamond", "RANDOM" : "mark=pentagon", "MAXconv" : "mark=+", "MINconv" : "mark=x", "MEANconv" : "mark=-", "MEDOIDconv" : "mark=|", "RANDOMconv" : "mark=star"}
    colors = {"BRANCH_FAST" : "draw=BRANCHFAST", "REFINE" : "draw=REFINE", "IPFP" : "draw=IPFP"} 
    num_inits = {"MAX" : 1, "MIN" : 1, "MEAN" : 1, "MEDOID" : 1, "RANDOMBRANCH_FAST" : 8, "RANDOMREFINE" : 2, "RANDOMIPFP" : 2}
    macros = {"BRANCH_FAST" : "BRANCHFAST", "REFINE" : "REFINE", "IPFP" : "IPFP"}
    tikz.write("\\nextgroupplot[\n")
    tikz.write("title=" + title + ",\n")
    tikz.write("xtick={10,20,30,40,50,60,70,80,90},\n")
    tikz.write("xticklabels={10,20,30,40,50,60,70,80,90},\n")
    tikz.write("xlabel={\\si{\\percent} of graphs in collection},\n")
    tikz.write("ylabel={max \# iterations},\n")
    tikz.write("ylabel style={yshift=-.1cm},\n")
    tikz.write("xlabel style={yshift=.1cm},\n")
    if legend:
        tikz.write("legend to name=grouplegend,\n")
    tikz.write("]\n")
    if legend:
        for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
            for init_type in ["MAX", "MIN", "MEAN", "MEDOID", "RANDOM"]:
                style = "only marks, mark size=3pt, thick, " + marks[init_type] + ", " + colors[algo]
                tikz.write("\\addlegendimage{" + style + "};\n")
                tikz.write("\\addlegendentry{\\" + macros[algo] + " (" + init_type.lower() + " init., w/ tightening)};\n")
                style = "only marks, mark size=3pt, " + marks[init_type + "conv"] + ", " + colors[algo]
                tikz.write("\\addlegendimage{" + style + "};\n")
                tikz.write("\\addlegendentry{\\" + macros[algo] + " (" + init_type.lower() + " init., w/o tightening)};\n")
    
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        for init_type in ["MAX", "MIN", "MEAN", "MEDOID", "RANDOM"]:
            num_inits_str = ""
            if init_type == "RANDOM":
                num_inits_str = str(num_inits[init_type + algo])
            else:
                num_inits_str = str(num_inits[init_type])
            col = algo + "_" + init_type + "_I" + num_inits_str + "_itrs_max"
            style = "mark size=3pt, " + marks[init_type + "conv"] + ", " + colors[algo]
            tikz.write("\\addplot[" + style + "] table [x=x, y=y, col sep=comma] {\n")
            tikz.write("x,y\n")
            for pos in range(9):
                tikz.write(str(10 * (pos + 1)) + "," + str(rs_vs_percent.at[pos, col]) + "\n")
            tikz.write("};\n")

def create_lineplots_inits_vs_sod(rs_vs_num_inits, tikz, title, legend = False):
    tikz.write("\\nextgroupplot[\n")
    tikz.write("title=" + title + ",\n")
    tikz.write("xtick={1,2,4,8,16,32},\n")
    tikz.write("xticklabels={1,2,4,8,16,32},\n")
    tikz.write("xlabel={initial solutions},\n")
    tikz.write("ylabel={mean \SOD at termination},\n")
    tikz.write("ylabel style={yshift=-.1cm},\n")
    tikz.write("xlabel style={yshift=.1cm},\n")
    tikz.write("xmode=log,\n")
    if legend:
        tikz.write("legend to name=grouplegend,\n")
    tikz.write("]\n")
    xshift = {"BRANCH_FAST" : -1, "REFINE" : 0, "IPFP" : 1}
    color = {"BRANCH_FAST" : "BRANCHFAST", "REFINE" : "REFINE", "IPFP" : "IPFP"}
    index = {2 ** i : i for i in range(6)}
    pos = {2 ** i : 4 * i + 2 for i in range(6)}
    if legend:
        for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
            tikz.write("\\addlegendimage{fill=" + color[algo] + ", area legend};\n")
            tikz.write("\\addlegendentry{\\" + color[algo] + " (converged)};\n")
        for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
            tikz.write("\\addlegendimage{fill=" + color[algo] + ", postaction={pattern=north east lines}, area legend};\n")
            tikz.write("\\addlegendentry{\\" + color[algo] + " (tightened)};\n")
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        col = algo + "_RANDOM_P90_sod";
        tikz.write("\\addplot[color=" + color[algo] + ", thick] table [x=x, y=y, col sep=comma] {\n")
        tikz.write("x,y\n")
        for num_inits in [2 ** i for i in range(6)]:
            tikz.write(str(num_inits) + "," + str(rs_vs_num_inits.at[index[num_inits], col]) + "\n")
        tikz.write("};\n")
        
def create_lineplots_inits_vs_max_itrs(rs_vs_num_inits, tikz, title, legend = False):
    tikz.write("\\nextgroupplot[\n")
    tikz.write("title=" + title + ",\n")
    tikz.write("xtick={1,2,4,8,16,32},\n")
    tikz.write("xticklabels={1,2,4,8,16,32},\n")
    tikz.write("xlabel={initial solutions},\n")
    tikz.write("ylabel={max \# iterations},\n")
    tikz.write("ylabel style={yshift=-.1cm},\n")
    tikz.write("xlabel style={yshift=.1cm},\n")
    tikz.write("xmode=log,\n")
    if legend:
        tikz.write("legend to name=grouplegend,\n")
    tikz.write("]\n")
    xshift = {"BRANCH_FAST" : -1, "REFINE" : 0, "IPFP" : 1}
    color = {"BRANCH_FAST" : "BRANCHFAST", "REFINE" : "REFINE", "IPFP" : "IPFP"}
    index = {2 ** i : i for i in range(6)}
    pos = {2 ** i : 4 * i + 2 for i in range(6)}
    if legend:
        for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
            tikz.write("\\addlegendimage{fill=" + color[algo] + ", area legend};\n")
            tikz.write("\\addlegendentry{\\" + color[algo] + " (converged)};\n")
            tikz.write("\\addlegendimage{fill=" + color[algo] + ", postaction={pattern=north east lines}, area legend};\n")
            tikz.write("\\addlegendentry{\\" + color[algo] + " (tightened)};\n")
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        col = algo + "_RANDOM_P90_itrs_max";
        tikz.write("\\addplot[color=" + color[algo] + ", thick] table [x=x, y=y, col sep=comma] {\n")
        tikz.write("x,y\n")
        for num_inits in [2 ** i for i in range(6)]:
            tikz.write(str(num_inits) + "," + str(rs_vs_num_inits.at[index[num_inits], col]) + "\n")
        tikz.write("};\n")
        
def create_lineplots_inits_vs_time(rs_vs_num_inits, tikz, title):
    tikz.write("\\nextgroupplot[\n")
    tikz.write("title=" + title + ",\n")
    tikz.write("xtick={1,2,4,8,16,32},\n")
    tikz.write("xticklabels={1,2,4,8,16,32},\n")
    tikz.write("xlabel={initial solutions},\n")
    tikz.write("ylabel={\SOD},\n")
    tikz.write("ylabel style={yshift=-.1cm},\n")
    tikz.write("xlabel style={yshift=.1cm},\n")
    tikz.write("xmode=log,\n")
    tikz.write("]\n")
    xshift = {"BRANCH_FAST" : -1, "REFINE" : 0, "IPFP" : 1}
    color = {"BRANCH_FAST" : "BRANCHFAST", "REFINE" : "REFINE", "IPFP" : "IPFP"}
    index = {2 ** i : i for i in range(6)}
    pos = {2 ** i : 4 * i + 2 for i in range(6)}
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        col = algo + "_RANDOM_P90_time";
        tikz.write("\\addplot[color=" + color[algo] + ", thick] table [x=x, y=y, col sep=comma] {\n")
        tikz.write("x,y\n")
        for num_inits in [2 ** i for i in range(6)]:
            tikz.write(str(num_inits) + "," + str(rs_vs_num_inits.at[index[num_inits], col]) + "\n")
        tikz.write("};\n")
    
def create_barplots_inits_vs_state(rs_vs_num_inits, tikz, title):
    tikz.write("\\nextgroupplot[\n")
    tikz.write("title=" + title + ",\n")
    tikz.write("xtick={2,6,10,14,18,22},\n")
    tikz.write("xticklabels={1,2,4,8,16,32},\n")
    tikz.write("xlabel={initial solutions},\n")
    tikz.write("ylabel={termination state in \\si{\\percent}},\n")
    tikz.write("ylabel style={yshift=-.1cm},\n")
    tikz.write("ymin=0,\n")
    tikz.write("xlabel style={yshift=.1cm},\n")
    tikz.write("cycle list={},\n")
    tikz.write("bar width=0.1cm,\n")
    tikz.write("]\n")
    xshift = {"BRANCH_FAST" : -1.1, "REFINE" : 0, "IPFP" : 1.1}
    color = {"BRANCH_FAST" : "BRANCHFAST", "REFINE" : "REFINE", "IPFP" : "IPFP"}
    index = {2 ** i : i for i in range(6)}
    pos = {2 ** i : 4 * i + 2 for i in range(6)}
    for algo in ["BRANCH_FAST", "REFINE", "IPFP"]:
        col_prefix = algo + "_RANDOM_P90_ratio";
        for num_inits in [2 ** i for i in range(6)]:
            ratio_refined = rs_vs_num_inits.at[index[num_inits], col_prefix + "_refined"] * 100
            ratio_converged = rs_vs_num_inits.at[index[num_inits], col_prefix + "_converged"] * 100
            ratio_init = rs_vs_num_inits.at[index[num_inits], col_prefix + "_initialized"] * 100
            if (ratio_converged > 0):
                tikz.write("\\addplot+[ybar,\n")
                tikz.write("fill=" + color[algo] + ",\n")
                tikz.write("draw=none,\n")
                tikz.write("] coordinates {(" + str(pos[num_inits] + xshift[algo]) +  "," + str(ratio_converged) + ")};\n")
            if (ratio_refined > 0):
                tikz.write("\\addplot+[ybar,\n")
                tikz.write("fill=" + color[algo] + ", postaction={pattern=north east lines},\n")
                tikz.write("draw=none,\n")
                tikz.write("] coordinates {(" + str(pos[num_inits] + xshift[algo]) +  "," + str(ratio_refined) + ")};\n")

agg_rs_letter = aggregate_results("Letter_RESULTS.csv")
agg_rs_aids = aggregate_results("AIDS_RESULTS.csv")
agg_rs_muta = aggregate_results("Mutagenicity_RESULTS.csv")
rs_vs_num_inits_letter = compute_rs_vs_num_inits(agg_rs_letter)
rs_vs_num_inits_aids = compute_rs_vs_num_inits(agg_rs_aids)
rs_vs_num_inits_muta = compute_rs_vs_num_inits(agg_rs_muta)
create_init_plots(rs_vs_num_inits_letter, rs_vs_num_inits_aids, rs_vs_num_inits_muta)
rs_vs_percent_letter = compute_rs_vs_percent(agg_rs_letter)
rs_vs_percent_aids = compute_rs_vs_percent(agg_rs_aids)
rs_vs_percent_muta = compute_rs_vs_percent(agg_rs_muta)
create_config_plots(rs_vs_percent_letter, rs_vs_percent_aids, rs_vs_percent_muta)
