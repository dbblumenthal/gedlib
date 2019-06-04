#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2018 by David B. Blumenthal                              #
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
# @brief Processes results of experiments for VLDB Journal paper.
#
# @details This file was used for visualizing the experiments in the following paper:
# - D. B. Blumenthal, N. Boria, J. Gamper, S. Bougleux, L. Brun:
#   &ldquo;Comparing heuristics for graph edit distance computation&rdquo;,
#   VLDB J. 2019
#
# Usage: 
# ```sh
# $ python install.py [--help] [-h] \<tikz_dir\> \<table_dir\> \<data-dir\> [--no_lsape] [--no_ls] [--no_lp] [--no_misc]
# ```
#
# For more information, execute `$ python process_results.py --help`.

'''Processes results of experiments for VLDB Journal paper.'''

import csv
from pickle import NONE
import argparse
from decimal import Decimal
import os.path

def computes_no_lb(method_name, method_config=""):
    if method_name == "BP":
        return True
    elif method_name == "SUBGRAPH":
        return True
    elif method_name == "WALKS":
        return True
    elif method_name == "RINGOPT":
        return True
    elif method_name == "RINGMS":
        return True
    elif method_name == "RINGMLDNN":
        return True
    elif method_name == "RINGMLSVM":
        return True
    elif method_name == "PREDICTDNN":
        return True
    elif method_name == "PREDICTSVM":
        return True
    elif method_name == "REFINE":
        return True
    elif method_name == "KREFINE":
        return True
    elif method_name == "BPBEAM":
        return True
    elif method_name == "IBPBEAM":
        return True
    elif method_name == "IPFP":
        return True
    elif method_name == "SA":
        return True
    elif (method_config != "") and (int(method_config[2:-2].split(",")[0]) > 1):
        return True
    else:
        return False

def computes_no_ub(method_name):
    if method_name == "HED":
        return True
    elif method_name == "BRANCHCOMPACT":
        return True
    elif method_name == "PARTITION":
        return True
    elif method_name == "HYBRID":
        return True
    else:
        return False

def is_ls_based(method):
    if method.name == "REFINE":
        return True
    elif method.name == "KREFINE":
        return True
    elif method.name == "BPBEAM":
        return True
    elif method.name == "IBPBEAM":
        return True
    elif method.name == "IPFP":
        return True
    else:
        return False

def is_lsape_based(method):
    if method.name == "BP":
        return True
    elif method.name == "BRANCH":
        return True
    elif method.name == "BRANCHFAST":
        return True
    elif method.name == "BRANCHUNI":
        return True
    elif method.name == "STAR":
        return True
    elif method.name == "NODE":
        return True
    elif method.name == "SUBGRAPH":
        return True
    elif method.name == "WALKS":
        return True
    elif method.name == "RINGOPT":
        return True
    elif method.name == "RINGMS":
        return True
    elif method.name == "RINGMLDNN":
        return True
    elif method.name == "RINGMLSVM":
        return True
    elif method.name == "PREDICTDNN":
        return True
    elif method.name == "PREDICTSVM":
        return True
    else:
        return False
    
def uses_randpost(method):
    if is_ls_based(method):
        return (int(method.config[2:-2].split(",")[2]) > 0)
    else:
        return False
    
def uses_multi_start(method):
    if is_ls_based(method):
        return (int(method.config[2:-2].split(",")[0]) > 1)
    else:
        return False

def uses_multi_sol(method):
    if is_lsape_based(method):
        return (int(method.config[2:-2].split(",")[0]) > 1)
    else:
        return False

def uses_centralities(method):
    if is_lsape_based(method):
        return (float(method.config[2:-2].split(",")[1]) > 0)
    else:
        return False

class Method:
    
    def __init__(self, name, lb, ub, t, coeff_lb, coeff_ub):
        self.consider_lb = True
        self.name = name[0]
        if self.name == "SUBGGRAPH":
            self.name = "SUBGRAPH"
        self.config = name[1]
        self.lb = float("{0:.2f}".format(Decimal(lb)))
        self.ub = float("{0:.2f}".format(Decimal(ub)))
        self.t = float("{:.6f}".format(Decimal(t)))
        self.coeff_lb = float("{:.2f}".format(Decimal(coeff_lb)))
        self.precise_coeff_lb = coeff_lb
        self.coeff_ub = float("{:.2f}".format(Decimal(coeff_ub)))
        self.precise_coeff_ub = coeff_ub
        self.is_fastest_lb = not computes_no_lb(self.name, self.config)
        self.is_fastest_ub = not computes_no_ub(self.name)
        self.has_tightest_lb = not computes_no_lb(self.name, self.config)
        self.has_tightest_ub = not computes_no_ub(self.name)
        self.has_best_coeff_lb = not computes_no_lb(self.name, self.config)
        self.has_best_coeff_ub = not computes_no_ub(self.name)
        self.is_maximum_lb = not computes_no_lb(self.name, self.config)
        self.is_maximum_ub = not computes_no_ub(self.name)
        self.discard_for_lb = computes_no_lb(self.name, self.config)
        self.discard_for_ub = computes_no_ub(self.name)
        self.score_lb = 0
        self.score_ub = 0
        self.adj_list_lb = []
        self.adj_list_ub = []
            
    def stats(self):
        method_stats = "\\begin{tiny}$\left(\\begin{smallmatrix}\\\\"
        if self.consider_lb:
            method_stats = method_stats + "t\\text{ in \si{\second}} & d_{\LB} & c_{\LB} & s_{LB}\\\\"
        else:
            method_stats = method_stats + "t\\text{ in \si{\second}} & d_{\UB} & c_{\UB} & s_{UB}\\\\"
        method_stats = method_stats + "\\num{" + "{:.2E}".format(Decimal(str(self.t))) + "} & "
        if self.consider_lb:
            method_stats = method_stats + "\\num{" + str(self.lb) + "} & "
            method_stats = method_stats + "\\num{" + str(self.coeff_lb) + "} & "
            method_stats = method_stats + "\\num{" + "{:.2f}".format(Decimal(str(self.score_lb))) + "}"
        else:
            method_stats = method_stats + "\\num{" + str(self.ub) + "} &"
            method_stats = method_stats + "\\num{" + str(self.coeff_ub) + "} &"
            method_stats = method_stats + "\\num{" + "{:.2f}".format(Decimal(str(self.score_ub))) + "}"
        method_stats = method_stats + "\\end{smallmatrix}\\right)$\end{tiny}"
        return method_stats
    
    def tikz_descriptor(self):
        descriptor = "\\" + self.name
        if self.is_maximum() and (self.config != ""):
            descriptor = descriptor + " " + self.config
        if not self.is_maximum():
            descriptor = descriptor + "\\\\" + self.config
        if self.consider_lb and self.is_maximum_lb:
            descriptor = descriptor + "\\\\" + self.stats()
        if (not self.consider_lb) and self.is_maximum_ub:
            descriptor = descriptor + "\\\\" + self.stats()
        return descriptor
    
    def label(self):
        labels = []
        if self.consider_lb and self.has_tightest_lb:
            labels.append("\\textcolor{Blue}{$d^\star_{\LB}$}")
        if (not self.consider_lb) and self.has_tightest_ub:
            labels.append("\\textcolor{Blue}{$d^\star_{\UB}$}")
        if self.consider_lb and self.is_fastest_lb:
            labels.append("\\textcolor{Red}{$t^\star_{\LB}$}")
        if (not self.consider_lb) and self.is_fastest_ub:
            labels.append("\\textcolor{Red}{$t^\star_{\UB}$}")
        if self.consider_lb and self.has_best_coeff_lb:
            labels.append("\\textcolor{Green}{$c^\star_{\LB}$}")
        if (not self.consider_lb) and self.has_best_coeff_ub:
            labels.append("\\textcolor{Green}{$c^\star_{\UB}$}")
        if len(labels) == 0:
            return ""
        label = labels[0]
        for index in range(1, len(labels)):
            label = label + " \\\\ " + labels[index]
        return label
    
    def compare_tightness(self, other):
        if self.consider_lb:
            if self.lb > other.lb:
                return 1
            elif self.lb == other.lb:
                return 0
            else:
                return -1
        else:
            if self.ub < other.ub:
                return 1
            elif self.ub == other.ub:
                return 0
            else:
                return -1
    
    def compare_time(self, other):
        if self.t < other.t:
            return 1
        elif self.t == other.t:
            return 0
        else:
            return -1
        
    def compare_coeff(self, other):
        if self.consider_lb:
            if self.coeff_lb > other.coeff_lb:
                return 1
            elif self.coeff_lb == other.coeff_lb:
                return 0
            else:
                return -1
        else:
            if self.coeff_ub > other.coeff_ub:
                return 1
            elif self.coeff_ub == other.coeff_ub:
                return 0
            else:
                return -1
    
    def get_edge_label(self, other):
        is_better_or_equal = True
        if self.compare_tightness(other) < 0:
            is_better_or_equal = False
            if self.consider_lb:
                self.has_tightest_lb = False
            else:
                self.has_tightest_ub = False
        if self.compare_time(other) < 0:
            is_better_or_equal = False
            if self.consider_lb:
                self.is_fastest_lb = False
            else:
                self.is_fastest_ub = False
        if self.compare_coeff(other) < 0:
            is_better_or_equal = False
            if self.consider_lb:
                self.has_best_coeff_lb = False
            else:
                self.has_best_coeff_ub = False
        label = ""
        if self.compare_tightness(other) > 0:
            if self.consider_lb:
                other.has_tightest_lb = False
            else:
                other.has_tightest_ub = False
            label = label + "d"
        if self.compare_time(other) > 0:
            if self.consider_lb:
                other.is_fastest_lb = False
            else:
                other.is_fastest_ub = False
            label = label + "t"
        if self.compare_coeff(other) > 0:
            if self.consider_lb:
                other.has_best_coeff_lb = False
            else:
                other.has_best_coeff_ub = False
            label = label + "c"
        if is_better_or_equal and (label != ""):
            if self.consider_lb:
                other.is_maximum_lb = False
            else:
                other.is_maximum_ub = False
        if is_better_or_equal:
            return label
        else:
            return ""
    
    def as_table_row(self):
        table_row = "\\" + self.name
        if self.is_maximum_lb:
            table_row = table_row + " $\LB^\star$"
        if self.is_maximum_ub:
            table_row = table_row + " $\UB^\star$"
        if self.config != "":
            table_row = table_row + " & " + self.config
        else:
            table_row = table_row + " & {--}"
        if not computes_no_lb(self.name):
            table_row = table_row + " & " + str(self.lb)
        else:
            table_row = table_row + " & {--}"
        if not computes_no_ub(self.name):
            table_row = table_row + " & " + str(self.ub)
        else:
            table_row = table_row + " & {--}"
        table_row = table_row + " & " + "{:.2E}".format(Decimal(str(self.t)))
        if not computes_no_lb(self.name):
            table_row = table_row + " & " + str(self.coeff_lb)
        else:
            table_row = table_row + " & {--}"
        if not computes_no_ub(self.name):
            table_row = table_row + " & " + str(self.coeff_ub)
        else:
            table_row = table_row + " & {--}"
        if not computes_no_lb(self.name):
            table_row = table_row + " & " + "{:.2}".format(Decimal(str(self.score_lb)))
        else:
            table_row = table_row + " & {--}"
        if not computes_no_ub(self.name):
            table_row = table_row + " & " + "{:.2}".format(Decimal(str(self.score_ub)))
        else:
            table_row = table_row + " & {--}"
        table_row = table_row + " \\\\\n"
        return table_row
    
    def is_maximum(self):
        if self.consider_lb:
            return self.is_maximum_lb
        else:
            return self.is_maximum_ub
    
    def dist(self):
        if self.consider_lb:
            return self.lb
        else:
            return self.ub
        
    def score(self):
        if self.consider_lb:
            return self.score_lb
        else:
            return self.score_ub
    
    def coeff(self):
        if self.consider_lb:
            return self.coeff_lb
        else:
            return self.coeff_ub
    
    def has_tightest_dist(self):
        if self.consider_lb:
            return self.has_tightest_lb
        else:
            return self.has_tightest_ub
        
    def is_fastest(self):
        if self.consider_lb:
            return self.is_fastest_lb
        else:
            return self.is_fastest_ub
        
    def has_best_coeff(self):
        if self.consider_lb:
            return self.has_best_coeff_lb
        else:
            return self.has_best_coeff_ub
    
    def discard(self):
        if self.consider_lb:
            return self.discard_for_lb
        else:
            return self.discard_for_ub
        
    def do_discard(self):
        if self.consider_lb:
            self.discard_for_lb = True 
        else:
            self.discard_for_ub = True
    
    def get_adj_list(self):
        if self.consider_lb:
            return self.adj_list_lb
        else:
            return self.adj_list_ub
    
    def set_adj_list(self, new_adj_list):
        if self.consider_lb:
            self.adj_list_lb = new_adj_list
        else:
            self.adj_list_ub = new_adj_list
    
    def set_score(self, best_dist, best_t, best_coeff):
        if self.consider_lb:
            self.score_lb = ((self.lb / best_dist) + (best_t / self.t) + (self.coeff_lb / best_coeff)) / 3.0
        else:
            self.score_ub = ((best_dist / self.ub) + (best_t / self.t) + (self.coeff_ub / best_coeff)) / 3.0
            
def parse_method_name(method_name):
    method_name_list = method_name.split(",", 1)
    if (len(method_name_list) == 1):
        method_name_list.append("")
    return method_name_list

def dfs(methods, is_discarded_edge, parent_id, child_id, seen):
    if seen[child_id]:
        return
    for edge in methods[child_id].get_adj_list():
        is_discarded_edge[parent_id][edge[0]] = True;
        dfs(methods, is_discarded_edge, parent_id, edge[0], seen)
    seen[child_id] = True
    
def read_results_from_csv_files(dataset, args):
    methods = []
    result_file_names = []
    prefix = os.path.join("results", dataset) + "__"
    if not args.no_lsape:
        result_file_names.append(prefix + "lsape_based_methods.csv")
    if not args.no_ls:
        result_file_names.append(prefix + "ls_based_methods.csv")
    if not args.no_lp:
        result_file_names.append(prefix + "lp_based_methods.csv")
    if not args.no_misc:
        result_file_names.append(prefix + "misc_methods.csv")
    for result_file_name in result_file_names:
        with open(result_file_name, "r") as result_file:  
            csv_reader = csv.reader(result_file,delimiter=";")
            next(csv_reader, NONE)
            for row in csv_reader:
                methods.append(Method(parse_method_name(row[0]), row[1], row[2], row[3], row[4], row[5]))
    return methods

def build_dependency_graph(methods, consider_lb):
    # set consider_lb for all methods
    for method in methods:
        method.consider_lb = consider_lb
    # construct dominance graph
    num_methods = len(methods)
    adj_list = [[] for method_id in range(0,num_methods)]
    for id_1 in range(0, num_methods):
        method_1 = methods[id_1]
        new_adj_list = []
        if not method_1.discard():
            for id_2 in range(0, num_methods):
                method_2 = methods[id_2]
                if not method_2.discard():
                    edge_label = method_1.get_edge_label(method_2)
                    if edge_label != "":
                        new_adj_list.append((id_2, edge_label))
        method_1.set_adj_list(new_adj_list)
    # compute scores    
    best_dist = 0
    best_t = 0
    best_coeff = 0
    for method in methods:
        if method.discard():
            continue
        if method.has_tightest_dist():
            best_dist = method.dist()
        if method.is_fastest():
            best_t = method.t
        if method.has_best_coeff():
            best_coeff = method.coeff()
    for method in methods:
        method.set_score(best_dist, best_t, best_coeff)
    # discard methods such that are not maximal w.r.t. their score within their heuristic
    method_names = set()
    for method in methods:
        method_names.add(method.name)
    best_scores = {method_name : -1 for method_name in method_names}
    for method in methods:
        if method.score() > best_scores[method.name]:
            best_scores[method.name] = method.score()
    for method in methods:
        if (method.score() < best_scores[method.name]) and (not method.is_maximum()):
            method.do_discard()
    # select undiscarded methods
    undiscarded_method_ids = []
    for id_1 in range(0, num_methods):
        method_1 = methods[id_1]
        new_adj_list = []
        if not method_1.discard():
            undiscarded_method_ids.append(id_1)
            old_adj_list = method_1.get_adj_list()
            for edge in old_adj_list:
                if not methods[edge[0]].discard():
                    new_adj_list.append(edge)
        method_1.set_adj_list(new_adj_list)
    # compute transitive reduction of undiscarded methods
    is_discarded_edge = {id_1 : {id_2 : False for id_2 in undiscarded_method_ids} for id_1 in undiscarded_method_ids}
    for method_id in undiscarded_method_ids:
        method = methods[method_id]
        seen = {child_id : False for child_id in undiscarded_method_ids}
        for edge in method.get_adj_list():
            dfs(methods, is_discarded_edge, method_id, edge[0], seen)
    # discard edges that are not contained in transitive reduction
    for id_1 in range(0, num_methods):
        method_1 = methods[id_1]
        new_adj_list = []
        if not method_1.discard():
            old_adj_list = method_1.get_adj_list()
            for edge in old_adj_list:
                if not is_discarded_edge[id_1][edge[0]]:
                    new_adj_list.append(edge)
        method_1.set_adj_list(new_adj_list)
    return methods

def infix(args):
    the_infix = ""
    if args.no_lsape:
        the_infix = the_infix + "_no-lsape"
    if args.no_ls:
        the_infix = the_infix + "_no-ls"
    if args.no_lp:
        the_infix = the_infix + "_no-lp"
    if args.no_misc:
        the_infix = the_infix + "_no-misc"
    return the_infix

class AggregatedScores:
    
    def __init__(self, method_names, scores_lb, scores_ub, has_tightest_lb, is_fastest_lb, has_best_coeff_lb, has_tightest_ub, is_fastest_ub, has_best_coeff_ub):
        self.method_names = method_names
        self.scores_lb = scores_lb
        self.scores_ub = scores_ub
        self.has_tightest_lb = has_tightest_lb
        self.is_fastest_lb = is_fastest_lb
        self.has_best_coeff_lb = has_best_coeff_lb
        self.has_tightest_ub = has_tightest_ub
        self.is_fastest_ub = is_fastest_ub
        self.has_best_coeff_ub = has_best_coeff_ub
    
    def chi_lb(self, method_name):
        chi = "$("
        if self.has_tightest_lb[method_name]:
            chi = chi + "\mathbf{1},"
        else:
            chi = chi + "0,"
        if self.is_fastest_lb[method_name]:
            chi = chi + "\mathbf{1},"
        else:
            chi = chi + "0,"
        if self.has_best_coeff_lb[method_name]:
            chi = chi + "\mathbf{1})$"
        else:
            chi = chi + "0)$"
        return chi
    
    def chi_ub(self, method_name):
        chi = "$("
        if self.has_tightest_ub[method_name]:
            chi = chi + "\mathbf{1},"
        else:
            chi = chi + "0,"
        if self.is_fastest_ub[method_name]:
            chi = chi + "\mathbf{1},"
        else:
            chi = chi + "0,"
        if self.has_best_coeff_ub[method_name]:
            chi = chi + "\mathbf{1})$"
        else:
            chi = chi + "0)$"
        return chi
    
    def write_to_csv_file(self, args):
        csv_file_name = os.path.join(args.table_dir, args.dataset) + infix(args) + "_scores.csv"
        csv_file = open(csv_file_name, "w")
        csv_file.write("heuristic;chi_lb;score_lb;chi_ub;score_ub\n")
        for method_name in self.method_names:
            csv_file.write(method_name + ";")
            if computes_no_lb(method_name):
                csv_file.write("na;na;")
            else:
                csv_file.write(self.chi_lb(method_name) + ";" + str(self.scores_lb[method_name]) + ";")
            if computes_no_ub(method_name):
                csv_file.write("na;na\n")
            else:
                csv_file.write(self.chi_ub(method_name) + ";" + str(self.scores_ub[method_name]) + "\n")
        ext_names = ["MULTISOL", "CENTRALITIES", "MULTISTART", "RANDPOST"]
        for ext_name in ext_names:
            csv_file.write(ext_name + ";na;na;")
            csv_file.write(self.chi_ub(ext_name) + ";" + str(self.scores_ub[ext_name]) + "\n")
        csv_file.close()

def create_latex_tables(args, datasets, aggregated_scores, lsape_based_method_names, lp_based_method_names, ls_based_method_names, misc_method_names, lsape_ext_names, ls_ext_names):
    sum_best_coeff_lb = {}
    table_ub_file_name = os.path.join(args.table_dir, "results_UB.tex")
    table_ub = open(table_ub_file_name, "w")
    table_ub.write("%!TEX root = ../root.tex\n")
    table_ub.write("\\begin{tabular}{@{}lcS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]@{}}\n")
    table_ub.write("\\toprule\n")
    table_ub.write("heuristic & \multicolumn{2}{c}{\\letter} & \multicolumn{2}{c}{\\mutagenicity} & \multicolumn{2}{c}{\\aids} & \multicolumn{2}{c}{\\protein} & \multicolumn{2}{c}{\\fingerprint} & \multicolumn{2}{c}{\\grec} \\\\\n")
    table_ub.write("\midrule\n")
    table_ub.write("& {$\chi_\UB$} & {$\widehat{s_\UB}$} & {$\chi_\UB$} & {$\widehat{s_\UB}$} & {$\chi_\UB$} & {$\widehat{s_\UB}$} & {$\chi_\UB$} & {$\widehat{s_\UB}$} & {$\chi_\UB$} & {$\widehat{s_\UB}$} & {$\chi_\UB$} & {$\widehat{s_\UB}$} \\\\\n")
    table_ub.write("\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7} \cmidrule(lr){8-9} \cmidrule(lr){10-11} \cmidrule(l){12-13}\n")
    table_ub.write("\multicolumn{13}{@{}l}{\emph{instantiations of the paradigm \LSAPEGED}} \\\\\n")
    for method_name in lsape_based_method_names:
        table_ub.write("\\" + method_name)
        for dataset in datasets:
            if aggregated_scores[dataset].scores_ub[method_name] > 0:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
            else:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
        table_ub.write(" \\\\\n")
    table_ub.write("\multicolumn{13}{@{}l}{\emph{extensions of the paradigm \LSAPEGED}} \\\\\n")
    for ext_name in lsape_ext_names:
        table_ub.write("\\" + ext_name)
        for dataset in datasets:
            if aggregated_scores[dataset].scores_ub[ext_name] > 0:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(ext_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[ext_name])))
            else:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(ext_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[ext_name])))
        table_ub.write(" \\\\\n")
    table_ub.write("\midrule\n")
    table_ub.write("\multicolumn{13}{@{}l}{\emph{instantiations of the paradigm \LPGED}} \\\\\n")
    for method_name in lp_based_method_names:
        table_ub.write("\\" + method_name)
        for dataset in datasets:
            if aggregated_scores[dataset].scores_ub[method_name] > 0:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
            else:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
        table_ub.write(" \\\\\n")
    table_ub.write("\midrule\n")
    table_ub.write("\multicolumn{13}{@{}l}{\emph{instantiations of the paradigm \LSGED}} \\\\\n")
    for method_name in ls_based_method_names:
        table_ub.write("\\" + method_name)
        for dataset in datasets:
            if aggregated_scores[dataset].scores_ub[method_name] > 0:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
            else:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
        table_ub.write(" \\\\\n")
    table_ub.write("\multicolumn{13}{@{}l}{\emph{extensions of the paradigm \LSGED}} \\\\\n")
    for ext_name in ls_ext_names:
        table_ub.write("\\" + ext_name)
        for dataset in datasets:
            if aggregated_scores[dataset].scores_ub[ext_name] > 0:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(ext_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[ext_name])))
            else:
                table_ub.write(" & " + aggregated_scores[dataset].chi_ub(ext_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[ext_name])))
        table_ub.write(" \\\\\n")
    table_ub.write("\midrule\n")
    table_ub.write("\multicolumn{13}{@{}l}{\emph{miscellaneous heuristics}} \\\\\n")
    for method_name in misc_method_names:
        if not computes_no_ub(method_name):
            table_ub.write("\\" + method_name)
            for dataset in datasets:
                if aggregated_scores[dataset].scores_ub[method_name] > 0:
                    table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
                else:
                    table_ub.write(" & " + aggregated_scores[dataset].chi_ub(method_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_ub[method_name])))
            table_ub.write(" \\\\\n")
    table_ub.write("\\bottomrule\n")
    table_ub.write("\end{tabular}\n")
    table_ub.close()
    table_lb_file_name = os.path.join(args.table_dir, "results_LB.tex")
    table_lb = open(table_lb_file_name, "w")
    table_lb.write("%!TEX root = ../root.tex\n")
    table_lb.write("\\begin{tabular}{@{}lcS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]cS[table-format=1.2]@{}}\n")
    table_lb.write("\\toprule\n")
    table_lb.write("heuristic & \multicolumn{2}{c}{\\letter} & \multicolumn{2}{c}{\\mutagenicity} & \multicolumn{2}{c}{\\aids} & \multicolumn{2}{c}{\\protein} & \multicolumn{2}{c}{\\fingerprint} & \multicolumn{2}{c}{\\grec} \\\\\n")
    table_lb.write("\midrule\n")
    table_lb.write("& {$\chi_\LB$} & {$\widehat{s_\LB}$} & {$\chi_\LB$} & {$\widehat{s_\LB}$} & {$\chi_\LB$} & {$\widehat{s_\LB}$} & {$\chi_\LB$} & {$\widehat{s_\LB}$} & {$\chi_\LB$} & {$\widehat{s_\LB}$} & {$\chi_\LB$} & {$\widehat{s_\LB}$} \\\\\n")
    table_lb.write("\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7} \cmidrule(lr){8-9} \cmidrule(lr){10-11} \cmidrule(l){12-13}\n")
    table_lb.write("\multicolumn{13}{@{}l}{\emph{instantiations of the paradigm \LSAPEGED}} \\\\\n")
    for method_name in lsape_based_method_names:
        if not computes_no_lb(method_name):
            table_lb.write("\\" + method_name)
            for dataset in datasets:
                if aggregated_scores[dataset].scores_lb[method_name]> 0:
                    table_lb.write(" & " + aggregated_scores[dataset].chi_lb(method_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_lb[method_name])))
                else:
                    table_lb.write(" & " + aggregated_scores[dataset].chi_lb(method_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_lb[method_name])))
            table_lb.write(" \\\\\n")
    table_lb.write("\midrule\n")
    table_lb.write("\multicolumn{13}{@{}l}{\emph{instantiations of the paradigm \LPGED}} \\\\\n")
    for method_name in lp_based_method_names:
        table_lb.write("\\" + method_name)
        for dataset in datasets:
            if aggregated_scores[dataset].scores_lb[method_name]> 0:
                table_lb.write(" & " + aggregated_scores[dataset].chi_lb(method_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_lb[method_name])))
            else:
                table_lb.write(" & " + aggregated_scores[dataset].chi_lb(method_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_lb[method_name])))
        table_lb.write(" \\\\\n")
    table_lb.write("\midrule\n")
    table_lb.write("\multicolumn{13}{@{}l}{\emph{miscellaneous heuristics}} \\\\\n")
    for method_name in misc_method_names:
        if not computes_no_lb(method_name):
            table_lb.write("\\" + method_name)
            for dataset in datasets:
                if aggregated_scores[dataset].scores_lb[method_name]> 0:
                    table_lb.write(" & " + aggregated_scores[dataset].chi_lb(method_name) + " & \\bfseries " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_lb[method_name])))
                else:
                    table_lb.write(" & " + aggregated_scores[dataset].chi_lb(method_name) + " & " + "{0:.2f}".format(Decimal(aggregated_scores[dataset].scores_lb[method_name])))
            table_lb.write(" \\\\\n")
    table_lb.write("\\bottomrule\n")
    table_lb.write("\end{tabular}\n")
    table_lb.close()
    
def create_pareto_data(args, dataset, methods, consider_lb, method_names):
    selected_method_names = [method_name for method_name in method_names if not computes_no_ub(method_name)]
    infix = "UB"
    bound_column_name = "avg_ub"
    runtime_column_name = "avg_runtime"
    if consider_lb:
        selected_method_names = [method_name for method_name in method_names if not computes_no_lb(method_name)]
        infix = "LB"
        bound_column_name = "avg_lb"
    max_score_methods = {method_name : None for method_name in selected_method_names}
    for method in methods:
        if not method.discard():
            if max_score_methods[method.name] is None:
                max_score_methods[method.name] = method
            elif method.score() > max_score_methods[method.name].score():
                max_score_methods[method.name] = method
    pareto_data_file_name = os.path.join(args.data_dir, dataset + "_Pareto_" + infix + ".csv")
    pareto_data_file = open(pareto_data_file_name, "w")
    for method_name in max_score_methods:
        pareto_data_file.write(method_name + "_" + runtime_column_name + "," + method_name + "_" + bound_column_name + ",")
    pareto_data_file.write("\n")
    for method_name in max_score_methods:
        pareto_data_file.write(str(max_score_methods[method_name].t) + "," + str(max_score_methods[method_name].dist()) + ",")
    pareto_data_file.write("\n")
    pareto_data_file.close()
        
def create_barplots(args, datasets, aggregated_scores, method_names, lsape_ext_names, ls_ext_names):
    lb_method_names = [method_name for method_name in method_names if not computes_no_lb(method_name)]
    ub_method_names = [method_name for method_name in method_names if not computes_no_ub(method_name)]
    # sum scores and chi for lower bounds
    scores_lb = {method_name : 0 for method_name in lb_method_names}
    has_tightest_lb = {method_name : 0 for method_name in lb_method_names}
    is_fastest_lb = {method_name : 0 for method_name in lb_method_names}
    has_best_coeff_lb = {method_name : 0 for method_name in lb_method_names}
    for method_name in lb_method_names:
        for dataset in datasets:
            scores_lb[method_name] = scores_lb[method_name] + aggregated_scores[dataset].scores_lb[method_name]
            if aggregated_scores[dataset].has_tightest_lb[method_name]:
                has_tightest_lb[method_name] = has_tightest_lb[method_name] + 1
            if aggregated_scores[dataset].is_fastest_lb[method_name]:
                is_fastest_lb[method_name] = is_fastest_lb[method_name] + 1 
            if aggregated_scores[dataset].has_best_coeff_lb[method_name]:
                has_best_coeff_lb[method_name] = has_best_coeff_lb[method_name] + 1
    scores_lb = [(k, v / 6.0) for k, v in sorted(scores_lb.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_tightest_lb = [(k, v) for k, v in sorted(has_tightest_lb.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    is_fastest_lb = [(k, v) for k, v in sorted(is_fastest_lb.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_best_coeff_lb = [(k, v) for k, v in sorted(has_best_coeff_lb.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    # sum scores and chi for upper bounds of heuristics
    scores_ub = {method_name : 0 for method_name in ub_method_names}
    has_tightest_ub = {method_name : 0 for method_name in ub_method_names}
    is_fastest_ub = {method_name : 0 for method_name in ub_method_names}
    has_best_coeff_ub = {method_name : 0 for method_name in ub_method_names}
    for method_name in ub_method_names:
        for dataset in datasets:
            scores_ub[method_name] = scores_ub[method_name] + aggregated_scores[dataset].scores_ub[method_name]
            if aggregated_scores[dataset].has_tightest_ub[method_name]:
                has_tightest_ub[method_name] = has_tightest_ub[method_name] + 1
            if aggregated_scores[dataset].is_fastest_ub[method_name]:
                is_fastest_ub[method_name] = is_fastest_ub[method_name] + 1 
            if aggregated_scores[dataset].has_best_coeff_ub[method_name]:
                has_best_coeff_ub[method_name] = has_best_coeff_ub[method_name] + 1
    scores_ub = [(k, v / 6.0) for k, v in sorted(scores_ub.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_tightest_ub = [(k, v) for k, v in sorted(has_tightest_ub.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    is_fastest_ub = [(k, v) for k, v in sorted(is_fastest_ub.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_best_coeff_ub = [(k, v) for k, v in sorted(has_best_coeff_ub.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    # sum scores and chi for upper bounds of LSAPE extensions
    scores_ub_lsape_ext = {method_name : 0 for method_name in lsape_ext_names}
    has_tightest_ub_lsape_ext = {method_name : 0 for method_name in lsape_ext_names}
    has_best_coeff_ub_lsape_ext = {method_name : 0 for method_name in lsape_ext_names}
    for ext_name in lsape_ext_names:
        for dataset in datasets:
            scores_ub_lsape_ext[ext_name] = scores_ub_lsape_ext[ext_name] + aggregated_scores[dataset].scores_ub[ext_name]
            if aggregated_scores[dataset].has_tightest_ub[ext_name]:
                has_tightest_ub_lsape_ext[ext_name] = has_tightest_ub_lsape_ext[ext_name] + 1
            if aggregated_scores[dataset].has_best_coeff_ub[ext_name]:
                has_best_coeff_ub_lsape_ext[ext_name] = has_best_coeff_ub_lsape_ext[ext_name] + 1
    scores_ub_lsape_ext = [(k, v / 6.0) for k, v in sorted(scores_ub_lsape_ext.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_tightest_ub_lsape_ext = [(k, v) for k, v in sorted(has_tightest_ub_lsape_ext.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_best_coeff_ub_lsape_ext = [(k, v) for k, v in sorted(has_best_coeff_ub_lsape_ext.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    # sum scores and chi for upper bounds of LS extensions
    scores_ub_ls_ext = {method_name : 0 for method_name in ls_ext_names}
    has_tightest_ub_ls_ext = {method_name : 0 for method_name in ls_ext_names}
    has_best_coeff_ub_ls_ext = {method_name : 0 for method_name in ls_ext_names}
    for ext_name in ls_ext_names:
        for dataset in datasets:
            scores_ub_ls_ext[ext_name] = scores_ub_ls_ext[ext_name] + aggregated_scores[dataset].scores_ub[ext_name]
            if aggregated_scores[dataset].has_tightest_ub[ext_name]:
                has_tightest_ub_ls_ext[ext_name] = has_tightest_ub_ls_ext[ext_name] + 1
            if aggregated_scores[dataset].has_best_coeff_ub[ext_name]:
                has_best_coeff_ub_ls_ext[ext_name] = has_best_coeff_ub_ls_ext[ext_name] + 1
    scores_ub_ls_ext = [(k, v / 6.0) for k, v in sorted(scores_ub_ls_ext.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_tightest_ub_ls_ext = [(k, v) for k, v in sorted(has_tightest_ub_ls_ext.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    has_best_coeff_ub_ls_ext = [(k, v) for k, v in sorted(has_best_coeff_ub_ls_ext.items(), key=lambda kv: kv[1], reverse=True) if v > 0]
    barplot_file_name = os.path.join(args.tikz_dir, "scores_LB.tex")
    write_barplot(barplot_file_name, "$\myavg_\mathcal{D}\widehat{s_\LB}$", 1.1, 1, scores_lb)
    barplot_file_name = os.path.join(args.tikz_dir, "has_tightest_LB.tex")
    write_barplot(barplot_file_name, "$(\sum_\mathcal{D}{\chi_\LB})_1$", 6.6, 0.5, has_tightest_lb)
    barplot_file_name = os.path.join(args.tikz_dir, "is_fastest_LB.tex")
    write_barplot(barplot_file_name, "$(\sum_\mathcal{D}{\chi_\LB})_2$", 6.6, 0.5, is_fastest_lb)
    barplot_file_name = os.path.join(args.tikz_dir, "has_best_coeff_LB.tex")
    write_barplot(barplot_file_name, "$(\sum_\mathcal{D}{\chi_\LB})_3$", 6.6, 1, has_best_coeff_lb)
    barplot_file_name = os.path.join(args.tikz_dir, "scores_UB.tex")
    write_barplot(barplot_file_name, "$\myavg_\mathcal{D}\widehat{s_\UB}$", 1.1, 1, scores_ub, scores_ub_lsape_ext, scores_ub_ls_ext)
    barplot_file_name = os.path.join(args.tikz_dir, "has_tightest_UB.tex")
    write_barplot(barplot_file_name, "$(\sum_\mathcal{D}{\chi_\UB})_1$", 6.6, 0.5, has_tightest_ub, has_tightest_ub_lsape_ext, has_tightest_ub_ls_ext)
    barplot_file_name = os.path.join(args.tikz_dir, "is_fastest_UB.tex")
    write_barplot(barplot_file_name, "$(\sum_\mathcal{D}{\chi_\UB})_2$", 6.6, 0.5, is_fastest_ub)
    barplot_file_name = os.path.join(args.tikz_dir, "has_best_coeff_UB.tex")
    write_barplot(barplot_file_name, "$(\sum_\mathcal{D}{\chi_\UB})_3$", 6.6, 1, has_best_coeff_ub, has_best_coeff_ub_lsape_ext, has_best_coeff_ub_ls_ext)
    
def write_barplot(barplot_file_name, ylabel, ymax, width, method_bars, lsape_ext_bars = [], ls_ext_bars = []):
    num_bars = len(method_bars + lsape_ext_bars + ls_ext_bars)
    xtick = ",".join(str(i * 0.5) for i in range(1, num_bars + 1))
    xticklabels = ",".join("\\" + kv[0] for kv in (method_bars + lsape_ext_bars + ls_ext_bars))
    barplot_file = open(barplot_file_name, "w")
    barplot_file.write("%!TEX root = ../root.tex\n")
    barplot_file.write("\\begin{tikzpicture}\n")
    barplot_file.write("\\begin{axis}[\n")
    barplot_file.write("height = 0.4\linewidth,\n")
    barplot_file.write("width = " + str(width) + "\linewidth,\n")
    barplot_file.write("xmin = 0,\n")
    barplot_file.write("xmax = " + str(0.5*(num_bars+1)) + ",\n")
    barplot_file.write("ymin = 0,\n")
    barplot_file.write("ymax = " + str(ymax) + ",\n")
    barplot_file.write("xtick = {" + xtick + "},\n")
    barplot_file.write("xticklabels = {" + xticklabels + "},\n")
    barplot_file.write("xticklabel style={rotate=67.5, anchor=east},\n")
    barplot_file.write("ylabel = " + ylabel + ",\n")
    barplot_file.write("every axis plot/.append style={ybar,bar width=.2, bar shift=0pt}]\n")
    xcoord = .5
    for kv in method_bars:
        barplot_file.write("\\addplot[" + kv[0] + "] coordinates {(" + str(xcoord) + "," + str(kv[1]) + ")};\n")
        xcoord = xcoord + .5
    if len(lsape_ext_bars) > 0:
        barplot_file.write("\draw[densely dotted] ({axis cs:" + str(xcoord - .25) + ",0}|-{rel axis cs:0,1}) -- ({axis cs:" + str(xcoord - .25) + ",0}|-{rel axis cs:0,0});\n")
    for kv in lsape_ext_bars:
        barplot_file.write("\\addplot[" + kv[0] + "] coordinates {(" + str(xcoord) + "," + str(kv[1]) + ")};\n")
        xcoord = xcoord + .5
    if len(ls_ext_bars) > 0:
        barplot_file.write("\draw[densely dotted] ({axis cs:" + str(xcoord - .25) + ",0}|-{rel axis cs:0,1}) -- ({axis cs:" + str(xcoord - .25) + ",0}|-{rel axis cs:0,0});\n")
    for kv in ls_ext_bars:
        barplot_file.write("\\addplot[" + kv[0] + "] coordinates {(" + str(xcoord) + "," + str(kv[1]) + ")};\n")
        xcoord = xcoord + .5
    barplot_file.write("\end{axis}\n")
    barplot_file.write("\end{tikzpicture}")
    barplot_file.close()
                  
def aggregate_scores(methods, method_names, method_ext_names):
    scores_lb = {method_name : 0.0 for method_name in method_names}
    scores_ub = {method_name : 0.0 for method_name in method_ext_names}
    has_tightest_lb = {method_name : False for method_name in method_names}
    is_fastest_lb = {method_name : False for method_name in method_names}
    has_best_coeff_lb = {method_name : False for method_name in method_names}
    has_tightest_ub = {method_name : False for method_name in method_ext_names}
    is_fastest_ub = {method_name : False for method_name in method_ext_names}
    has_best_coeff_ub = {method_name : False for method_name in method_ext_names}
    sum_scores_ub_lsape = 0.0
    sum_scores_ub_centralities = 0.0
    sum_scores_ub_multi_sol = 0.0
    sum_scores_ub_ls = 0.0
    sum_scores_ub_randpost = 0.0
    sum_scores_ub_multi_start = 0.0
    for method in methods:
        if method.is_maximum_lb:
            if method.score_lb > scores_lb[method.name]:
                scores_lb[method.name] = method.score_lb
        if method.is_maximum_ub:
            if method.score_ub > scores_ub[method.name]:
                scores_ub[method.name] = method.score_ub
            if is_lsape_based(method):
                sum_scores_ub_lsape = sum_scores_ub_lsape + method.score_ub
            if uses_centralities(method):
                sum_scores_ub_centralities = sum_scores_ub_centralities + method.score_ub
            if uses_multi_sol(method):
                sum_scores_ub_multi_sol = sum_scores_ub_multi_sol + method.score_ub
            if is_ls_based(method):
                sum_scores_ub_ls = sum_scores_ub_ls + method.score_ub
            if uses_multi_start(method):
                sum_scores_ub_multi_start = sum_scores_ub_multi_start + method.score_ub
            if uses_randpost(method):
                sum_scores_ub_randpost = sum_scores_ub_randpost + method.score_ub
        if method.has_tightest_lb:
            has_tightest_lb[method.name] = True
        if method.is_fastest_lb:
            is_fastest_lb[method.name] = True
        if method.has_best_coeff_lb:
            has_best_coeff_lb[method.name] = True
        if method.has_tightest_ub:
            has_tightest_ub[method.name] = True
            if uses_centralities(method):
                has_tightest_ub["CENTRALITIES"] = True
            if uses_multi_sol(method):
                has_tightest_ub["MULTISOL"] = True
            if uses_multi_start(method):
                has_tightest_ub["MULTISTART"] = True
            if uses_randpost(method):
                has_tightest_ub["RANDPOST"] = True
        if method.is_fastest_ub:
            is_fastest_ub[method.name] = True
            if uses_centralities(method):
                is_fastest_ub["CENTRALITIES"] = True
            if uses_multi_sol(method):
                is_fastest_ub["MULTISOL"] = True
            if uses_multi_start(method):
                is_fastest_ub["MULTISTART"] = True
            if uses_randpost(method):
                is_fastest_ub["RANDPOST"] = True
        if method.has_best_coeff_ub:
            has_best_coeff_ub[method.name] = True
            if uses_centralities(method):
                has_best_coeff_ub["CENTRALITIES"] = True
            if uses_multi_sol(method):
                has_best_coeff_ub["MULTISOL"] = True
            if uses_multi_start(method):
                has_best_coeff_ub["MULTISTART"] = True
            if uses_randpost(method):
                has_best_coeff_ub["RANDPOST"] = True
    if sum_scores_ub_lsape > 0:
        scores_ub["MULTISOL"] = sum_scores_ub_multi_sol / sum_scores_ub_lsape
        scores_ub["CENTRALITIES"] = sum_scores_ub_centralities / sum_scores_ub_lsape
    if sum_scores_ub_ls > 0:
        scores_ub["MULTISTART"] = sum_scores_ub_multi_start / sum_scores_ub_ls
        scores_ub["RANDPOST"] = sum_scores_ub_randpost / sum_scores_ub_ls
    return AggregatedScores(method_names, scores_lb, scores_ub, has_tightest_lb, is_fastest_lb, has_best_coeff_lb, has_tightest_ub, is_fastest_ub, has_best_coeff_ub)
                
def create_coeff_vs_dist_table(args, dataset, methods, consider_lb):
    table_file_name = os.path.join(args.data_dir, dataset) + infix(args)
    if consider_lb:
        table_file_name = table_file_name + "_LB.csv"
    else:
        table_file_name = table_file_name + "_UB.csv"
    table_file = open(table_file_name, "w")
    if consider_lb:
        table_file.write("avg_lb,coeff_lb\n")
    else:
        table_file.write("avg_ub,coeff_ub\n")
    for method in methods:
        if consider_lb:
            if not computes_no_lb(method.name, method.config):
                table_file.write(str(method.lb) + "," + str(method.coeff_lb) + "\n")
        else:
            if not computes_no_ub(method.name):
                table_file.write(str(method.ub) + "," + str(method.coeff_ub) + "\n")
    table_file.close()
    
def create_tikz_graph(args, dataset, methods, consider_lb):
    tikz_file_name = os.path.join(args.tikz_dir, dataset) + infix(args)
    if consider_lb:
        tikz_file_name = tikz_file_name + "_LB.tex"
    else:
        tikz_file_name = tikz_file_name + "_UB.tex"
    # construct tikz file
    tikz_file = open(tikz_file_name, "w")
    tikz_file.write("%!TEX root = ../root.tex\n")
    tikz_file.write("\\begin{tikzpicture}[rounded corners,every label/.style={align=center}]\n")
    tikz_file.write("\graph[layered layout,grow=right,\n")
    tikz_file.write("sibling distance=5pt,\n")
    tikz_file.write("part distance=5pt,\n")
    tikz_file.write("component distance=5pt,\n")
    tikz_file.write("sibling sep=5pt,\n")
    tikz_file.write("level sep=12pt,\n")
    tikz_file.write("part sep=5pt,\n")
    tikz_file.write("component sep=5pt,\n")
    tikz_file.write("component direction=up,\n")
    tikz_file.write("component align=first node,\n")
    tikz_file.write("nodes={minimum width=4mm, minimum height=4mm, align=center, font=\scriptsize},\n")
    tikz_file.write("edge quotes mid,\n")
    tikz_file.write("edges={nodes={font=\scriptsize, fill=white, inner sep=1.5pt}}] {\n")
    edge_id = 0
    tikz_file.write("{ [same layer]\n")
    for method_id in range(0, len(methods)):
        method = methods[method_id]
        if (not method.discard()) and method.is_maximum():
            tikz_file.write(str(method_id) + " [" + method.name + ", minimum width=2.5cm, as={" + method.tikz_descriptor() + "}, label=left:{\\footnotesize " + method.label() + "}],\n");
    tikz_file.write("};\n")
    for method_id in range(0, len(methods)):
        method = methods[method_id]
        if (not method.discard()) and (not method.is_maximum()):
            tikz_file.write(str(method_id) + " [" + method.name + ", as={" + method.tikz_descriptor() + "}];\n");
    for method_id in range(0, len(methods)):
        for edge in methods[method_id].get_adj_list():
            tikz_file.write(str(method_id) + " -> [" + edge[1] + "] " + str(edge[0]) + ";\n")
    tikz_file.write("};\n")
    tikz_file.write("\end{tikzpicture}")
    tikz_file.close()
       
parser = argparse.ArgumentParser(description="Generates TikZ dominance graph from CSV file.")
parser.add_argument("tikz_dir", help="name of output directory for TikZ files")
parser.add_argument("table_dir", help="name of output directory for LaTeX tables")
parser.add_argument("data_dir", help="name of output directory for csv tables")
parser.add_argument("--no_lsape", help="do not consider LSAPE based methods", action="store_true")
parser.add_argument("--no_ls", help="do not consider local search based methods", action="store_true")
parser.add_argument("--no_lp", help="do not consider LP based methods", action="store_true")
parser.add_argument("---no_misc", help="do not consider miscellaneous methods", action="store_true")

args = parser.parse_args()
datasets = ["Letter_HIGH", "Mutagenicity", "AIDS", "Protein", "Fingerprint", "GREC"]
aggregated_scores = {}
lsape_based_method_names = ["NODE", "BP", "BRANCH", "BRANCHFAST", "BRANCHUNI", "STAR", "SUBGRAPH", "WALKS", "RINGOPT", "RINGMS", "RINGMLSVM", "RINGMLDNN", "PREDICTSVM", "PREDICTDNN"]
lp_based_method_names = ["FONE", "FTWO", "COMPACTMIP", "JUSTICEIP"]
ls_based_method_names = ["REFINE", "KREFINE", "BPBEAM", "IBPBEAM", "IPFP"]
misc_method_names = ["HED", "BRANCHTIGHT", "SA", "BRANCHCOMPACT", "PARTITION", "HYBRID"]
lsape_ext_names = ["MULTISOL", "CENTRALITIES"]
ls_ext_names = ["MULTISTART", "RANDPOST"]
method_names = lsape_based_method_names + lp_based_method_names + ls_based_method_names + misc_method_names
method_ext_names = method_names + lsape_ext_names + ls_ext_names
for dataset in datasets:
    methods = read_results_from_csv_files(dataset, args)
    for consider_lb in [True, False]:
        methods = build_dependency_graph(methods, consider_lb)
        create_pareto_data(args, dataset, methods, consider_lb, method_names)
        create_tikz_graph(args, dataset, methods, consider_lb)
        create_coeff_vs_dist_table(args, dataset, methods, consider_lb)
    aggregated_scores[dataset] = aggregate_scores(methods, method_names, method_ext_names)
create_latex_tables(args, datasets, aggregated_scores, lsape_based_method_names, lp_based_method_names, ls_based_method_names, misc_method_names, lsape_ext_names, ls_ext_names)
create_barplots(args, datasets, aggregated_scores, method_names, lsape_ext_names, ls_ext_names)
