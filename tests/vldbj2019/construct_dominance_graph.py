import csv
from pickle import NONE
import argparse
from decimal import Decimal

class Method:
    
    def __init__(self, consider_lb, name, lb, ub, t, coeff_lb, coeff_ub):
        self.consider_lb = consider_lb
        self.name = name[0]
        self.config = name[1]
        self.lb = float("{:.2E}".format(Decimal(lb)))
        self.ub = float("{:.2E}".format(Decimal(ub)))
        self.t = float("{:.2E}".format(Decimal(t)))
        self.coeff_lb = float("{:.2E}".format(Decimal(coeff_lb)))
        self.coeff_ub = float("{:.2E}".format(Decimal(coeff_ub)))
        self.is_fastest = True
        self.is_tightest = True
        self.has_best_coeff = True
        self.is_maximum = True
        self.discard = False
            
    def stats(self):
        method_stats = ""
        if self.consider_lb:
            method_stats = method_stats + "$d=\\num{" + "{:.2E}".format(Decimal(str(self.lb))) + "}$\\\\"
        else:
            method_stats = method_stats + "$d=\\num{" + "{:.2E}".format(Decimal(str(self.ub))) + "}$\\\\"
        method_stats = method_stats + "$t=\\num{" + "{:.2E}".format(Decimal(str(self.t))) + "}$\\\\"
        if self.consider_lb:
            method_stats = method_stats + "$c=\\num{" + "{:.2}".format(Decimal(str(self.coeff_lb))) + "}$"
        else:
            method_stats = method_stats + "$c=\\num{" + "{:.2}".format(Decimal(str(self.coeff_ub))) + "}$"
        return method_stats
    
    def tikz_descriptor(self):
        if not self.is_maximum:
            return self.config
        else:
            if self.config == "":
                return self.stats()
            else:
                return (self.config + "\\\\" + self.stats())    
    
    def label(self):
        label = ""
        if self.is_tightest:
            label = "\\textcolor{blue}{$d$}"
        if self.is_fastest:
            if label == "":
                label = "\\textcolor{red}{$t$}"
            else:
                label = label + ",\\textcolor{red}{$t$}"
        if self.has_best_coeff:
            if label == "":
               label = "\\textcolor{green}{$c$}"
            else:
                label = label + ",\\textcolor{green}{$c$}"
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
            
    def __lt__(self, other):
        if self.compare_tightness(other) > 0:
            return True
        if self.compare_tightness(other) < 0:
            return False
        elif self.compare_time(other) > 0:
            return True
        elif self.compare_time(other) < 0:
            return False
        elif self.compare_coeff(other) > 0:
            return True
        else:
            return False
    
    def get_edge_label(self, other):
        is_better_or_equal = True
        if self.compare_tightness(other) < 0:
            is_better_or_equal = False
            self.is_tightest = False
        if self.compare_time(other) < 0:
            is_better_or_equal = False
            self.is_fastest = False
        if self.compare_coeff(other) < 0:
            is_better_or_equal = False
            self.has_best_coeff = False
        label = ""
        if self.compare_tightness(other) > 0:
            other.is_tightest = False
            label = label + "d"
        if self.compare_time(other) > 0:
            other.is_fastest = False
            label = label + "t"
        if self.compare_coeff(other) > 0:
            other.has_best_coeff = False
            label = label + "c"
        if is_better_or_equal and (label != ""):
            other.is_maximum = False
            if self.name == other.name:
                other.discard = True
        if is_better_or_equal:
            return label
        else:
            return ""
    
    def as_csv_row(self):
        csv_row = self.name
        if self.config != "":
            csv_row = csv_row + "," + self.config
        csv_row = csv_row + ";" + str(self.lb) + ";" + str(self.ub) + ";" + str(self.t) + ";" + str(self.coeff_lb) + ";" + str(self.coeff_ub) + "\n"
        return csv_row
        

def parse_method_name(method_name):
    method_name_list = method_name.split(",", 1)
    if (len(method_name_list) == 1):
        method_name_list.append("")
    return method_name_list

def dfs(adj_list, is_discarded_edge, method_id, child, seen):
    if seen[child]:
        return
    for grandchild in adj_list[child]:
        is_discarded_edge[method_id][grandchild] = True;
        dfs(adj_list, is_discarded_edge, method_id, grandchild, seen)
    seen[child] = True

def build_dependency_graph(result_file_name, consider_lb, tikz_file_name):
    # read csv file
    methods = []
    with open(result_file_name, "r") as result_file:  
        csv_reader = csv.reader(result_file,delimiter=";")
        next(csv_reader, NONE)
        for row in csv_reader:
            methods.append(Method(consider_lb, parse_method_name(row[0]), row[1], row[2], row[3], row[4], row[5]))
    # sort the methods in lexicographical order
    methods.sort()
    # construct dominance graph
    num_methods = len(methods)
    adj_list = [[] for method_id in range(0,num_methods)]
    adj_list_labels = [[] for method_id in range(0,num_methods)]
    for id_1 in range(0, num_methods):
        method_1 = methods[id_1]
        for id_2 in range(0, num_methods):
            method_2 = methods[id_2]
            edge_label = method_1.get_edge_label(method_2)
            if edge_label != "":
                adj_list[id_1].append(id_2)
                adj_list_labels[id_1].append(edge_label)
    # discard methods that are dominated by themselves with a different configuration
    undiscarded_method_ids = []
    for id_1 in range(0, num_methods):
        if not methods[id_1].discard:
            undiscarded_method_ids.append(id_1)
    undiscarded_adj_list = {method_id : [] for method_id in undiscarded_method_ids}
    undiscarded_adj_list_labels = {method_id : [] for method_id in undiscarded_method_ids}
    for id_1 in undiscarded_method_ids:
        for index in range(0, len(adj_list[id_1])):
            id_2 = adj_list[id_1][index]
            if not methods[id_2].discard:
                undiscarded_adj_list[id_1].append(id_2)
                undiscarded_adj_list_labels[id_1].append(adj_list_labels[id_1][index])
    # compute transitive reduction of undiscarded methods
    is_discarded_edge = {id_1 : {id_2 : False for id_2 in undiscarded_method_ids} for id_1 in undiscarded_method_ids}
    for method_id in undiscarded_method_ids:
        seen = {child : False for child in undiscarded_method_ids}
        for child in undiscarded_adj_list[method_id]:
            dfs(undiscarded_adj_list, is_discarded_edge, method_id, child, seen)
    # collect edges contained in transitive reduction
    edges = []
    edge_labels = []
    for id_1 in undiscarded_method_ids:
        nbs = undiscarded_adj_list[id_1]
        for index in range(0, len(nbs)):
            id_2 = nbs[index]
            if not is_discarded_edge[id_1][id_2]:
                edges.append((id_1, id_2))
                edge_labels.append(undiscarded_adj_list_labels[id_1][index])
    # write csv file containing only maxima
    maxima_file_name = result_file_name.split(".")[0]
    if consider_lb:
        maxima_file_name = maxima_file_name + "__LB_maxima.csv"
    else:
        maxima_file_name = maxima_file_name + "__UB_maxima.csv"
    maxima_file = open(maxima_file_name, "w")
    maxima_file.write("method;avg_lb;avg_ub;avg_runtime;classification_coefficient_lb;classification_coefficient_ub\n")
    for method in methods:
        if method.is_maximum:
            maxima_file.write(method.as_csv_row())
    maxima_file.close()
    # construct tikz file
    tikz_file = open(tikz_file_name, "w")
    tikz_file.write("%!TEX root = ../root.tex\n")
    tikz_file.write("\input{img/tikzstyles}\n")
    tikz_file.write("\\begin{tikzpicture}[rounded corners]\n")
    tikz_file.write("\graph[layered layout,\n")
    tikz_file.write("level distance=15mm,\n")
    tikz_file.write("nodes={minimum width=4mm, minimum height=4mm, align=center, font=\scriptsize},\n")
    tikz_file.write("edge quotes mid,\n")
    tikz_file.write("edges={nodes={font=\scriptsize, fill=white, inner sep=1.5pt}}] {\n")
    edge_id = 0
    tikz_file.write("{ [same layer]\n")
    for method_id in undiscarded_method_ids:
        method = methods[method_id]
        if method.is_maximum:
            tikz_file.write(str(method_id) + " [" + method.name + ", as={" + method.tikz_descriptor() + "}, label=above:{\\footnotesize " + method.label() + "}],\n");
    tikz_file.write("};\n")
    for method_id in undiscarded_method_ids:
        method = methods[method_id]
        if not method.is_maximum:
            tikz_file.write(str(method_id) + " [" + method.name + ", as={" + method.tikz_descriptor() + "}];\n");
    for edge in edges:
        tikz_file.write(str(edge[0]) + " -> [" + edge_labels[edge_id] + "] " + str(edge[1]) + ";\n")
    tikz_file.write("};\n")
    tikz_file.write("\end{tikzpicture}")
    tikz_file.close()
    
    
parser = argparse.ArgumentParser(description="Generates TikZ dominance graph from CSV file.")
parser.add_argument("result_file_name", help="path to existing CSV file")
parser.add_argument("tikz_file_name", help="path to existing CSV file")
parser.add_argument("--lb", help="generate dominance graph for lower bounds", action="store_true")
args = parser.parse_args()
build_dependency_graph(args.result_file_name, args.lb, args.tikz_file_name)