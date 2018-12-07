import csv
from pickle import NONE
import argparse

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

def build_dependency_graph(result_file_name, consider_lb):
    methods = []
    dists = []
    times = []
    coeffs = []
    with open(result_file_name, "r") as result_file:  
        csv_reader = csv.reader(result_file,delimiter=";")
        next(csv_reader, NONE)
        for row in csv_reader:
            methods.append(parse_method_name(row[0]))
            times.append(float(row[3]))
            if consider_lb:
                dists.append(float(row[1]))
                coeffs.append(float(row[4]))
            else:
                dists.append(float(row[2]))
                coeffs.append(float(row[5]))
    num_methods = len(methods)
    adj_list = [[] for method_id in range(0,num_methods)]
    adj_list_types = [[] for method_id in range(0,num_methods)]
    is_maximum = [True for method_id in range(0,num_methods)]
    for id_1 in range(0, num_methods):
        time_1 = times[id_1]
        dist_1 = dists[id_1]
        coeff_1 = coeffs[id_1]
        for id_2 in range(0, num_methods):
            time_2 = times[id_2]
            dist_2 = dists[id_2]
            coeff_2 = coeffs[id_2]
            method_1_better_equal_2 = (time_1 <= time_2) and (coeff_1 >= coeff_2)
            if (consider_lb):
                method_1_better_equal_2 = method_1_better_equal_2 and (dist_1 >= dist_2)
            else:
                method_1_better_equal_2 = method_1_better_equal_2 and (dist_1 <= dist_2)
            edge_type = ""
            if (consider_lb):
                if (dist_1 > dist_2):
                    edge_type = edge_type + "$d$"
            else:
                if (dist_1 < dist_2):
                    edge_type = edge_type + "$d$"
            if (time_1 < time_2):
                if (edge_type == ""):
                    edge_type = edge_type + "$t$"
                else:
                    edge_type = edge_type + ",$t$"
            if (coeff_1 > coeff_2):
                if (edge_type == ""):
                    edge_type = edge_type + "$c$"
                else:
                    edge_type = edge_type + ",$c$"
            if method_1_better_equal_2 and (edge_type != ""):
                adj_list[id_1].append(id_2)
                adj_list_types[id_1].append(edge_type)
                is_maximum[id_2] = False
    is_discard_edge = [[False for nb in range(0, num_methods)] for method_id in range(0,num_methods)]
    for method_id in range(0, num_methods):
        seen = [False for child in range(0, num_methods)]
        for child in adj_list[method_id]:
            dfs(adj_list, is_discard_edge, method_id, child, seen)
    edge_types = []
    edges = []
    for method_id in range(0, num_methods):
        nbs = adj_list[method_id]
        nbs_types = adj_list_types[method_id]
        for index in range(0, len(nbs)):
            if not is_discard_edge[method_id][nbs[index]]:
                edges.append((method_id, nbs[index]))
                edge_types.append(nbs_types[index])
    return methods, edges, edge_types, is_maximum

def write_dependency_graph(tikz_file_name, methods, edges, edge_types, is_maximum):
    tikz_file = open(tikz_file_name, "w")
    tikz_file.write("%!TEX root = ../root.tex\n")
    tikz_file.write("\input{img/tikzstyles}\n")
    tikz_file.write("\\begin{tikzpicture}[rounded corners]\n")
    tikz_file.write("\graph[layered layout,\n")
    tikz_file.write("nodes={minimum width=4mm, minimum height=4mm, font=\scriptsize},\n")
    tikz_file.write("edge quotes mid,\n")
    tikz_file.write("edges={nodes={font=\scriptsize, fill=white, inner sep=1.5pt}}] {\n")
    edge_id = 0
    tikz_file.write("{ [same layer]\n")
    for method_id in range(0, len(methods)):
        if is_maximum[method_id]:
            tikz_file.write(str(method_id) + " [" + methods[method_id][0] + ", as={" + methods[method_id][1] + "}],\n");
    tikz_file.write("};\n")
    for method_id in range(0, len(methods)):
        if not is_maximum[method_id]:
            tikz_file.write(str(method_id) + " [" + methods[method_id][0] + ", as={" + methods[method_id][1] + "}];\n");
    for edge in edges:
        tikz_file.write(str(edge[0]) + " -> [\"" + edge_types[edge_id] + "\"] " + str(edge[1]) + ";\n")
    tikz_file.write("};\n")
    tikz_file.write("\end{tikzpicture}")
    tikz_file.close()
    
parser = argparse.ArgumentParser(description="Generates TikZ dominance graph from CSV file.")
parser.add_argument("result_file_name", help="path to existing CSV file")
parser.add_argument("tikz_file_name", help="path to existing CSV file")
parser.add_argument("--lb", help="generate dominance graph for lower bounds", action="store_true")
args = parser.parse_args()
methods, edges, edge_types, is_maximum = build_dependency_graph(args.result_file_name, args.lb)
write_dependency_graph(args.tikz_file_name, methods, edges, edge_types, is_maximum)