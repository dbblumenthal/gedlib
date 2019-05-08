import csv
from pickle import NONE
import argparse
from decimal import Decimal
import os.path
from numpy import Infinity
from buildtools import process
    
class Method:
    
    def __init__(self, name, ub, t, ub_diag):
        self.consider_lb = True
        self.name = name[0]
        if self.name == "KREFINE":
            self.name = "3REFINE"
        elif self.name == "REFINE":
            self.name = "2REFINE"
        self.num_randpost_loops = name[1][2:-2].split(",")[2]
        self.randpost_penalty = name[1][2:-2].split(",")[3]
        self.ub = float("{0:.8f}".format(Decimal(ub)))
        self.t = float("{:.8f}".format(Decimal(t)))
        self.ub_diag = float("{0:.8f}".format(Decimal(ub_diag)))

def parse_method_name(method_name):
    method_name_list = method_name.split(",", 1)
    if (len(method_name_list) == 1):
        method_name_list.append("")
    return method_name_list
    
def read_results(result_file_name):
    methods = []
    with open(result_file_name, "r") as result_file:  
        csv_reader = csv.reader(result_file, delimiter=";")
        next(csv_reader, NONE)
        for row in csv_reader:
            methods.append(Method(parse_method_name(row[0]), row[2], row[3], row[4]))
    return methods

def find_best_randpost_penalties(methods):
    method_names = set()
    randpost_penalties = set()
    for method in methods:
        method_names.add(method.name)
        randpost_penalties.add(method.randpost_penalty)
    ub_sums = {method_name : { randpost_penalty : 0.0 for randpost_penalty in randpost_penalties } for method_name in method_names}
    for method in methods:
        if method.num_randpost_loops != "0":
            ub_sums[method.name][method.randpost_penalty] = ub_sums[method.name][method.randpost_penalty] + method.ub
    best_randpost_penalties = {method_name : "" for method_name in method_names}
    best_ub_sums = {method_name : Infinity for method_name in method_names}
    for method_name, method_dict in ub_sums.items():
        for randpost_penalty, ub_sum in method_dict.items():
            if ub_sum < best_ub_sums[method_name]:
                best_ub_sums[method_name] = ub_sum
                best_randpost_penalties[method_name] = randpost_penalty
    return best_randpost_penalties

def write_best_randpost_penalties(randpost_penalties_file_name, best_randpost_penalties):
    randpost_penalties_file = open(randpost_penalties_file_name, "w")
    col = 1
    num_cols = len(best_randpost_penalties)
    for method_name in best_randpost_penalties:
        randpost_penalties_file.write(method_name)
        if col < num_cols:
            randpost_penalties_file.write(",")
        else:
            randpost_penalties_file.write("\n")
        col = col + 1
    col = 1
    for method_name, randpost_penalty in best_randpost_penalties.items():
        randpost_penalties_file.write(randpost_penalty)
        if col < num_cols:
            randpost_penalties_file.write(",")
        else:
            randpost_penalties_file.write("\n")
        col = col + 1
    randpost_penalties_file.close()

def write_pgf_file(pgf_file_name, best_randpost_penalties, methods):
    method_ids = {num_randpost_loops : {method_name : -1 for method_name in best_randpost_penalties} for num_randpost_loops in ["0", "1", "3", "7"]}
    method_id = 0
    for method in methods:
        if (method.num_randpost_loops == "0") or (method.randpost_penalty == best_randpost_penalties[method.name]):
            method_ids[method.num_randpost_loops][method.name] = method_id
        method_id = method_id + 1
    pgf_file = open(pgf_file_name, "w")
    pgf_file.write("L")
    for method_name in best_randpost_penalties:
        pgf_file.write("," + method_name + "_avg_ub," + method_name + "_avg_runtime," + method_name + "_avg_ub_diag")
    pgf_file.write("\n")
    for num_randpost_loops in ["0", "1", "3", "7"]:
        pgf_file.write(num_randpost_loops)
        for method_name in best_randpost_penalties:
            method = methods[method_ids[num_randpost_loops][method_name]]
            pgf_file.write("," + str(method.ub) + "," + str(method.t) + "," + str(method.ub_diag))
        pgf_file.write("\n")
    pgf_file.close()
    
parser = argparse.ArgumentParser(description="Determines best RANDPOST penalties and outputs .")
parser.add_argument("result_file_name", help="name of input result file")
parser.add_argument("pgf_file_name", help="name of output file to be used for PGF plots")
parser.add_argument("randpost_penalties_file_name", help="name of output file listing best RANDPOST penalties")
args = parser.parse_args()
methods = read_results(args.result_file_name)
best_randpost_penalties = find_best_randpost_penalties(methods)
write_best_randpost_penalties(args.randpost_penalties_file_name, best_randpost_penalties)
write_pgf_file(args.pgf_file_name, best_randpost_penalties, methods)