import xml.etree.ElementTree as ET
import argparse
import os.path

def graph_size(dir, gxl_file):
    graph = ET.parse(os.path.join(dir, gxl_file)).getroot()
    num_nodes = 0
    for node in graph.findall("graph/node"):
        num_nodes = num_nodes + 1
    num_edges = 0
    for node in graph.findall("graph/edge"):
        num_edges = num_edges + 1
    return num_nodes, num_edges

# Parse the command line arguments.
parser = argparse.ArgumentParser(description="Creates collection of graphs with given number of nodes from directory containing GXL files.")
parser.add_argument("dataset", help="path to existing dataset file")
parser.add_argument("dir", help="path to directory containing GXL files")
args = parser.parse_args()
if not os.path.isdir(args.dir):
        raise Exception("Invalid argument \"" + dir + "\": not a directory. Usage: python graphs_of_given_size.py <dataset> <collection> <dir> <min-num-nodes> <max-num-nodes>")

# Parse the dataset file.
dataset = ET.parse(args.dataset).getroot()
graphs = [graph.attrib["file"] for graph in dataset]

# Compute the statistics.
max_num_nodes = 0
avg_num_nodes = 0.0
max_num_edges = 0
avg_num_edges = 0.0
for graph in graphs:
    num_nodes, num_edges = graph_size(args.dir, graph)
    max_num_nodes = max(max_num_nodes, num_nodes)
    avg_num_nodes = avg_num_nodes + num_nodes
    max_num_edges = max(max_num_edges, num_edges)
    avg_num_edges = avg_num_edges + num_edges
avg_num_nodes = avg_num_nodes / float(len(graphs))
avg_num_edges = avg_num_edges / float(len(graphs))

# Print the statistics.
print("max. num. nodes = %d, avg. num. nodes = %f" % (max_num_nodes, avg_num_nodes))
print("max. num. edges = %d, avg. num. edges = %f" % (max_num_edges, avg_num_edges))