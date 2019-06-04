import xml.etree.ElementTree as ET
import argparse
import os.path
import planarity
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import numpy as np

def is_planar(edge_list):
	return planarity.is_planar(edge_list)

def num_connected_components(adj_matrix):
	graph = csr_matrix(adj_matrix)
	num_components = connected_components(csgraph=graph, directed=False, return_labels=False)
	return num_components
	
def is_cyclic_util(adj_list, node_id, visited, parent_id):
	visited[node_id] = True
	for neighbor_id in adj_list[node_id]:
		if not visited[neighbor_id]:
			if is_cyclic_util(adj_list, neighbor_id, visited, node_id):
				return True
		elif neighbor_id != parent_id:
				return True
	return False
			
def is_cyclic(adj_list, num_nodes):
	visited = [False for node_id in range(num_nodes)]
	for node_id in range(num_nodes):
		if not visited[node_id]:
			if is_cyclic_util(adj_list, node_id, visited, -1):
				return True
	return False
    
def as_2d_hist(min_num_nodes, max_num_nodes, min_num_edges, max_num_edges, nums_nodes, nums_edges):
    hist = {num_nodes : {num_edges : 0 for num_edges in range(min_num_edges, max_num_edges + 1)} for num_nodes in range(min_num_nodes, max_num_nodes + 1)}
    num_graphs = len(nums_nodes)
    for i in range(num_graphs):
        hist[nums_nodes[i]][nums_edges[i]] = hist[nums_nodes[i]][nums_edges[i]] + 1
    return hist

def parse_graph(dir, gxl_file):
    num_nodes = 0
    graph = ET.parse(os.path.join(dir, gxl_file)).getroot()
    str_to_id = {}
    for node in graph.findall("graph/node"):
        str_to_id.update({node.attrib["id"] : num_nodes})	
        num_nodes = num_nodes + 1
    edge_list = []
    adj_list = {node_id : [] for node_id in range(num_nodes)}
    adj_matrix = [[0 for col in range(num_nodes)] for row in range(num_nodes)]
    for edge in graph.findall("graph/edge"):
		tail = str_to_id[edge.attrib["from"]]
		head = str_to_id[edge.attrib["to"]]
		if adj_matrix[tail][head] == 1:
			continue
		adj_matrix[tail][head] = 1
		adj_matrix[head][tail] = 1
		edge_list.append((tail,head))
		adj_list[tail].append(head)
		adj_list[head].append(tail)
    return edge_list, adj_list, adj_matrix, num_nodes, len(edge_list)

# Parse the command line arguments.
parser = argparse.ArgumentParser(description="Computes dataset statistics.")
parser.add_argument("dataset", help="path to existing dataset file")
parser.add_argument("dir", help="path to directory containing GXL files")
parser.add_argument("--max_size", help="maximal size", type=int)
parser.add_argument("--topology", help="compute topology avgerage number components and ratio of acyclic and planar graphs", action="store_true")
parser.add_argument("--distr", help="file to save distibution of number of nodes and edges")
args = parser.parse_args()
if not os.path.isdir(args.dir):
        raise Exception("Invalid argument \"" + dir + "\": not a directory. Usage: python analyze_dataset.py <dataset> <dir> [--max-size <arg>]")

# Parse the dataset file.
dataset = ET.parse(args.dataset).getroot()
graphs = [graph.attrib["file"] for graph in dataset]

# Analyze the dataset.
num_graphs = 0
total_num_graphs = 0
ratio_acyclic = 0.0
ratio_planar = 0.0
nums_nodes = np.array([], dtype=int)
nums_edges = np.array([], dtype=int)
nums_components = np.array([], dtype=int)
for gxl_file in graphs:
	total_num_graphs = total_num_graphs + 1
	edge_list, adj_list, adj_matrix, num_nodes, num_edges = parse_graph(args.dir, gxl_file)
	if args.max_size and num_nodes > args.max_size:
		continue
	num_graphs = num_graphs + 1
	nums_nodes = np.append(nums_nodes, [num_nodes])
	nums_edges = np.append(nums_edges, [num_edges])
	if args.topology:
	    num_components = num_nodes
	    if num_edges > 0:
		    num_components = num_connected_components(adj_matrix)
	    nums_components = np.append(nums_components, [components])
	    if num_edges == 0 or not is_cyclic(adj_list, num_nodes):
		    ratio_acyclic = ratio_acyclic + 1.0
	    if num_edges == 0 or is_planar(edge_list):
		    ratio_planar = ratio_planar + 1.0
ratio_acyclic = ratio_acyclic / float(num_graphs)
ratio_planar = ratio_planar / float(num_graphs)
print("graphs (total / not filtered): " + str(total_num_graphs) + " / " + str(num_graphs))
print("nodes (min / max / mean / std / median): " + str(np.min(nums_nodes)) + " & " + str(np.max(nums_nodes)) + " & " + str(np.mean(nums_nodes)) + " & " + str(np.std(nums_nodes, ddof=1)) + " & " + str(np.median(nums_nodes)))
print("edges (min / max / mean / std / median): " + str(np.min(nums_edges)) + " & " + str(np.max(nums_edges)) + " & " + str(np.mean(nums_edges)) + " & " + str(np.std(nums_edges, ddof=1)) + " & " + str(np.median(nums_edges)))
if args.topology:
    print("components (min / max / mean / std / median): " + str(np.min(nums_components)) + " & " + str(np.max(nums_components)) + " & " + str(np.mean(nums_components)) + " & " + str(np.std(nums_components, ddof=1)) + " & " + str(np.median(nums_components)))
    print("ratio planar graphs: " + str(ratio_planar))
    print("ratio acylcic graphs: " + str(ratio_acyclic))

if args.distr:
    hist = as_2d_hist(min_num_nodes, max_num_nodes, min_num_edges, max_num_edges, nums_nodes, nums_edges)
    file = open(args.distr + ".dat", "w")
    for num_nodes in range(min_num_nodes, max_num_nodes + 1):
        for num_edges in range(min_num_edges, max_num_edges + 1):
            file.write(str(num_nodes) + " " + str(num_edges) + " " + str(hist[num_nodes][num_edges]) + "\n")
        file.write("\n")
    file.close()
	