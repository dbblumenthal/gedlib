import xml.etree.ElementTree as ET
import argparse
import os.path
import planarity
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

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
parser = argparse.ArgumentParser(description="Creates collection of graphs with given number of nodes from directory containing GXL files.")
parser.add_argument("dataset", help="path to existing dataset file")
parser.add_argument("dir", help="path to directory containing GXL files")
parser.add_argument("--max_size", help="maximal size", type=int)
args = parser.parse_args()
if not os.path.isdir(args.dir):
        raise Exception("Invalid argument \"" + dir + "\": not a directory. Usage: python graphs_of_given_size.py <dataset> <collection> <dir> <min-num-nodes> <max-num-nodes>")

# Parse the dataset file.
dataset = ET.parse(args.dataset).getroot()
graphs = [graph.attrib["file"] for graph in dataset]

# Analyze the dataset.
num_graphs = 0.0
total_num_graphs = 0
avg_num_nodes = 0.0
max_num_nodes = 0
min_num_nodes = 100000000
avg_num_edges = 0.0
max_num_edges = 0
min_num_edges = 100000000
avg_num_components = 0.0
max_num_components = 0
min_num_components = 100000000
ratio_acyclic = 0.0
ratio_planar = 0.0
for gxl_file in graphs:
	total_num_graphs = total_num_graphs + 1
	edge_list, adj_list, adj_matrix, num_nodes, num_edges = parse_graph(args.dir, gxl_file)
	if args.max_size and num_nodes > args.max_size:
		continue
	num_graphs = num_graphs + 1.0
	avg_num_nodes = avg_num_nodes + num_nodes
	max_num_nodes = max(max_num_nodes, num_nodes)
	min_num_nodes = min(min_num_nodes, num_nodes)
	avg_num_edges = avg_num_edges + num_edges
	max_num_edges = max(max_num_edges, num_edges)
	min_num_edges = min(min_num_edges, num_edges)
	num_components = num_nodes
	if num_edges > 0:
		num_components = num_connected_components(adj_matrix)
	avg_num_components = avg_num_components + num_components
	max_num_components = max(max_num_components, num_components)
	min_num_components = min(min_num_components, num_components)
	if num_edges == 0 or not is_cyclic(adj_list, num_nodes):
		ratio_acyclic = ratio_acyclic + 1.0
	if num_edges == 0 or is_planar(edge_list):
		ratio_planar = ratio_planar + 1.0
avg_num_nodes = avg_num_nodes / num_graphs
avg_num_edges = avg_num_edges / num_graphs
avg_num_components = avg_num_components / num_graphs
ratio_acyclic = ratio_acyclic / num_graphs
ratio_planar = ratio_planar / num_graphs
print("graphs (total / not filtered): " + str(total_num_graphs) + " / " + str(num_graphs))
print("nodes (avg / max / min): " + str(avg_num_nodes) + " / " + str(max_num_nodes) + " / " + str(min_num_nodes))
print("edges (avg / max / min): " + str(avg_num_edges) + " / " + str(max_num_edges) + " / " + str(min_num_edges))
print("components (avg / max / min): " + str(avg_num_components) + " / " + str(max_num_components) + " / " + str(min_num_components))
print("ratio planar graphs: " + str(ratio_planar))
print("ratio acylcic graphs: " + str(ratio_acyclic))
	