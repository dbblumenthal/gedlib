from random import randint
from random import shuffle
from os.path import join

class Tree:
    
    def __init__(self, num_nodes, edge_list):
        self.num_nodes = num_nodes
        self.node_labels = [0 for node in range(num_nodes)]
        self.edge_list = edge_list
    
    def __repr__(self):
        string = "num_nodes =  " + str(self.num_nodes) + "\n"
        string = string + "node_labels =  " + str(self.node_labels) + "\n"
        string = string + "edge_list =  " + str(self.edge_list)
        return string
    
    def generate_node_labels(self, num_labels):
        for node in range(self.num_nodes):
            self.node_labels[node] = randint(0, num_labels - 1)
            
    def write_to_gxl(self, directory, file_name):
        gxl_file_name = join(directory, file_name)
        gxl_file = open(gxl_file_name, "w")
        gxl_file.write("<?xml version=\"1.0\"?>\n")
        gxl_file.write("<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n")
        gxl_file.write("<gxl>\n")
        gxl_file.write("<graph id=\"" + file_name + "\" edgeids=\"false\" edgemode=\"undirected\">\n")
        for node in range(self.num_nodes):
            gxl_file.write("<node id=\"_" + str(node) + "\">\n")
            gxl_file.write("<attr name=\"chem\"><int>" + str(self.node_labels[node]) + "</int></attr>\n")
            gxl_file.write("</node>\n")
        for edge in self.edge_list:
            gxl_file.write("<edge from=\"_" + str(edge[0]) + "\" to=\"_" + str(edge[1]) + "\">\n")
            gxl_file.write("<attr name=\"valence\"><int>1</int></attr>\n")
            gxl_file.write("</edge>\n")
        gxl_file.write("</graph>\n")
        gxl_file.write("</gxl>\n")
        gxl_file.close()
        

def generate_canonical_pruefer_seq(num_nodes):
    # generate Pruefer sequence
    pruefer_sec = []
    for i in range(num_nodes - 2):
        pruefer_sec.append(randint(0, num_nodes - 1))
    # compute canonical form
    old_to_new_id = {}
    new_id = 0
    for i in range(num_nodes - 2):
        old_id = pruefer_sec[i]
        if not old_id in old_to_new_id:
            old_to_new_id[old_id] = new_id
            new_id = new_id + 1
        pruefer_sec[i] = old_to_new_id[old_id]
    # return canonical form
    return pruefer_sec
    
def pruefer_seq_to_tree(pruefer_seq):
    # collect the degrees
    num_nodes = len(pruefer_seq) + 2
    deg = [1 for node in range(num_nodes)]
    for node in pruefer_seq:
        deg[node] = deg[node] + 1
    # collect all edges except for the last
    edge_list = []
    for tail in pruefer_seq:
        for head in range(num_nodes):
            if deg[head] == 1:
                edge_list.append((tail, head))
                deg[tail] = deg[tail] - 1
                deg[head] = deg[head] - 1
                break
    # collect last edge
    u = num_nodes
    v = num_nodes
    for node in range(num_nodes):
        if deg[node] == 1:
            if u == num_nodes:
                u = node
            else:
                v = node
                break
    edge_list.append((u, v))
    # return tree
    return Tree(num_nodes, edge_list)
    
def generate_molecules(num_molecules, min_num_nodes, max_num_nodes, max_num_trials = 10):
    # create non-isomorphic Pruefer sequences
    seqs = []
    for i in range(num_molecules):
        found_new_pruefer_seq = False
        num_trials = 0
        while (not found_new_pruefer_seq) and (num_trials < max_num_trials):
            num_nodes = randint(min_num_nodes, max_num_nodes)
            new_pruefer_seq = generate_canonical_pruefer_seq(num_nodes)
            found_new_pruefer_seq = True
            for old_pruefer_seq in seqs:
                if old_pruefer_seq == new_pruefer_seq:
                    found_new_pruefer_seq = False
                    break;
            if found_new_pruefer_seq:
                seqs.append(new_pruefer_seq)
            else:
                num_trials = num_trials + 1
        if num_trials == max_num_trials:
            raise Exception("Cannot generate new Pruefer sequence. Maximum number of trials reached.")
    # shuffle the Pruefer sequences
    for pruefer_seq in seqs:
        num_nodes = len(pruefer_seq) + 2
        permutation = range(num_nodes)
        shuffle(permutation)
        for i in range(num_nodes - 2):
            pruefer_seq[i] = permutation[pruefer_seq[i]]
    # return trees obtained from Pruefer sequences
    return [pruefer_seq_to_tree(pruefer_seq) for pruefer_seq in seqs]
    
molecules = generate_molecules(150, 8, 12)
file_names = ["molecule_" + str(i) + ".gxl" for i in range(150)]
num_labels = [1, 4, 7, 10]
dirs = ["datasets/S-MOL/NL01", "datasets/S-MOL/NL04", "datasets/S-MOL/NL07", "datasets/S-MOL/NL10"]
for dir_id in range(4):
    for mol_id in range(150):
        molecules[mol_id].generate_node_labels(num_labels[dir_id])
        molecules[mol_id].write_to_gxl(dirs[dir_id], file_names[mol_id])
collection = open("collections/S-MOL.xml", "w")
collection.write("<?xml version=\"1.0\"?>\n")
collection.write("<!DOCTYPE GraphCollection SYSTEM \"http://www.inf.unibz.it/~blumenthal/dtd/GraphCollection.dtd\">\n")
collection.write("<GraphCollection>\n")
for file_name in file_names:
    collection.write("<graph file=\"" + file_name + "\" class=\"a\"/>\n")
collection.write("</GraphCollection>\n")
collection.close()