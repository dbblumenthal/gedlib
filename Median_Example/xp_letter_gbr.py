import pdb
import sys
import pathlib
import numpy as np
import networkx as nx

sys.path.insert(0, "../")

from ctypes import *
lib1 = cdll.LoadLibrary('../lib/fann/libdoublefann.so')
lib2 = cdll.LoadLibrary('../lib/libsvm.3.22/libsvm.so')
lib3 = cdll.LoadLibrary('../lib/nomad/libnomad.so')
lib4 = cdll.LoadLibrary('../lib/nomad/libsgtelib.so')

import gedlibpy

sys.path.insert(0, "optim-graphes/")
import pygraph

from median import draw_Letter_graph, compute_median, compute_median_set

gedlibpy.load_GXL_graphs('../include/gedlib-master/data/datasets/Letter/HIGH/', '../include/gedlib-master/data/collections/Letter.xml')
gedlibpy.set_edit_cost("LETTER")
gedlibpy.init()
gedlibpy.set_method("IPFP", "")
gedlibpy.init_method()

dataset,my_y = pygraph.utils.graphfiles.loadDataset("../include/gedlib-master/data/datasets/Letter/HIGH/Letter_all.cxl")
print("toto")

letters = np.unique(my_y)

listID = gedlibpy.get_all_graph_ids()
medians = {}
set_medians = {}
sods = {}
for f in letters:
    print(f)
    sublist_id = [listID[i] for i in range(0,len(listID)) if my_y[i] == f]
    sub_dataset = [dataset[i] for i  in range(0,len(listID)) if my_y[i] == f]
    median, sod, sods_path,set_median = compute_median(gedlibpy,sublist_id,sub_dataset,verbose=True)
    medians[f] = median
    set_medians[f] = set_median
    sods[f] = sods_path
    #draw_Letter_graph(median)  
  
import pickle

# write a file
f = open("medians", "wb")
pickle.dump(medians, f)
pickle.dump(set_medians, f)
pickle.dump(sods, f)
pickle.dump(letters, f)
f.close()
