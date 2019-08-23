# distutils: language = c++

"""
    Python GedLib module
    ======================
    
    This module allow to use a C++ library for edit distance between graphs (GedLib) with Python.

    
    Authors
    -------------------
 
    David Blumenthal
    Natacha Lambert

    Copyright (C) 2019 by all the authors

    Classes & Functions
    -------------------
 
"""

################################
##DEFINITIONS OF C++ FUNCTIONS##
################################


#Types imports for C++ compatibility
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.list cimport list

#Long unsigned int equivalent
cimport numpy as np
ctypedef np.npy_uint32 UINT32_t
from cpython cimport array

#Functions importation
cdef extern from "src/GedLibBind.h" :
    cdef vector[string] getEditCostStringOptions()
    cdef vector[string] getMethodStringOptions()
    cdef vector[string] getInitStringOptions()
    cdef bool isInitialized()
    cdef void restartEnv()
    cdef void loadGXLGraph(string pathFolder, string pathXML)
    cdef pair[size_t,size_t] getGraphIds()
    cdef vector[size_t] getAllGraphIds()
    cdef string getGraphClass(size_t id)
    cdef string getGraphName(size_t id)
    cdef size_t addGraph(string name, string classe)
    cdef void addNode(size_t graphId, string nodeId, map[string,string] nodeLabel)
    cdef void addEdge(size_t graphId, string tail, string head, map[string,string] edgeLabel, bool ignoreDuplicates)
    cdef void clearGraph(size_t graphId)
    cdef size_t getGraphInternalId(size_t graphId)
    cdef size_t getGraphNumNodes(size_t graphId)
    cdef size_t getGraphNumEdges(size_t graphId)
    cdef vector[string] getGraphOriginalNodeIds(size_t graphId)
    cdef vector[map[string, string]] getGraphNodeLabels(size_t graphId)
    cdef map[pair[size_t,size_t], map[string,string]] getGraphEdges(size_t graphId)
    cdef vector[vector[size_t]] getGraphAdjacenceMatrix(size_t graphId)
    cdef void setEditCost(string editCost, vector[double] editCostConstant)
    cdef void setPersonalEditCost(vector[double] editCostConstant)
    cdef void initEnv(string initOption)
    cdef void setMethod(string method, string options)
    cdef void initMethod()
    cdef double getInitime()
    cdef void runMethod(size_t g, size_t h)
    cdef double getUpperBound(size_t g, size_t h)
    cdef double getLowerBound(size_t g, size_t h)
    cdef vector[np.npy_uint64] getForwardMap(size_t g, size_t h)
    cdef vector[np.npy_uint64] getBackwardMap(size_t g, size_t h)
    cdef size_t getNodeImage(size_t g, size_t h, size_t nodeId)
    cdef size_t getNodePreImage(size_t g, size_t h, size_t nodeId)
    cdef size_t getDummyNode()
    cdef vector[pair[size_t,size_t]] getNodeMap(size_t g, size_t h)
    cdef vector[vector[int]] getAssignmentMatrix(size_t g, size_t h)
    cdef vector[vector[np.npy_uint64]] getAllMap(size_t g, size_t h)
    cdef double getRuntime(size_t g, size_t h)
    cdef bool quasimetricCosts()
    cdef vector[vector[size_t]] hungarianLSAP(vector[vector[size_t]] matrixCost);
    cdef vector[vector[double]] hungarianLSAPE(vector[vector[double]] matrixCost);


    
###########################################
##REDEFINITION OF C++ FUNCTIONS IN PYTHON##
###########################################

def is_initialized() :
    """
        Checks and returns if the computation environment is initialized or not.
 
        :return: True if it's initialized, False otherwise
        :rtype: bool
        
        .. note:: This function exists for internals verifications but you can use it for your code. 
    """
    return isInitialized()

def get_edit_cost_options() :
    """
        Searchs the differents edit cost functions and returns the result.
 
        :return: The list of edit cost functions
        :rtype: list[string]
 
        .. warning:: This function is useless for an external use. Please use directly list_of_edit_cost_options. 
        .. note:: Prefer the list_of_edit_cost_options attribute of this module.
    """
    
    return getEditCostStringOptions()

def get_method_options() :
    """
        Searchs the differents method for edit distance computation between graphs and returns the result.
 
        :return: The list of method to compute the edit distance between graphs
        :rtype: list[string]
 
        .. warning:: This function is useless for an external use. Please use directly list_of_method_options.
        .. note:: Prefer the list_of_method_options attribute of this module.
    """
    return getMethodStringOptions()

def get_init_options() :
    """
        Searchs the differents initialization parameters for the environment computation for graphs and returns the result.
 
        :return: The list of options to initialize the computation environment
        :rtype: list[string]
 
        .. warning:: This function is useless for an external use. Please use directly list_of_init_options.
        .. note:: Prefer the list_of_init_options attribute of this module.
    """
    return getInitStringOptions()

def restart_env() :
    """
        Restarts the environment variable. All data related to it will be delete. 
 
        .. warning:: This function deletes all graphs, computations and more so make sure you don't need anymore your environment. 
        .. note:: You can now delete and add somes graphs after initialization so you can avoid this function. 
    """
    restartEnv()

def load_GXL_graphs(path_folder, path_XML) :
    """
        Loads some GXL graphes on the environment which is in a same folder, and present in the XMLfile. 
        
        :param path_folder: The folder's path which contains GXL graphs
        :param path_XML: The XML's path which indicates which graphes you want to load
        :type path_folder: string
        :type path_XML: string
 
        .. note:: You can call this function multiple times if you want, but not after an init call. 
    """
    loadGXLGraph(path_folder.encode('utf-8'), path_XML.encode('utf-8'))

def graph_ids() :
    """
        Searchs the first and last IDs of the loaded graphs in the environment. 
 
        :return: The pair of the first and the last graphs Ids
        :rtype: tuple(size_t, size_t)
        
        .. note:: Prefer this function if you have huges structures with lots of graphs.  
    """
    return getGraphIds()

def get_all_graph_ids() :
    """
        Searchs all the IDs of the loaded graphs in the environment. 
 
        :return: The list of all graphs's Ids 
        :rtype: list[size_t]
        
        .. note:: The last ID is equal to (number of graphs - 1). The order correspond to the loading order. 
    """
    return getAllGraphIds()

def get_graph_class(id) :
    """
        Returns the class of a graph with its ID.

        :param id: The ID of the wanted graph
        :type id: size_t
        :return: The class of the graph which correpond to the ID
        :rtype: string
        
        .. seealso:: get_graph_class()
        .. note:: An empty string can be a class. 
    """
    return getGraphClass(id)

def get_graph_name(id) :
    """
        Returns the name of a graph with its ID. 

        :param id: The ID of the wanted graph
        :type id: size_t
        :return: The name of the graph which correpond to the ID
        :rtype: string
        
        .. seealso:: get_graph_class()
        .. note:: An empty string can be a name. 
    """
    return getGraphName(id)

def add_graph(name="", classe="") :
    """
        Adds a empty graph on the environment, with its name and its class. Nodes and edges will be add in a second time. 

        :param name: The name of the new graph, an empty string by default
        :param classe: The class of the new graph, an empty string by default
        :type name: string
        :type classe: string
        :return: The ID of the newly graphe
        :rtype: size_t
        
        .. seealso::add_node(), add_edge() , add_symmetrical_edge()
        .. note:: You can call this function without parameters. You can also use this function after initialization, call init() after you're finished your modifications. 
    """
    return addGraph(name.encode('utf-8'),classe.encode('utf-8'))

def add_node(graph_id, node_id, node_label):
    """
        Adds a node on a graph selected by its ID. A ID and a label for the node is required. 

        :param graph_id: The ID of the wanted graph
        :param node_id: The ID of the new node
        :param node_label: The label of the new node
        :type graph_id: size_t
        :type node_id: string
        :type node_label: dict{string : string}
        
        .. seealso:: add_graph(), add_edge(), add_symmetrical_edge()
        .. note:: You can also use this function after initialization, but only on a newly added graph. Call init() after you're finished your modifications. 
    """
    addNode(graph_id, node_id.encode('utf-8'), encode_your_map(node_label))

def add_edge(graph_id, tail, head, edge_label, ignore_duplicates = True) :
    """
        Adds an edge on a graph selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :param tail: The ID of the tail node for the new edge
        :param head: The ID of the head node for the new edge
        :param edge_label: The label of the new edge
        :param ignore_duplicates: If True, duplicate edges are ignored, otherwise it's raise an error if an existing edge is added. True by default
        :type graph_id: size_t
        :type tail: string
        :type head: string
        :type edge_label: dict{string : string}
        :type ignore_duplicates: bool
        
        .. seealso:: add_graph(), add_node(), add_symmetrical_edge()
        .. note:: You can also use this function after initialization, but only on a newly added graph. Call init() after you're finished your modifications. 
    """
    addEdge(graph_id, tail.encode('utf-8'), head.encode('utf-8'), encode_your_map(edge_label), ignore_duplicates)

def add_symmetrical_edge(graph_id, tail, head, edge_label) :
    """
        Adds a symmetrical edge on a graph selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :param tail: The ID of the tail node for the new edge
        :param head: The ID of the head node for the new edge
        :param edge_label: The label of the new edge
        :type graph_id: size_t
        :type tail: string
        :type head: string
        :type edge_label: dict{string : string}
        
        .. seealso:: add_graph(), add_node(), add_edge()
        .. note:: You can also use this function after initialization, but only on a newly added graph. Call init() after you're finished your modifications. 
    """
    tailB = tail.encode('utf-8')
    headB = head.encode('utf-8')
    edgeLabelB = encode_your_map(edge_label)
    addEdge(graph_id, tailB, headB, edgeLabelB, True)
    addEdge(graph_id, headB, tailB, edgeLabelB, True)

def clear_graph(graph_id) :
    """
        Deletes a graph, selected by its ID, to the environment.

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        
        .. note:: Call init() after you're finished your modifications. 
    """
    clearGraph(graph_id)

def get_graph_internal_id(graph_id) :
    """
        Searchs and returns the internal Id of a graph, selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        :return: The internal ID of the selected graph
        :rtype: size_t
        
        .. seealso:: get_graph_num_nodes(), get_graph_num_edges(), get_original_node_ids(), get_graph_node_labels(), get_graph_edges(), get_graph_adjacence_matrix()
        .. note:: These functions allow to collect all the graph's informations.
    """
    return getGraphInternalId(graph_id)

def get_graph_num_nodes(graph_id) :
    """
        Searchs and returns the number of nodes on a graph, selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        :return: The number of nodes on the selected graph
        :rtype: size_t
        
        .. seealso:: get_graph_internal_id(), get_graph_num_edges(), get_original_node_ids(), get_graph_node_labels(), get_graph_edges(), get_graph_adjacence_matrix()
        .. note:: These functions allow to collect all the graph's informations.
    """
    return getGraphNumNodes(graph_id)

def get_graph_num_edges(graph_id) :
    """
        Searchs and returns the number of edges on a graph, selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        :return: The number of edges on the selected graph
        :rtype: size_t
        
        .. seealso:: get_graph_internal_id(), get_graph_num_nodes(), get_original_node_ids(), get_graph_node_labels(), get_graph_edges(), get_graph_adjacence_matrix()
        .. note:: These functions allow to collect all the graph's informations.
    """
    return getGraphNumEdges(graph_id)

def get_original_node_ids(graph_id) :
    """
        Searchs and returns all th Ids of nodes on a graph, selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        :return: The list of IDs's nodes on the selected graph
        :rtype: list[string]
        
        .. seealso::get_graph_internal_id(), get_graph_num_nodes(), get_graph_num_edges(), get_graph_node_labels(), get_graph_edges(), get_graph_adjacence_matrix()
        .. note:: These functions allow to collect all the graph's informations.
    """
    return getGraphOriginalNodeIds(graph_id)

def get_graph_node_labels(graph_id) :
    """
        Searchs and returns all the labels of nodes on a graph, selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        :return: The list of labels's nodes on the selected graph
        :rtype: list[dict{string : string}]
        
        .. seealso:: get_graph_internal_id(), get_graph_num_nodes(), get_graph_num_edges(), get_original_node_ids(), get_graph_edges(), get_graph_adjacence_matrix()
        .. note:: These functions allow to collect all the graph's informations.
    """
    return getGraphNodeLabels(graph_id)

def get_graph_edges(graph_id) :
    """
        Searchs and returns all the edges on a graph, selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        :return: The list of edges on the selected graph
        :rtype: dict{tuple(size_t,size_t) : dict{string : string}}
        
        .. seealso::get_graph_internal_id(), get_graph_num_nodes(), get_graph_num_edges(), get_original_node_ids(), get_graph_node_labels(), get_graph_adjacence_matrix()
        .. note:: These functions allow to collect all the graph's informations.
    """
    return getGraphEdges(graph_id)

def get_graph_adjacence_matrix(graph_id) :
    """
        Searchs and returns the adjacence list of a graph, selected by its ID. 

        :param graph_id: The ID of the wanted graph
        :type graph_id: size_t
        :return: The adjacence list of the selected graph
        :rtype: list[list[size_t]]
        
        .. seealso:: get_graph_internal_id(), get_graph_num_nodes(), get_graph_num_edges(), get_original_node_ids(), get_graph_node_labels(), get_graph_edges()
        .. note:: These functions allow to collect all the graph's informations.
    """
    return getGraphAdjacenceMatrix(graph_id)

def set_edit_cost(edit_cost, edit_cost_constant = []) :
    """
        Sets an edit cost function to the environment, if its exists. 

        :param edit_cost: The name of the edit cost function
        :type edit_cost: string
        :param edi_cost_constant: The parameters you will add to the editCost, empty by default
        :type edit_cost_constant: list
        
        .. seealso:: list_of_edit_cost_options
        .. note:: Try to make sure the edit cost function exists with list_of_edit_cost_options, raise an error otherwise. 
    """
    editCostB = edit_cost.encode('utf-8')
    if editCostB in list_of_edit_cost_options : 
        setEditCost(editCostB, edit_cost_constant)
    else :
        raise EditCostError("This edit cost function doesn't exist, please see list_of_edit_cost_options for selecting a edit cost function")

def set_personal_edit_cost(edit_cost_constant = []) :
    """
        Sets an personal edit cost function to the environment.

        :param edit_cost_constant: The parameters you will add to the editCost, empty by default
        :type edit_cost_constant: list

        .. seealso:: list_of_edit_cost_options, set_edit_cost()
        .. note::You have to modify the C++ function to use it. Please see the documentation to add your Edit Cost function. 
    """
    setPersonalEditCost(edit_cost_constant)

def init(init_option = "EAGER_WITHOUT_SHUFFLED_COPIES") :
    """
        Initializes the environment with the chosen edit cost function and graphs.

        :param init_option: The name of the init option, "EAGER_WITHOUT_SHUFFLED_COPIES" by default
        :type init_option: string
        
        .. seealso:: list_of_init_options
        .. warning:: No modification were allowed after initialization. Try to make sure your choices is correct. You can though clear or add a graph, but recall init() after that. 
        .. note:: Try to make sure the option exists with list_of_init_options or choose no options, raise an error otherwise.
    """
    initB = init_option.encode('utf-8')
    if initB in list_of_init_options : 
        initEnv(initB)
    else :
        raise InitError("This init option doesn't exist, please see list_of_init_options for selecting an option. You can choose any options.")

def set_method(method, options="") :
    """
        Sets a computation method to the environment, if its exists. 

        :param method: The name of the computation method
        :param options: The options of the method (like bash options), an empty string by default
        :type method: string
        :type options: string
        
        .. seealso:: init_method(), list_of_method_options
        .. note:: Try to make sure the edit cost function exists with list_of_method_options, raise an error otherwise. Call init_method() after your set. 
    """
    methodB = method.encode('utf-8')
    if methodB in list_of_method_options :
        setMethod(methodB, options.encode('utf-8'))
    else :
        raise MethodError("This method doesn't exist, please see list_of_method_options for selecting a method")

def init_method() :
    """
        Inits the environment with the set method.

        .. seealso:: set_method(), list_of_method_options
        .. note:: Call this function after set the method. You can't launch computation or change the method after that. 
    """
    initMethod()

def get_init_time() :
    """
        Returns the initialization time.

        :return: The initialization time
        :rtype: double
    """
    return getInitime()

def run_method(g, h) :
    """
        Computes the edit distance between two graphs g and h, with the edit cost function and method computation selected.  

        :param g: The Id of the first graph to compare
        :param h: The Id of the second graph to compare
        :type g: size_t
        :type h: size_t
        
        .. seealso:: get_upper_bound(), get_lower_bound(),  get_forward_map(), get_backward_map(), get_runtime(), quasimetric_cost()
        .. note:: This function only compute the distance between two graphs, without returning a result. Use the differents function to see the result between the two graphs.  
    """
    runMethod(g,h)

def get_upper_bound(g,h) :
    """
        Returns the upper bound of the edit distance cost between two graphs g and h. 

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The upper bound of the edit distance cost
        :rtype: double
        
        .. seealso:: run_method(), get_lower_bound(),  get_forward_map(), get_backward_map(), get_runtime(), quasimetric_cost()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: The upper bound is equivalent to the result of the pessimist edit distance cost. Methods are heuristics so the library can't compute the real perfect result because it's NP-Hard problem.
    """
    return getUpperBound(g,h)

def get_lower_bound(g,h) :
    """
         Returns the lower bound of the edit distance cost between two graphs g and h. 

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The lower bound of the edit distance cost
        :rtype: double
        
        .. seealso:: run_method(), get_upper_bound(),  get_forward_map(), get_backward_map(), get_runtime(), quasimetric_cost()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: This function can be ignored, because lower bound doesn't have a crucial utility.    
    """
    return getLowerBound(g,h)

def get_forward_map(g,h) :
    """
        Returns the forward map (or the half of the adjacence matrix) between nodes of the two indicated graphs. 

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The forward map to the adjacence matrix between nodes of the two graphs
        :rtype: list[npy_uint32]
        
        .. seealso:: run_method(), get_upper_bound(), get_lower_bound(), get_backward_map(), get_runtime(), quasimetric_cost(), get_node_map(), get_assignment_matrix()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: I don't know how to connect the two map to reconstruct the adjacence matrix. Please come back when I know how it's work ! 
    """
    return getForwardMap(g,h)

def get_backward_map(g,h) :
    """
        Returns the backward map (or the half of the adjacence matrix) between nodes of the two indicated graphs. 

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The backward map to the adjacence matrix between nodes of the two graphs
        :rtype: list[npy_uint32]
        
        .. seealso:: run_method(), get_upper_bound(), get_lower_bound(),  get_forward_map(), get_runtime(), quasimetric_cost(), get_node_map(), get_assignment_matrix()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: I don't know how to connect the two map to reconstruct the adjacence matrix. Please come back when I know how it's work ! 
    """
    return getBackwardMap(g,h)

def get_node_image(g,h,node_id) :
    """
        Returns the node's image in the adjacence matrix, if it exists.   

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :param node_id: The ID of the node which you want to see the image
        :type g: size_t
        :type h: size_t
        :type node_id: size_t
        :return: The ID of the image node
        :rtype: size_t
        
        .. seealso:: run_method(), get_forward_map(), get_backward_map(), get_node_pre_image(), get_node_map(), get_assignment_matrix()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: Use BackwardMap's Node to find its images ! You can also use get_forward_map() and get_backward_map().     

    """
    return getNodeImage(g, h, node_id)

def get_node_pre_image(g,h,node_id) :
    """
        Returns the node's preimage in the adjacence matrix, if it exists.   

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :param node_id: The ID of the node which you want to see the preimage
        :type g: size_t
        :type h: size_t
        :type node_id: size_t
        :return: The ID of the preimage node
        :rtype: size_t
        
        .. seealso:: run_method(), get_forward_map(), get_backward_map(), get_node_image(), get_node_map(), get_assignment_matrix()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: Use ForwardMap's Node to find its images ! You can also use get_forward_map() and get_backward_map().     

    """
    return getNodePreImage(g, h, node_id)

def get_dummy_node() :
    """
        Returns the ID of a dummy node.

        :return: The ID of the dummy node (18446744073709551614 for my computer, the hugest number possible)
        :rtype: size_t
        
        .. note:: A dummy node is used when a node isn't associated to an other node.      
    """
    return getDummyNode()

def get_node_map(g,h) : 
    """
        Returns the Node Map, like C++ NodeMap.   

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The Node Map between the two selected graph. 
        :rtype: list[tuple(size_t, size_t)]
        
        .. seealso:: run_method(), get_forward_map(), get_backward_map(), get_node_image(), get_node_pre_image(), get_assignment_matrix()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: This function creates datas so use it if necessary, however you can understand how assignement works with this example.     
    """
    return getNodeMap(g, h)

def get_assignment_matrix(g,h) :
    """
        Returns the Assignment Matrix between two selected graphs g and h.   

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The Assignment Matrix between the two selected graph. 
        :rtype: list[list[int]]
        
        .. seealso:: run_method(), get_forward_map(), get_backward_map(), get_node_image(), get_node_pre_image(), get_node_map()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: This function creates datas so use it if necessary.     
    """
    return getAssignmentMatrix(g, h)
        

def get_all_map(g,h) :
    """
         Returns a vector which contains the forward and the backward maps between nodes of the two indicated graphs. 

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The forward and backward maps to the adjacence matrix between nodes of the two graphs
        :rtype: list[list[npy_uint32]]
        
        .. seealso:: run_method(), get_upper_bound(), get_lower_bound(),  get_forward_map(), get_backward_map(), get_runtime(), quasimetric_cost()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: This function duplicates data so please don't use it. I also don't know how to connect the two map to reconstruct the adjacence matrix. Please come back when I know how it's work !  
    """
    return getAllMap(g,h)

def get_runtime(g,h) :
    """
        Returns the runtime to compute the edit distance cost between two graphs g and h  

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: The runtime of the computation of edit distance cost between the two selected graphs
        :rtype: double
        
        .. seealso:: run_method(), get_upper_bound(), get_lower_bound(),  get_forward_map(), get_backward_map(), quasimetric_cost()
        .. warning:: run_method() between the same two graph must be called before this function. 
        .. note:: Python is a bit longer than C++ due to the functions's encapsulate.    
    """
    return getRuntime(g,h)

def quasimetric_cost() :
    """
        Checks and returns if the edit costs are quasimetric. 

        :param g: The Id of the first compared graph 
        :param h: The Id of the second compared graph
        :type g: size_t
        :type h: size_t
        :return: True if it's verified, False otherwise
        :rtype: bool
        
        .. seealso:: run_method(), get_upper_bound(), get_lower_bound(),  get_forward_map(), get_backward_map(), get_runtime()
        .. warning:: run_method() between the same two graph must be called before this function. 
    """
    return quasimetricCosts()

def hungarian_LSAP(matrix_cost) :
    """
        Applies the hungarian algorithm (LSAP) on a matrix Cost. 

        :param matrix_cost: The matrix Cost  
        :type matrix_cost: vector[vector[size_t]]
        :return: The values of rho, varrho, u and v, in this order
        :rtype: vector[vector[size_t]]
        
        .. seealso:: hungarian_LSAPE() 
    """
    return hungarianLSAP(matrix_cost)

def hungarian_LSAPE(matrix_cost) :
    """
        Applies the hungarian algorithm (LSAPE) on a matrix Cost. 

        :param matrix_cost: The matrix Cost 
        :type matrix_cost: vector[vector[double]]
        :return: The values of rho, varrho, u and v, in this order
        :rtype: vector[vector[double]]
        
        .. seealso:: hungarian_LSAP() 
    """
    return hungarianLSAPE(matrix_cost)

#####################################################################
##LISTS OF EDIT COST FUNCTIONS, METHOD COMPUTATION AND INIT OPTIONS##
#####################################################################

list_of_edit_cost_options = get_edit_cost_options()
list_of_method_options = get_method_options()
list_of_init_options = get_init_options()


#####################
##ERRORS MANAGEMENT##
#####################

class Error(Exception):
    """
        Class for error's management. This one is general. 
    """
    pass

class EditCostError(Error) :
    """
        Class for Edit Cost Error. Raise an error if an edit cost function doesn't exist in the library (not in list_of_edit_cost_options).

        :attribute message: The message to print when an error is detected.
        :type message: string
    """
    def __init__(self, message):
        """
            Inits the error with its message. 

            :param message: The message to print when the error is detected
            :type message: string
        """
        self.message = message
    
class MethodError(Error) :
    """
        Class for Method Error. Raise an error if a computation method doesn't exist in the library (not in list_of_method_options).

        :attribute message: The message to print when an error is detected.
        :type message: string
    """
    def __init__(self, message):
        """
            Inits the error with its message. 

            :param message: The message to print when the error is detected
            :type message: string
        """
        self.message = message

class InitError(Error) :
    """
        Class for Init Error. Raise an error if an init option doesn't exist in the library (not in list_of_init_options).

        :attribute message: The message to print when an error is detected.
        :type message: string
    """
    def __init__(self, message):
        """
            Inits the error with its message. 

            :param message: The message to print when the error is detected
            :type message: string
        """
        self.message = message


#########################################
##PYTHON FUNCTIONS FOR SOME COMPUTATION##
#########################################

def encode_your_map(map) :
    """
        Encodes a string dictionnary to utf-8 for C++ functions

        :param map: The map to encode
        :type map: dict{string : string}
        :return: The encoded map
        :rtype: dict{'b'string : 'b'string}

        .. note:: This function is used for type connection.  
        
    """
    res = {}
    for key, value in map.items():
        res[key.encode('utf-8')] = value.encode('utf-8')
    return res

def add_random_graph(name, classe, list_of_nodes, list_of_edges, ignore_duplicates=True) :
    """
        Add a Graph (not GXL) on the environment. Be careful to respect the same format as GXL graphs for labelling nodes and edges. 

        :param name: The name of the graph to add, can be an empty string
        :param classe: The classe of the graph to add, can be an empty string
        :param list_of_nodes: The list of nodes to add
        :param list_of_edges: The list of edges to add
        :param ignore_duplicates: If True, duplicate edges are ignored, otherwise it's raise an error if an existing edge is added. True by default
        :type name: string
        :type classe: string
        :type list_of_nodes: list[tuple(size_t, dict{string : string})]
        :type list_of_edges: list[tuple(tuple(size_t,size_t), dict{string : string})]
        :type ignore_duplicates: bool
        :return: The ID of the newly added graphe
        :rtype: size_t

        .. note:: The graph must respect the GXL structure. Please see how a GXL graph is construct.  
        
    """
    id = add_graph(name, classe)
    for node in list_of_nodes :
        add_node(id, node[0], node[1])
    for edge in list_of_edges :
        add_edge(id, edge[0], edge[1], edge[2], ignore_duplicates)
    return id

def add_nx_graph(g, classe, ignore_duplicates=True) :
    """
        Add a Graph (made by networkx) on the environment. Be careful to respect the same format as GXL graphs for labelling nodes and edges. 

        :param g: The graph to add (networkx graph)
        :param ignore_duplicates: If True, duplicate edges are ignored, otherwise it's raise an error if an existing edge is added. True by default
        :type g: networkx.graph
        :type ignore_duplicates: bool
        :return: The ID of the newly added graphe
        :rtype: size_t

        .. note:: The NX graph must respect the GXL structure. Please see how a GXL graph is construct.  
        
    """
    id = add_graph(g.name, classe)
    for node in g.nodes :
        add_node(id, str(node), g.node[node])
    for edge in g.edges :
        add_edge(id, str(edge[0]), str(edge[1]), g.get_edge_data(edge[0],edge[1]), ignore_duplicates)
    return id


def compute_ged_on_two_graphs(g1,g2, edit_cost, method, options, init_option = "EAGER_WITHOUT_SHUFFLED_COPIES") :
    """
        Computes the edit distance between two NX graphs. 
        
        :param g1: The first graph to add and compute
        :param g2: The second graph to add and compute
        :param edit_cost: The name of the edit cost function
        :param method: The name of the computation method
        :param options: The options of the method (like bash options), an empty string by default
        :param init_option:  The name of the init option, "EAGER_WITHOUT_SHUFFLED_COPIES" by default
        :type g1: networksx.graph
        :type g2: networksx.graph
        :type edit_cost: string
        :type method: string
        :type options: string
        :type init_option: string
        :return: The edit distance between the two graphs and the nodeMap between them. 
        :rtype: double, list[tuple(size_t, size_t)]

        .. seealso:: list_of_edit_cost_options, list_of_method_options, list_of_init_options 
        .. note:: Make sure each parameter exists with your architecture and these lists :  list_of_edit_cost_options, list_of_method_options, list_of_init_options. The structure of graphs must be similar as GXL. 
        
    """
    if is_initialized() :
        restart_env()

    g = add_nx_graph(g1, "")
    h = add_nx_graph(g2, "")

    set_edit_cost(edit_cost)
    init(init_option)
    
    set_method(method, options)
    init_method()

    resDistance = 0
    resMapping = []
    run_method(g,h)
    resDistance = get_upper_bound(g,h)
    resMapping = get_node_map(g,h)

    return resDistance, resMapping

def compute_edit_distance_on_nx_graphs(dataset, classes, edit_cost, method, options, init_option = "EAGER_WITHOUT_SHUFFLED_COPIES") :
    """

        Computes all the edit distance between each NX graphs on the dataset. 
        
        :param dataset: The list of graphs to add and compute
        :param classes: The classe of all the graph, can be an empty string
        :param edit_cost: The name of the edit cost function
        :param method: The name of the computation method
        :param options: The options of the method (like bash options), an empty string by default
        :param init_option:  The name of the init option, "EAGER_WITHOUT_SHUFFLED_COPIES" by default
        :type dataset: list[networksx.graph]
        :type classes: string
        :type edit_cost: string
        :type method: string
        :type options: string
        :type init_option: string
        :return: Two matrix, the first with edit distances between graphs and the second the nodeMap between graphs. The result between g and h is one the [g][h] coordinates.
        :rtype: list[list[double]], list[list[list[tuple(size_t, size_t)]]]

        .. seealso:: list_of_edit_cost_options, list_of_method_options, list_of_init_options
        .. note:: Make sure each parameter exists with your architecture and these lists :  list_of_edit_cost_options, list_of_method_options, list_of_init_options. The structure of graphs must be similar as GXL. 
        
    """
    if is_initialized() :
        restart_env()

    print("Loading graphs in progress...")
    for graph in dataset :
        add_nx_graph(graph, classes)
    listID = graph_ids()
    print("Graphs loaded ! ")
    print("Number of graphs = " + str(listID[1]))

    set_edit_cost(edit_cost)
    print("Initialization in progress...")
    init(init_option)
    print("Initialization terminated !")
    
    set_method(method, options)
    init_method()

    resDistance = [[]]
    resMapping = [[]]
    for g in range(listID[0], listID[1]) :
        print("Computation between graph " + str(g) + " with all the others including himself.")
        for h in range(listID[0], listID[1]) :
            #print("Computation between graph " + str(g) + " and graph " + str(h))
            run_method(g,h)
            resDistance[g][h] = get_upper_bound(g,h)
            resMapping[g][h] = get_node_map(g,h)

    print("Finish ! The return contains edit distances and NodeMap but you can check the result with graphs'ID until you restart the environment")
    return resDistance, resMapping
    
    
    
def compute_edit_distance_on_GXl_graphs(path_folder, path_XML, edit_cost, method, options="", init_option = "EAGER_WITHOUT_SHUFFLED_COPIES") :
    """
        Computes all the edit distance between each GXL graphs on the folder and the XMl file. 
        
        :param path_folder: The folder's path which contains GXL graphs
        :param path_XML: The XML's path which indicates which graphes you want to load
        :param edit_cost: The name of the edit cost function
        :param method: The name of the computation method
        :param options: The options of the method (like bash options), an empty string by default
        :param init_option:  The name of the init option, "EAGER_WITHOUT_SHUFFLED_COPIES" by default
        :type path_folder: string
        :type path_XML: string
        :type edit_cost: string
        :type method: string
        :type options: string
        :type init_option: string
        :return: The list of the first and last-1 ID of graphs
        :rtype: tuple(size_t, size_t)

        .. seealso:: list_of_edit_cost_options, list_of_method_options, list_of_init_options
        .. note:: Make sure each parameter exists with your architecture and these lists : list_of_edit_cost_options, list_of_method_options, list_of_init_options. 
        
    """

    if is_initialized() :
        restart_env()

    print("Loading graphs in progress...")
    load_GXL_graphs(path_folder, path_XML)
    listID = graph_ids()
    print("Graphs loaded ! ")
    print("Number of graphs = " + str(listID[1]))

    set_edit_cost(edit_cost)
    print("Initialization in progress...")
    init(init_option)
    print("Initialization terminated !")
    
    set_method(method, options)
    init_method()

    #res = []
    for g in range(listID[0], listID[1]) :
        print("Computation between graph " + str(g) + " with all the others including himself.")
        for h in range(listID[0], listID[1]) :
            #print("Computation between graph " + str(g) + " and graph " + str(h))
            run_method(g,h)
            #res.append((get_upper_bound(g,h), get_node_map(g,h), get_runtime(g,h)))
            
    #return res

    print ("Finish ! You can check the result with each ID of graphs ! There are in the return")
    print ("Please don't restart the environment or recall this function, you will lose your results !")
    return listID


    
