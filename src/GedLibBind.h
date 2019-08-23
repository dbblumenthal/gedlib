/****************************************************************************
 *                                                                          *
 *   Copyright (C) 2019 by Natacha Lambert and David B. Blumenthal          *
 *                                                                          *
 *   This file should be used by Python.                                    *
 * 	 Please call the Python module if you want to use GedLib with this code.* 
 *                                                                          *
 * 	 Otherwise, you can directly use GedLib for C++.                        *
 *                                                                          *
 ***************************************************************************/
 
/*!
 * @file GedLibBind.h
 * @brief Functions declaration to call easly GebLib in Python without Gedlib's types
 */
 
//Include standard libraries.
#include <string>
#include <vector>
#include <map>
#include <list>

std::vector<std::string> getEditCostStringOptions();  //!< List of available edit cost functions readable by Python.  
std::vector<std::string> getMethodStringOptions();  //!< List of available computation methods readable by Python.  
std::vector<std::string> getInitStringOptions();  //!<List of available initilaization options readable by Python.

/*!
 * @brief Tests if the environment is initialized or not. 
 * @return Boolean @p true if the environment is initialized and @p false otherwise.
 */
bool isInitialized();

/*!
 * @brief Restart the environment (recall a new empty environment).
 */
void restartEnv();

/*!
 * @brief Loads graphs given in the [GXL file format](http://www.gupro.de/GXL/).
 * @param[in] pathFolder The path to the directory containing the graphs.
 * @param[in] pathXML The path to a XML file thats lists the graphs contained in @p pathFolder that should be loaded.
 */
void loadGXLGraph(std::string pathFolder, std::string pathXML);

/*!
 * @brief Provides access to the IDs of the graphs contained in the environment.
 * @return Pair <tt>(ID of first graphs, ID of last graph + 1)</tt> of graph IDs.
 * If both entries equal 0, the environment does not contain any graphs.
 */
std::pair<std::size_t,std::size_t> getGraphIds();

/*!
 * @brief Returns the list of graphs IDs which are loaded in the environment. 
 * @return A vector which contains all the graphs Ids. 
 */
std::vector<std::size_t> getAllGraphIds();

/*!
 * @brief Returns the graph class.
 * @param[in] id ID of an input graph that has been added to the environment.
 * @return Class of the input graph.
 */
std::string getGraphClass(std::size_t id);

/*!
 * @brief Returns the graph name.
 * @param[in] id ID of an input graph that has been added to the environment.
 * @return Name of the input graph.
 */
std::string getGraphName(std::size_t id);

/*!
 * @brief Adds a new uninitialized graph to the environment. Call initEnv() after calling this method.
 * @param[in] name The name of the added graph. Empty if not specified.
 * @param[in] class The class of the added graph. Empty if not specified.
 * @return The ID of the newly added graph.
 */
std::size_t addGraph(std::string name, std::string classe);

/*!
 * @brief Adds a labeled node.
 * @param[in] graphId ID of graph that has been added to the environment.
 * @param[in] nodeId The user-specific ID of the vertex that has to be added.
 * @param[in] nodeLabel The label of the vertex that has to be added.
 */
void addNode(std::size_t graphId, std::string nodeId, std::map<std::string, std::string> nodeLabel);

/*!
 * @brief Adds a labeled edge.
 * @param[in] graphId ID of graph that has been added to the environment.
 * @param[in] tail The user-specific ID of the tail of the edge that has to be added.
 * @param[in] head The user-specific ID of the head of the edge that has to be added.
 * @param[in] edgeLabel The label of the vertex that has to be added. 
 * @param[in] ignoreDuplicates If @p true, duplicate edges are ignores. Otherwise, an exception is thrown if an existing edge is added to the graph.
 */
void addEdge(std::size_t graphId, std::string tail, std::string head, std::map<std::string, std::string> edgeLabel, bool ignoreDuplicates = true);

/*!
 * @brief Clears and de-initializes a graph that has previously been added to the environment. Call initEnv() after calling this method.
 * @param[in] graphId ID of graph that has to be cleared.
 */
void clearGraph(std::size_t graphId);

/*!
 * @brief Returns the internal Id of a graph, selected by its ID.
 * @param[in] graphId ID of an input graph that has been added to the environment.
 * @return The internal ID of the selected graph
 */
std::size_t getGraphInternalId(std::size_t graphId);

/*!
 * @brief Returns all the number of nodes on a graph, selected by its ID.
 * @param[in] graphId ID of an input graph that has been added to the environment.
 * @return The number of nodes on the selected graph
 */
std::size_t getGraphNumNodes(std::size_t graphId);

/*!
 * @brief Returns all the number of edges on a graph, selected by its ID.
 * @param[in] graphId ID of an input graph that has been added to the environment.
 * @return The number of edges on the selected graph
 */
std::size_t getGraphNumEdges(std::size_t graphId);

/*!
 * @brief Returns all th Ids of nodes on a graph, selected by its ID.
 * @param[in] graphId ID of an input graph that has been added to the environment.
 * @return The list of IDs's nodes on the selected graph
 */
std::vector<std::string> getGraphOriginalNodeIds(std::size_t graphId);

/*!
 * @brief Returns all the labels of nodes on a graph, selected by its ID.
 * @param[in] graphId ID of an input graph that has been added to the environment.
 * @return The list of labels's nodes on the selected graph
 */
std::vector<std::map<std::string, std::string>> getGraphNodeLabels(std::size_t graphId);

/*!
 * @brief Returns all the edges on a graph, selected by its ID.
 * @param[in] graphId ID of an input graph that has been added to the environment.
 * @return The list of edges on the selected graph
 */
std::map<std::pair<std::size_t, std::size_t>, std::map<std::string, std::string>> getGraphEdges(std::size_t graphId);

/*!
 * @brief Returns the adjacence list of a graph, selected by its ID.
 * @param[in] graphId ID of an input graph that has been added to the environment.
 * @return The adjacence list of the selected graph
 */
std::vector<std::vector<std::size_t>> getGraphAdjacenceMatrix(std::size_t graphId);

/*!
 * @brief Sets the edit costs to one of the predefined edit costs.
 * @param[in] editCost Select one of the predefined edit costs.
 * @param[in] editCostConstants Parameters for the edit cost, empty by default.
 */
void setEditCost(std::string editCost, std::vector<double> editCostConstants = {});

/*!
 * @brief Sets the edit costs to a personal Edit Cost Class.
 * @param[in] editCostConstants Parameters for the edit cost, empty by default.
 * @note You have to add your class, which should inherit from EditCost class, in the function. After that, you can compile and use it in Python
 */
void setPersonalEditCost(std::vector<double> editCostConstants = {});

/*!
 * @brief Initializes the environment.
 * @param[in] initOption Select initialization options.
 */
void initEnv(std::string initOption);

/*!
 * @brief Sets the GEDMethod to be used by run_method().
 * @param[in] method Select the method that is to be used.
 * @param[in] options An options string of the form @"[--@<option@> @<arg@>] [...]@" passed to the selected method.
 */
void setMethod(std::string method, std::string options);

/*!
 * @brief Initializes the method specified by call to set_method().
 */
void initMethod();

/*!
 * @brief Returns initialization time.
 * @return Runtime of the last call to init_method().
 */
double getInitime();

/*!
 * @brief Runs the GED method specified by call to set_method() between the graphs with IDs @p g and @p h.
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 */
void runMethod(std::size_t g, std::size_t h );

/*!
 * @brief Returns upper bound for edit distance between the input graphs.
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return Upper bound computed by the last call to run_method() with arguments @p g and @p h.
 */
double getUpperBound(std::size_t g, std::size_t h);

/*!
 * @brief Returns lower bound for edit distance between the input graphs.
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return Lower bound computed by the last call to run_method() with arguments @p g and @p h.
 */
double getLowerBound(std::size_t g,std::size_t h);

/*!
 * @brief  Returns the forward map between nodes of the two indicated graphs. 
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return The forward map to the adjacence matrix computed by the last call to run_method() with arguments @p g and @p h.
 */
std::vector<long unsigned int> getForwardMap(std::size_t g, std::size_t h);

/*!
 * @brief  Returns the backward map between nodes of the two indicated graphs. 
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return The backward map to the adjacence matrix computed by the last call to run_method() with arguments @p g and @p h.
 */
std::vector<long unsigned int> getBackwardMap(std::size_t g, std::size_t h);

/*!
 * @brief Returns image of a node.
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @param[in] nodeId Node whose image is to be returned.
 * @return Node to which node @p node is assigned.
 */
std::size_t getNodeImage(std::size_t g, std::size_t h, std::size_t nodeId);

/*!
 * @brief Returns pre-image of a node.
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @param[in] nodeId Node whose pre-image is to be returned.
 * @return Node to which node @p node is assigned.
 */
std::size_t getNodePreImage(std::size_t g, std::size_t h, std::size_t nodeId);

/*!
 * @brief Returns a dummy node.
 * @return ID of dummy node.
 */
std::size_t getDummyNode();

/*!
 * @brief Returns node map between the input graphs. This function duplicates datas. 
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return Node map computed by the last call to run_method() with arguments @p g and @p h.
 */
std::vector<std::pair<std::size_t, std::size_t>> getNodeMap(std::size_t g, std::size_t h);

/*!
 * @brief Returns assignment matrix between the input graphs. This function duplicates datas. 
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return Assignment matrix computed by the last call to run_method() with arguments @p g and @p h.
 */
std::vector<std::vector<int>> getAssignmentMatrix(std::size_t g, std::size_t h);

/*!
 * @brief  Returns a vector which contains the forward and the backward maps between nodes of the two indicated graphs. 
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return The forward and backward maps to the adjacence matrix computed by the last call to run_method() with arguments @p g and @p h.
 */
std::vector<std::vector<unsigned long int>> getAllMap(std::size_t g, std::size_t h);

/*!
 * @brief Returns runtime.
 * @param[in] g ID of an input graph that has been added to the environment.
 * @param[in] h ID of an input graph that has been added to the environment.
 * @return Runtime of last call to run_method() with arguments @p g and @p h.
 */
double getRuntime(std::size_t g, std::size_t h);

/*!
 * @brief Checks if the edit costs are quasimetric.
 * @return Boolean @p true if the edit costs are quasimetric and @p false, otherwise.
 */
bool quasimetricCosts();

/*!
 * @brief Applies the hungarian algorithm (LSAP) to a matrix cost.
 * @param[in] matrixCost The matrix cost.
 * @return the values of rho, varrho, u and v, in this order.
 */
std::vector<std::vector<size_t>> hungarianLSAP(std::vector<std::vector<std::size_t>> matrixCost);

/*!
 * @brief Applies the hungarian algorithm (LSAPE) to a matrix cost.
 * @param[in] matrixCost The matrix cost.
 * @return the values of rho, varrho, u and v, in this order.
 */
std::vector<std::vector<double>>  hungarianLSAPE(std::vector<std::vector<double>> matrixCost);

