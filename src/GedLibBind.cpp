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
 * @file GedLibBind.cpp
 * @brief Functions definition to call easly GebLib in Python without Gedlib's types
 */

//Include standard libraries + GedLib library
#include <iostream>
#include "GedLibBind.h"
#include "../include/gedlib-master/src/env/ged_env.hpp"
//#include "../include/gedlib-master/median/src/median_graph_estimator.hpp"

using namespace std;

//Definition of types and templates used in this code for my human's memory :). 
//ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> env;
//template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel> struct ExchangeGraph

//typedef std::map<std::string, std::string> GXLLabel;
//typedef std::string GXLNodeID;

ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env; //Environment variable


bool initialized = false; //Initialization boolean (because Env has one but not accessible). 

bool isInitialized(){
	return initialized;
}

//!< List of available edit cost functions readable by Python.  
std::vector<std::string> editCostStringOptions = { 
	"CHEM_1", 
	"CHEM_2", 
	"CMU", 
	"GREC_1", 
	"GREC_2", 
	"LETTER", 
	"FINGERPRINT", 
	"PROTEIN", 
	"CONSTANT" 
};

std::vector<std::string> getEditCostStringOptions(){
	return editCostStringOptions;
}

//!< Map of available edit cost functions between enum type in C++ and string in Python  
std::map<std::string, ged::Options::EditCosts> editCostOptions = {
	{"CHEM_1", ged::Options::EditCosts::CHEM_1},
	{"CHEM_2", ged::Options::EditCosts::CHEM_2},
	{"CMU", ged::Options::EditCosts::CMU},
	{"GREC_1", ged::Options::EditCosts::GREC_1},
	{"GREC_2", ged::Options::EditCosts::GREC_2},
	{"LETTER", ged::Options::EditCosts::LETTER},
	{"FINGERPRINT", ged::Options::EditCosts::FINGERPRINT},
	{"PROTEIN", ged::Options::EditCosts::PROTEIN},
	{"CONSTANT", ged::Options::EditCosts::CONSTANT}	
};

 //!< List of available computation methods readable by Python.  
std::vector<std::string> methodStringOptions = { 
	"BRANCH",
	"BRANCH_FAST",
	"BRANCH_TIGHT",
	"BRANCH_UNIFORM",
	"BRANCH_COMPACT",
	"PARTITION",
	"HYBRID",
	"RING",
	"ANCHOR_AWARE_GED",
	"WALKS",
	"IPFP",
	"BIPARTITE",
	"SUBGRAPH",
	"NODE",
	"RING_ML",
	"BIPARTITE_ML",
	"REFINE",
	"BP_BEAM",
	"SIMULATED_ANNEALING",
	"HED",
	"STAR"				 
}; 

std::vector<std::string> getMethodStringOptions(){
	return methodStringOptions;
}

//!< Map of available computation methods readables between enum type in C++ and string in Python  
std::map<std::string, ged::Options::GEDMethod> methodOptions = {
	{"BRANCH", ged::Options::GEDMethod::BRANCH},
	{"BRANCH_FAST", ged::Options::GEDMethod::BRANCH_FAST},
	{"BRANCH_TIGHT", ged::Options::GEDMethod::BRANCH_TIGHT},
	{"BRANCH_UNIFORM", ged::Options::GEDMethod::BRANCH_UNIFORM},
	{"BRANCH_COMPACT", ged::Options::GEDMethod::BRANCH_COMPACT},
	{"PARTITION", ged::Options::GEDMethod::PARTITION},
	{"HYBRID", ged::Options::GEDMethod::HYBRID},
	{"RING", ged::Options::GEDMethod::RING},
	{"ANCHOR_AWARE_GED", ged::Options::GEDMethod::ANCHOR_AWARE_GED},
	{"WALKS", ged::Options::GEDMethod::WALKS},
	{"IPFP", ged::Options::GEDMethod::IPFP},
	{"BIPARTITE", ged::Options::GEDMethod::BIPARTITE},
	{"SUBGRAPH", ged::Options::GEDMethod::SUBGRAPH},
	{"NODE", ged::Options::GEDMethod::NODE},
	{"RING_ML", ged::Options::GEDMethod::RING_ML},
	{"BIPARTITE_ML",ged::Options::GEDMethod::BIPARTITE_ML},
	{"REFINE",ged::Options::GEDMethod::REFINE},
	{"BP_BEAM", ged::Options::GEDMethod::BP_BEAM},
	{"SIMULATED_ANNEALING", ged::Options::GEDMethod::SIMULATED_ANNEALING},
	{"HED", ged::Options::GEDMethod::HED},
	{"STAR"	, ged::Options::GEDMethod::STAR},	
};

//!<List of available initilaization options readable by Python.
std::vector<std::string> initStringOptions = { 
	"LAZY_WITHOUT_SHUFFLED_COPIES", 
	"EAGER_WITHOUT_SHUFFLED_COPIES", 
	"LAZY_WITH_SHUFFLED_COPIES", 
	"EAGER_WITH_SHUFFLED_COPIES"
};

std::vector<std::string> getInitStringOptions(){
	return initStringOptions;
}

//!< Map of available initilaization options readables between enum type in C++ and string in Python 
std::map<std::string, ged::Options::InitType> initOptions = {
	{"LAZY_WITHOUT_SHUFFLED_COPIES", ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES},
	{"EAGER_WITHOUT_SHUFFLED_COPIES", ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES},
	{"LAZY_WITH_SHUFFLED_COPIES", ged::Options::InitType::LAZY_WITH_SHUFFLED_COPIES},
	{"EAGER_WITH_SHUFFLED_COPIES", ged::Options::InitType::EAGER_WITH_SHUFFLED_COPIES}
};

void restartEnv(){
	env = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>();
	initialized = false;
}

void loadGXLGraph(std::string pathFolder, std::string pathXML){
	 std::vector<ged::GEDGraph::GraphID> tmp_graph_ids(env.load_gxl_graphs(pathFolder, pathXML));
}

std::pair<std::size_t,std::size_t> getGraphIds(){
	return env.graph_ids();
}

std::vector<std::size_t> getAllGraphIds(){
	std::vector<std::size_t> listID;
	for (std::size_t i = env.graph_ids().first; i != env.graph_ids().second; i++){
		listID.push_back(i);
    }
	return listID;
}

std::string getGraphClass(std::size_t id){
	return env.get_graph_class(id);
}

std::string getGraphName(std::size_t id){
	return env.get_graph_name(id);
}

std::size_t addGraph(std::string name, std::string classe){
	ged::GEDGraph::GraphID newId = env.add_graph(name, classe); 
	initialized = false;
	return std::stoi(std::to_string(newId));
}

void addNode(std::size_t graphId, std::string nodeId, std::map<std::string, std::string> nodeLabel){
	env.add_node(graphId, nodeId, nodeLabel);
	initialized = false;
}

/*void addEdge(std::size_t graphId, ged::GXLNodeID tail, ged::GXLNodeID head, ged::GXLLabel edgeLabel){
	env.add_edge(graphId, tail, head, edgeLabel);
}*/

void addEdge(std::size_t graphId, std::string tail, std::string head, std::map<std::string, std::string> edgeLabel, bool ignoreDuplicates){
	env.add_edge(graphId, tail, head, edgeLabel, ignoreDuplicates);
	initialized = false;
}

void clearGraph(std::size_t graphId){
	env.clear_graph(graphId);
	initialized = false;
}

/*!
 * @brief Returns ged::ExchangeGraph representation.
 * @param graphId ID of the selected graph.
 * @return ged::ExchangeGraph representation of the selected graph.
 */
ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> getGraph(std::size_t graphId){
	return env.get_graph(graphId);
}

std::size_t getGraphInternalId(std::size_t graphId){
	return getGraph(graphId).id;
}

std::size_t getGraphNumNodes(std::size_t graphId){
	return getGraph(graphId).num_nodes;
}

std::size_t getGraphNumEdges(std::size_t graphId){
	return getGraph(graphId).num_edges;
}

std::vector<std::string> getGraphOriginalNodeIds(std::size_t graphId){
	return getGraph(graphId).original_node_ids;
}

std::vector<std::map<std::string, std::string>> getGraphNodeLabels(std::size_t graphId){
	return getGraph(graphId).node_labels;
}

std::map<std::pair<std::size_t, std::size_t>, std::map<std::string, std::string>> getGraphEdges(std::size_t graphId){
	return getGraph(graphId).edge_labels;
}

std::vector<std::vector<std::size_t>> getGraphAdjacenceMatrix(std::size_t graphId){
	return getGraph(graphId).adj_matrix;
}

/*!
 * @brief Returns the enum EditCost which correspond to the string parameter
 * @param editCost Select one of the predefined edit costs in the list.
 * @return The edit cost function which correspond in the edit cost functions map. 
 */
ged::Options::EditCosts translateEditCost(std::string editCost){
	 for (int i = 0; i != editCostStringOptions.size(); i++){
		 if (editCostStringOptions[i] == editCost){
			 return editCostOptions[editCostStringOptions[i]];
		 } 
	 }
	 return ged::Options::EditCosts::CONSTANT;
}

void setEditCost(std::string editCost, std::vector<double> editCostConstants){
	env.set_edit_costs(translateEditCost(editCost), editCostConstants);
}

void setPersonalEditCost(std::vector<double> editCostConstants){
	//env.set_edit_costs(Your EditCost Class(editCostConstants));
}

void initEnv(){
	env.init();
	initialized = true;
}

/*!
 * @brief Returns the enum IniType which correspond to the string parameter
 * @param initOption Select initialization options.
 * @return The init Type which correspond in the init options map. 
 */
ged::Options::InitType translateInitOptions(std::string initOption){
	 for (int i = 0; i != initStringOptions.size(); i++){
		 if (initStringOptions[i] == initOption){
			 return initOptions[initStringOptions[i]];
		 } 
	 }
	 return ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES;
}

void initEnv(std::string initOption){
	env.init(translateInitOptions(initOption));
	initialized = true;
}

/*!
 * @brief Returns the enum Method which correspond to the string parameter
 * @param method Select the method that is to be used.
 * @return The computation method which correspond in the edit cost functions map. 
 */
ged::Options::GEDMethod translateMethod(std::string method){
	 for (int i = 0; i != methodStringOptions.size(); i++){
		 if (methodStringOptions[i] == method){
			 return methodOptions[methodStringOptions[i]];
		 } 
	 }
	 return ged::Options::GEDMethod::STAR;
}

void setMethod(std::string method, std::string options){
	env.set_method(translateMethod(method),options);
}

void initMethod(){
	env.init_method();
}

double getInitime(){
	return env.get_init_time();
}

void runMethod(std::size_t g, std::size_t h){
	env.run_method(g, h);
}

double getUpperBound(std::size_t g, std::size_t h){
	return env.get_upper_bound(g, h);
}

double getLowerBound(std::size_t g, std::size_t h){
	return env.get_lower_bound(g, h);
}

std::vector<long unsigned int> getForwardMap(std::size_t g, std::size_t h){
	return env.get_node_map(g,h).get_forward_map(); 
}

std::vector<long unsigned int> getBackwardMap(std::size_t g, std::size_t h){
	return env.get_node_map(g,h).get_backward_map(); 
}

std::size_t getNodeImage(std::size_t g, std::size_t h, std::size_t nodeId){
	return env.get_node_map(g,h).image(nodeId);
}

std::size_t getNodePreImage(std::size_t g, std::size_t h, std::size_t nodeId){
	return env.get_node_map(g,h).pre_image(nodeId);
}

std::size_t getDummyNode(){
	return ged::GEDGraph::dummy_node();
}

std::vector<pair<std::size_t, std::size_t>> getNodeMap(std::size_t g, std::size_t h){
	std::vector<pair<std::size_t, std::size_t>> res; 
	std::vector<ged::NodeMap::Assignment> relation;
	env.get_node_map(g,h).as_relation(relation);
	for (const auto & assignment : relation) {
		res.push_back(std::make_pair(assignment.first, assignment.second)); 
	}
	return res;
}

std::vector<std::vector<int>> getAssignmentMatrix(std::size_t g, std::size_t h){
	std::vector<std::vector<int>> res;
	for(std::size_t i = 0; i != getForwardMap(g, h).size(); i++){
		std::vector<int> newLine;
		bool have1 = false;
		for(std::size_t j = 0; j != getBackwardMap(g, h).size(); j++){
			if (getNodeImage(g, h, i) == j){
				newLine.push_back(1);
				have1 = true;
			}
			else{
				newLine.push_back(0);
			}
		}
		if(have1){
			newLine.push_back(0);
		}
		else{
			newLine.push_back(1);
		}
		res.push_back(newLine);
	}
	std::vector<int> lastLine;
	for (size_t k = 0; k != getBackwardMap(g,h).size(); k++){
		if (getBackwardMap(g,h)[k] ==  ged::GEDGraph::dummy_node()){
			lastLine.push_back(1);
		}
		else{
			lastLine.push_back(0);
		}
	}
	res.push_back(lastLine);
	return res;
}

std::vector<std::vector<unsigned long int>> getAllMap(std::size_t g, std::size_t h){
	std::vector<std::vector<unsigned long int>> res; 
	res.push_back(getForwardMap(g, h));
	res.push_back(getBackwardMap(g,h)); 
	return res;
}

double getRuntime(std::size_t g, std::size_t h){
	return env.get_runtime(g, h);
}

bool quasimetricCosts(){
	return env.quasimetric_costs();
}

/*!
 * @brief Returns the vector of values which correspond to the pointer parameter.
 * @param pointer The size_t pointer to convert. 
 * @return The vector which contains the pointer's values. 
 */
std::vector<size_t> translatePointer(std::size_t* pointer, std::size_t dataSize ){
	std::vector<size_t> res;
	for(int i = 0; i < dataSize; i++){
		res.push_back(pointer[i]);
	}
	return res;
}

/*!
 * @brief Returns the vector of values which correspond to the pointer parameter.
 * @param pointer The double pointer to convert. 
 * @return The vector which contains the pointer's values. 
 */
std::vector<double> translatePointer(double* pointer, std::size_t dataSize ){
	std::vector<double> res;
	for(std::size_t i = 0; i < dataSize; i++){
		res.push_back(pointer[i]);
	}
	return res;
}

/*!
 * @brief Returns the vector of values which correspond to the pointer parameter.
 * @param pointer The size_t pointer to convert. 
 * @return The vector which contains the pointer's values, with double type. 
 */
std::vector<double> translateAndConvertPointer(std::size_t* pointer, std::size_t dataSize ){
	std::vector<double> res;
	for(std::size_t i = 0; i < dataSize; i++){
		res.push_back((double)pointer[i]);
	}
	return res;
}

std::vector<std::vector<size_t>> hungarianLSAP(std::vector<std::vector<std::size_t>> matrixCost){
	std::size_t nrows = matrixCost.size();
	std::size_t ncols = matrixCost[0].size();
	std::size_t *rho = new std::size_t[nrows], *varrho = new std::size_t[ncols];
	std::size_t *u = new std::size_t[nrows], *v = new std::size_t[ncols];
	std::size_t *C = new std::size_t[nrows*ncols];
	std::size_t i = 0, j;
	for (std::size_t i = 0; i < nrows; i++){
		for (std::size_t j = 0; j < ncols; j++){
			C[j*nrows+i] = matrixCost[i][j];
		}
	}
	lsape::hungarianLSAP<std::size_t>(C,nrows,ncols,rho,u,v,varrho);
	std::vector<std::vector<size_t>> res;
	res.push_back(translatePointer(rho, nrows));
	res.push_back(translatePointer(varrho, ncols));
	res.push_back(translatePointer(u, nrows));
	res.push_back(translatePointer(v, ncols));
	return res;
}

std::vector<std::vector<double>> hungarianLSAPE(std::vector<std::vector<double>> matrixCost){
	std::size_t nrows = matrixCost.size();
	std::size_t ncols = matrixCost[0].size();
	std::size_t *rho = new std::size_t[nrows-1], *varrho = new std::size_t[ncols-1];
	double *u = new double[nrows], *v = new double[ncols];
	double *C = new double[nrows*ncols];
	for (std::size_t i = 0; i < nrows; i++){
		for (std::size_t j = 0; j < ncols; j++){
			C[j*nrows+i] = matrixCost[i][j];
		}
	}
	lsape::hungarianLSAPE<double,std::size_t>(C,nrows,ncols,rho,varrho,u,v);
	std::vector<std::vector<double>> res;
	res.push_back(translateAndConvertPointer(rho, nrows-1));
	res.push_back(translateAndConvertPointer(varrho, ncols-1));
	res.push_back(translatePointer(u, nrows));
	res.push_back(translatePointer(v, ncols));
	return res;
}

/*void medianLetter(pathFolder, pathXML, editCost, method, options="", initOption = "EAGER_WITHOUT_SHUFFLED_COPIES"){
	
	if(isInitialized()){
		restartEnv();
	}
	setEditCost(editCost);*/
	
	/*std::string letter_class("A");
	if (argc > 1) {
		letter_class = std::string(argv[1]);
	}*/
	//std::string seed("0");
	/*if (argc > 2) {
		seed = std::string(argv[2]);
	}*/
	
	/*loadGXLGraph(pathFolder, pathXML);
	std::vector<std::size_t> graph_ids = getAllGraphIds();
	std::size_t median_id = env.add_graph("median", "");
	
	initEnv(initOption);
	
	setMethod(method);
	
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median_estimator(&env, false);
	median_estimator.set_options("--init-type RANDOM --randomness PSEUDO --seed " + seed);
	median_estimator.run(graph_ids, median_id);
	std::string gxl_file_name("../output/gen_median_Letter_HIGH_" + letter_class + ".gxl");
	env.save_as_gxl_graph(median_id, gxl_file_name);*/
	
	/*std::string tikz_file_name("../output/gen_median_Letter_HIGH_" + letter_class + ".tex");
	save_letter_graph_as_tikz_file(env.get_graph(median_id), tikz_file_name);*/
//}

/*!
 * @brief Returns the string which contains all element of a int list. 
 * @param vector The vector to translate. 
 * @return The string which contains all elements separated with a blank space. 
 */
std::string toStringVectorInt(std::vector<int> vector){
	std::string res = "";

    for (int i = 0; i != vector.size(); i++)
    {
       res += std::to_string(vector[i]) + " ";
    }
    
    return res;
}

/*!
 * @brief Returns the string which contains all element of a unsigned long int list. 
 * @param vector The vector to translate. 
 * @return The string which contains all elements separated with a blank space. 
 */
std::string toStringVectorInt(std::vector<unsigned long int> vector){
	std::string res = "";

    for (int i = 0; i != vector.size(); i++)
    {
        res += std::to_string(vector[i]) + " ";
    }
    
    return res;
}
