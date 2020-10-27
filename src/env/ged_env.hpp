/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2018 by David B. Blumenthal                              *
 *                                                                          *
 *   This file is part of GEDLIB.                                           *
 *                                                                          *
 *   GEDLIB is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 3 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   GEDLIB is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                          *
 ***************************************************************************/

/*!
 * @file  ged_env.hpp
 * @brief ged::GEDEnv class declaration.
 */

#ifndef SRC_ENV_GED_ENV_HPP_
#define SRC_ENV_GED_ENV_HPP_

#include "common_types.hpp"
#include "ged_graph.hpp"
#include "node_map.hpp"
#include "../methods/all_methods.hpp"
#include "ged_data.hpp"

/*!
 * @namespace ged
 * @brief Global namespace for GEDLIB.
 */
namespace ged {

/*!
 * @brief Provides the API of GEDLIB.
 * @tparam UserNodeID Class of user-specific node IDs.
 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
 */
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
class GEDEnv {
public:

	/*!
	 * @brief Destructor.
	 */
	~GEDEnv();

	/*!
	 * @brief Constructor.
	 */
	GEDEnv();

	/*!
	 * @brief Sets the edit costs to one of the predefined edit costs.
	 * @param[in] edit_costs Select one of the predefined edit costs.
	 * @param[in] edit_cost_constants Constants passed to the constructor of the edit cost class selected by @p edit_costs.
	 */
	void set_edit_costs(Options::EditCosts edit_costs, std::vector<double> edit_cost_constants = {});

	/*!
	 * @brief Sets the edit costs to user defined edit costs.
	 * @param[in] edit_costs Pointer to user defined edit costs. Must be freed by the user.
	 */
	void set_edit_costs(EditCosts<UserNodeLabel, UserEdgeLabel> * edit_costs);

	/*!
	 * @brief Adds a new uninitialized graph to the environment. Call init() after calling this method.
	 * @param[in] graph_name The name of the added graph. Empty if not specified.
	 * @param[in] graph_class The class of the added graph. Empty if not specified.
	 * @return The ID of the newly added graph.
	 */
	GEDGraph::GraphID add_graph(const std::string & graph_name = "", const std::string & graph_class = "");

	/*!
	 * @brief Clears and de-initializes a graph that has previously been added to the environment. Call init() after calling this method.
	 * @param[in] graph_id ID of graph that has to be cleared.
	 */
	void clear_graph(GEDGraph::GraphID graph_id);

	/*!
	 * @brief Adds a labeled node.
	 * @param[in] graph_id ID of graph that has been added to the environment.
	 * @param[in] node_id The user-specific ID of the vertex that has to be added.
	 * @param[in] node_label The label of the vertex that has to be added. Set to ged::NoLabel() if template parameter @p UserNodeLabel equals ged::NoLabel.
	 */
	void add_node(GEDGraph::GraphID graph_id, const UserNodeID & node_id, const UserNodeLabel & node_label);

	/*!
	 * @brief Adds a labeled edge.
	 * @param[in] graph_id ID of graph that has been added to the environment.
	 * @param[in] tail The user-specific ID of the tail of the edge that has to be added.
	 * @param[in] head The user-specific ID of the head of the edge that has to be added.
	 * @param[in] edge_label The label of the vertex that has to be added. Set to ged::NoLabel() if template parameter @p UserEdgeLabel equals ged::NoLabel.
	 * @param[in] ignore_duplicates If @p true, duplicate edges are ignores. Otherwise, an exception is thrown if an existing edge is added to the graph.
	 */
	void add_edge(GEDGraph::GraphID graph_id, const UserNodeID & tail, const UserNodeID & head, const UserEdgeLabel & edge_label, bool ignore_duplicates = true);

	/*!
	 * @brief Loads ged::ExchangeGraph into the environment.
	 * @param[in] exchange_graph The ged::ExchangeGraph graph that should be loaded.
	 * @param[in] graph_id The ID of a graph contained the environment (overwrite existing graph) or ged::undefined() (add new graph).
	 * @param[in] graph_name The name of newly added graph. Empty if not specified. Has no effect unless @p graph_id equals ged::undefined().
	 * @param[in] graph_class The class of newly added graph. Empty if not specified. Has no effect unless @p graph_id equals ged::undefined().
	 * @return The ID of the newly loaded graph.
	 * @note If @p graph_id is set to ged::undefined(), all graphs contained in the environment have to be re-initialized.
	 */
	GEDGraph::GraphID load_exchange_graph(const ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & exchange_graph,
			GEDGraph::GraphID graph_id = ged::undefined(), Options::ExchangeGraphType exchange_graph_type = Options::ExchangeGraphType::ADJ_MATRIX,
			const std::string & graph_name = "", const std::string & graph_class = "");

	/*!
	 * @brief Load graph given in the [GXL file format](http://www.gupro.de/GXL/).
	 * @param[in] graph_dir The path to the directory containing the graphs.
	 * @param[in] gxl_file_name Name of the GXL file.
	 * @param[in] node_type Select if nodes are labeled or unlabeled.
	 * @param[in] edge_type Select if edges are labeled or unlabeled.
	 * @param[in] irrelevant_node_attributes Set of node attributes that are irrelevant for the selected edit costs.
	 * @param[in] irrelevant_edge_attributes Set of edge attributes that are irrelevant for the selected edit costs.
	 * @param[in] graph_id The ID of a graph contained the environment (overwrite existing graph) or ged::undefined() (add new graph).
	 * @param[in] graph_class The class of the added graph. Empty if not specified. Has no effect unless @p graph_id equals ged::undefined().
	 * @return The ID of the newly loaded graph.
	 * @note If @p graph_id is set to ged::undefined(), all graphs contained in the environment have to be re-initialized.
	 * @warning Calls to this method create a compiler error unless the template parameters @p UserNodeID is set to ged::GXLUserNodeID
	 * and the template parameters @p UserNodeLabel and @p UserEdgeLabel are set to ged::GXLLabel.
	 */
	GEDGraph::GraphID load_gxl_graph(const std::string & graph_dir, const std::string & gxl_file_name, Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
			const std::unordered_set<std::string> & irrelevant_node_attributes, const std::unordered_set<std::string> & irrelevant_edge_attributes,
			GEDGraph::GraphID graph_id = ged::undefined(), const std::string & graph_class = "");

	/*!
	 * @brief Loads graphs given in the [GXL file format](http://www.gupro.de/GXL/).
	 * @param[in] graph_dir The path to the directory containing the graphs.
	 * @param[in] collection_file The path to a XML file thats lists the graphs contained in @p graph_dir that should be loaded.
	 * @param[in] node_type Select if nodes are labeled or unlabeled.
	 * @param[in] edge_type Select if edges are labeled or unlabeled.
	 * @param[in] irrelevant_node_attributes Set of node attributes that are irrelevant for the selected edit costs.
	 * @param[in] irrelevant_edge_attributes Set of edge attributes that are irrelevant for the selected edit costs.
	 * @return A vector containing the IDs of the newly added graphs.
	 * @warning Calls to this method create a compiler error unless the template parameters @p UserNodeID is set to ged::GXLUserNodeID
	 * and the template parameters @p UserNodeLabel and @p UserEdgeLabel are set to ged::GXLLabel.
	 */
	std::vector<GEDGraph::GraphID> load_gxl_graphs(const std::string & graph_dir, const std::string & collection_file,
			Options::GXLNodeEdgeType node_type = Options::GXLNodeEdgeType::LABELED, Options::GXLNodeEdgeType edge_type = Options::GXLNodeEdgeType::LABELED,
			const std::unordered_set<std::string> & irrelevant_node_attributes = std::unordered_set<std::string>(), const std::unordered_set<std::string> & irrelevant_edge_attributes = std::unordered_set<std::string>());

	/*!
	 * @brief Saves graph contained in the environment in the [GXL file format](http://www.gupro.de/GXL/).
	 * @param[in] graph_id ID of the graph that should be saved as GXL file.
	 * @param[in,out] gxl_file_name Name of GXL file.
	 * @warning Calls to this method create a compiler error unless the template parameters @p UserNodeID is set to ged::GXLUserNodeID
	 * and the template parameters @p UserNodeLabel and @p UserEdgeLabel are set to ged::GXLLabel.
	 */
	void save_as_gxl_graph(GEDGraph::GraphID graph_id, const std::string & gxl_file_name) const;

	/*!
	 * @brief Creates an XML GraphCollection file.
	 * @param[in] xml_file_name Name of the created XML file.
	 * @param[in] gxl_file_names Vector that contains the names of the GXL files that contain the graphs listed in the collection.
	 * @param[in] graph_classes Vector that contains the classes of the graphs listed in the collection.
	 * If left empty, all graphs are assigned the class "no_class". Otherwise, it must have the same size as @p gxl_file_names.
	 */
	void save_graph_collection(const std::string & xml_file_name, const std::vector<std::string> & gxl_file_names, const std::vector<std::string> & graph_classes = {}) const;

	/*!
	 * @brief Initializes the environment.
	 * @param[in] init_type Select initialization type.
	 * @param[in] print_to_stdout If set to @p true, the progress of the initialization is printed to std::out.
	 */
	void init(Options::InitType init_type = Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES, bool print_to_stdout = false);

	/*!
	 * @brief Check if the environment is initialized.
	 * @return True if the environment is initialized.
	 */
	bool initialized() const;

	/*!
	 * @brief Returns the initialization type of the last initialization.
	 * @return Initialization type.
	 */
	Options::InitType get_init_type() const;

	/*!
	 * @brief Sets the GEDMethod to be used by run_method().
	 * @param[in] method Select the method that is to be used.
	 * @param[in] options An options string of the form @"[--@<option@> @<arg@>] [...]@" passed to the selected method.
	 */
	void set_method(Options::GEDMethod method, const std::string & options = std::string(""));

	/*!
	 * @brief Runs the GED method specified by call to set_method() between the graphs with IDs @p g_id and @p h_id.
	 * @param[in] g_id ID of an input graph that has been added to the environment.
	 * @param[in] h_id ID of an input graph that has been added to the environment.
	 * @param[in] use_shuffled_graphs If true, the shuffled copies of the graphs are used if available.
	 */
	void run_method(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id, bool use_shuffled_graphs = false);

	/*!
	 * @brief Initializes the method specified by call to set_method().
	 */
	void init_method();

	/*!
	 * @brief Provides access to the IDs of the graphs contained in the environment.
	 * @return Pair <tt>(ID of first graphs, ID of last graph + 1)</tt> of graph IDs.
	 * If both entries equal 0, the environment does not contain any graphs.
	 */
	std::pair<GEDGraph::GraphID, GEDGraph::GraphID> graph_ids() const;

	/*!
	 * @brief The number of graphs contained in the environment.
	 * @return The number of graphs (without shuffled copies).
	 */
	std::size_t num_graphs() const;

	/*!
	 * @brief Returns the number of node labels.
	 * @return Number of pairwise different node labels contained in the environment.
	 * @note If @p 1 is returned, the nodes are unlabeled.
	 */
	std::size_t num_node_labels() const;

	/*!
	 * @brief Returns node label.
	 * @param[in] label_id ID of node label that should be returned. Must be between 1 and num_node_labels().
	 * @return Node label for selected label ID.
	 */
	UserNodeLabel get_node_label(LabelID label_id) const;

	/*!
	 * @brief Returns the number of edge labels.
	 * @return Number of pairwise different edge labels contained in the environment.
	 * @note If @p 1 is returned, the edges are unlabeled.
	 */
	std::size_t num_edge_labels() const;

	/*!
	 * @brief Returns edge label.
	 * @param[in] label_id ID of edge label that should be returned. Must be between 1 and num_node_labels().
	 * @return Edge label for selected label ID.
	 */
	UserEdgeLabel get_edge_label(LabelID label_id) const;

	/*!
	 * @brief Returns lower bound for edit distance between the input graphs.
	 * @param[in] g_id ID of an input graph that has been added to the environment.
	 * @param[in] h_id ID of an input graph that has been added to the environment.
	 * @return Lower bound computed by the last call to run_method() with arguments @p g_id and @p h_id.
	 */
	double get_lower_bound(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const;

	/*!
	 * @brief Returns upper bound for edit distance between the input graphs.
	 * @param[in] g_id ID of an input graph that has been added to the environment.
	 * @param[in] h_id ID of an input graph that has been added to the environment.
	 * @return Upper bound computed by the last call to run_method() with arguments @p g_id and @p h_id.
	 */
	double get_upper_bound(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const;

	/*!
	 * @brief Returns node map between the input graphs.
	 * @param[in] g_id ID of an input graph that has been added to the environment.
	 * @param[in] h_id ID of an input graph that has been added to the environment.
	 * @return Node map computed by the last call to run_method() with arguments @p g_id and @p h_id.
	 */
	const NodeMap & get_node_map(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const;

	/*!
	 * @brief Returns runtime.
	 * @param[in] g_id ID of an input graph that has been added to the environment.
	 * @param[in] h_id ID of an input graph that has been added to the environment.
	 * @return Runtime of last call to run_method() with arguments @p g_id and @p h_id.
	 */
	double get_runtime(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const;

	/*!
	 * @brief Returns initialization time.
	 * @return Runtime of the last call to init_method().
	 */
	double get_init_time() const;

	/*!
	 * @brief Computes the edit cost between two graphs induced by a node map.
	 * @param[in] g_id ID of input graph.
	 * @param[in] h_id ID of input graph.
	 * @param[in,out] node_map Node map whose induced edit cost is to be computed.
	 */
	void compute_induced_cost(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id, NodeMap & node_map) const;

	/*!
	 * @brief Returns ged::ExchangeGraph representation.
	 * @param[in] graph_id ID of the selected graph.
	 * @param[in] adj_matrix Set to @p true to construct and adjacency matrix @p adj_matrix and a hash-map
	 * @p edge_labels, which has a key for each pair @p (i,j) such that @p adj_matrix[i][j] equals 1.
	 * @param[in] adj_lists Construct adjacency lists for all nodes, i.e., a vector @p adj_lists that, for each node
	 * @p i and each edge @p (i,j), contains a vector with entries of the form @p (j,edge_label).
	 * @param[in] edge_list Construct a list of all edges that contains entries of the form @p ((i,j),edge_label).
	 * @return ged::ExchangeGraph representation of the selected graph.
	 */
	ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> get_graph(GEDGraph::GraphID graph_id, bool adj_matrix = true,
			bool adj_lists = false, bool edge_list = false) const;

	/*!
	 * @brief Returns the graph class.
	 * @param[in] graph_id ID of an input graph that has been added to the environment.
	 * @return Class of the input graph.
	 */
	const std::string & get_graph_class(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Returns the graph name.
	 * @param[in] graph_id ID of an input graph that has been added to the environment.
	 * @return Name of the input graph.
	 */
	const std::string & get_graph_name(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Checks if the edit costs are quasimetric.
	 * @return Boolean @p true if the edit costs are quasimetric and @p false, otherwise.
	 */
	bool quasimetric_costs() const;

	/*!
	 * @brief Returns the number of nodes.
	 * @param[in] graph_id ID of an input graph that has been added to the environment.
	 * @return Number of nodes in the graph.
	 */
	std::size_t get_num_nodes(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Returns average number of nodes.
	 * @return Average number of nodes of the graphs contained in the environment.
	 */
	double get_avg_num_nodes() const;

	/*!
	 * @brief Returns node relabeling cost.
	 * @param[in] node_label_1 First node label.
	 * @param[in] node_label_2 Second node label.
	 * @return Node relabeling cost for the given node labels.
	 */
	double node_rel_cost(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const;

	/*!
	 * @brief Returns node deletion cost.
	 * @param[in] node_label Node label.
	 * @return Cost of deleting node with given label.
	 */
	double node_del_cost(const UserNodeLabel & node_label) const;

	/*!
	 * @brief Returns node insertion cost.
	 * @param[in] node_label Node label.
	 * @return Cost of inserting node with given label.
	 */
	double node_ins_cost(const UserNodeLabel & node_label) const;

	/*!
	 * @brief Computes median node label.
	 * @param[in] node_labels The node labels whose median should be computed.
	 * @return Median of the given node labels.
	 */
	UserNodeLabel median_node_label(const std::vector<UserNodeLabel> & node_labels) const;

	/*!
	 * @brief Returns edge relabeling cost.
	 * @param[in] edge_label_1 First edge label.
	 * @param[in] edge_label_2 Second edge label.
	 * @return Edge relabeling cost for the given edge labels.
	 */
	double edge_rel_cost(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const;

	/*!
	 * @brief Returns edge deletion cost.
	 * @param[in] edge_label Edge label.
	 * @return Cost of deleting edge with given label.
	 */
	double edge_del_cost(const UserEdgeLabel & edge_label) const;

	/*!
	 * @brief Returns edge insertion cost.
	 * @param[in] edge_label Edge label.
	 * @return Cost of inserting edge with given label.
	 */
	double edge_ins_cost(const UserEdgeLabel & edge_label) const;

	/*!
	 * @brief Computes median edge label.
	 * @param[in] edge_labels The edge labels whose median should be computed.
	 * @return Median of the given edge labels.
	 */
	UserEdgeLabel median_edge_label(const std::vector<UserEdgeLabel> & edge_labels) const;

private:

	bool initialized_;

	std::vector<GEDGraph::GraphID> new_graph_ids_;

	GEDData<UserNodeLabel, UserEdgeLabel> ged_data_;

	// Variables needed for approximating ged_instance_.

	std::map<std::pair<GEDGraph::GraphID, GEDGraph::GraphID>, double> lower_bounds_;

	std::map<std::pair<GEDGraph::GraphID, GEDGraph::GraphID>, double> upper_bounds_;

	std::map<std::pair<GEDGraph::GraphID, GEDGraph::GraphID>, Seconds> runtimes_;

	std::map<std::pair<GEDGraph::GraphID, GEDGraph::GraphID>, NodeMap> node_maps_;

	std::vector<std::map<UserNodeID, GEDGraph::NodeID>> original_to_internal_node_ids_;

	std::vector<std::map<GEDGraph::NodeID, UserNodeID>> internal_to_original_node_ids_;

	GEDMethod<UserNodeLabel, UserEdgeLabel> * ged_method_;

	void read_gxl_label_from_ptree_(const boost::property_tree::ptree::value_type & node_or_edge, const std::unordered_set<std::string> & irrelevant_attributes, const std::string & file, GXLLabel & label);

	void construct_shuffled_graph_copy_(GEDGraph::GraphID graph_id);

	GEDGraph::GraphID add_or_clear_shuffled_graph_copy_(GEDGraph::GraphID graph_id);

	std::string to_string_(UserNodeID node_id) const;

	std::string gxl_label_to_string_(GXLLabel gxl_label) const;

};

}

#include "ged_env.ipp"

namespace ged {

#ifdef GXL_GEDLIB_SHARED
#ifndef SRC_ENV_GED_ENV_GXL_CPP_
extern template class GEDEnv<GXLNodeID, GXLLabel, GXLLabel>;
#endif /* SRC_ENV_GED_ENV_GXL_CPP_ */
#endif /* GXL_GEDLIB_SHARED */

}

#endif /* SRC_ENV_GED_ENV_HPP_ */
