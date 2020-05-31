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
 * @file  ged_data.hpp
 * @brief ged::GEDData class declaration.
 */

#ifndef SRC_ENV_GED_DATA_HPP_
#define SRC_ENV_GED_DATA_HPP_

#include "ged_graph.hpp"
#include "result.hpp"
#include "matrix.hpp"
#include "common_types.hpp"
#include "node_map.hpp"
#include "../edit_costs/edit_costs.hpp"
#include "../edit_costs/chem_1.hpp"
#include "../edit_costs/chem_2.hpp"
#include "../edit_costs/cmu.hpp"
#include "../edit_costs/grec_1.hpp"
#include "../edit_costs/grec_2.hpp"
#include "../edit_costs/protein.hpp"
#include "../edit_costs/fingerprint.hpp"
#include "../edit_costs/letter.hpp"
#include "../edit_costs/constant.hpp"

namespace ged {

template<class, class, class> class GEDEnv;

/*!
 * @brief Contains the standardized input data along with basic functionality.
 * @details All derived classes of ged::GEDMethod access the input graphs and their edit costs via this class.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class GEDData {

	template<class, class, class> friend class GEDEnv;

public:

	/*!
	 * @brief Destructor.
	 */
	~GEDData();

	/*!
	 * @brief Returns the number of graphs.
	 * @return Number of graphs in the instance.
	 */
	std::size_t num_graphs() const;

	/*!
	 * @brief Provides access to a graph.
	 * @param[in] graph_id The ID of the graph.
	 * @return Constant reference to the graph with ID @p graph_id.
	 */
	const GEDGraph & graph(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Checks if shuffled graph copies are available.
	 * @return Boolean @p true if shuffled graph copies are available.
	 */
	bool shuffled_graph_copies_available() const;

	/*!
	 * @brief Returns the number of graphs in the instance without the shuffled copies.
	 * @return Number of graphs without shuffled copies contained in the instance.
	 */
	std::size_t num_graphs_without_shuffled_copies() const;

	/*!
	 * @brief Returns ID of a graph's shuffled copy.
	 * @param[in] graph_id The ID of the graph.
	 * @return The ID of the shuffled copy of the graph with ID @p graph_id.
	 */
	GEDGraph::GraphID id_shuffled_graph_copy(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Checks if a graph is a shuffled copy of another graph.
	 * @param[in] graph_id The ID of the graph.
	 * @return Boolean @p true if the graph with ID @p graph_id is the shuffled copy of another graph.
	 */
	bool is_shuffled_graph_copy(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Provides access to the graphs contained in the instance.
	 * @return Constant iterator pointing to the first graph in the instance.
	 * @note The member functions begin() and end() provide access to both the original graphs and the shuffled copies.
	 * If necessary, use the member function is_shuffled_graph_copy() to check if a graph is a shuffle copy of another graph.
	 */
	std::vector<GEDGraph>::const_iterator begin() const;

	/*!
	 * @brief Provides access to the graphs contained in the instance.
	 * @return Constant iterator pointing to one position after the last graph in the instance.
	 * @note The member functions begin() and end() provide access to both the original graphs and the shuffled copies.
	 * If necessary, use the member function is_shuffled_graph_copy() to check if a graph is a shuffle copy of another graph.
	 */
	std::vector<GEDGraph>::const_iterator end() const;

	/*!
	 * @brief Returns the number of node labels.
	 * @return Number of different node labels contained in the instance.
	 */
	std::size_t num_node_labels() const;

	/*!
	 * @brief Returns the number of edge labels.
	 * @return Number of different edge labels contained in the instance.
	 */
	std::size_t num_edge_labels() const;

	/*!
	 * @brief Returns maximal number of nodes.
	 * @return Maximal number of nodes contained in the instance's graphs.
	 */
	std::size_t max_num_nodes() const;

	/*!
	 * @brief Returns maximal number of nodes.
	 * @return Maximal number of nodes contained in the instance's graphs.
	 */
	std::size_t max_num_edges() const;

	/*!
	 * @brief Returns node relabeling, insertion, or deletion cost.
	 * @param[in] label1 First node label.
	 * @param[in] label2 Second node label.
	 * @return Node relabeling cost if @p label1 and @p label2 are both different from ged::dummy_label(),
	 * node insertion cost if @p label1 equals ged::dummy_label and @p label2 does not,
	 * node deletion cost if @p label1 does not equal ged::dummy_label and @p label2 does,
	 * and 0 otherwise.
	 */
	double node_cost(LabelID label1, LabelID label2) const;

	/*!
	 * @brief Computes a node label's representation as a real-valued vector.
	 * @param[in] node_label A node label.
	 * @param[out] vector_representation The node label's vector representation.
	 * @note If the selected edit costs do not override ged::EditCosts::vectorize_node_label(), @p vector_representation is empty.
	 */
	void vectorize_node_label(LabelID node_label, std::vector<double> & vector_representation) const;

	/*!
	 * @brief Returns edge relabeling, insertion, or deletion cost.
	 * @param[in] label1 First edge label.
	 * @param[in] label2 Second edge label.
	 * @return Edge relabeling cost if @p label1 and @p label2 are both different from ged::dummy_label(),
	 * edge insertion cost if @p label1 equals ged::dummy_label and @p label2 does not,
	 * edge deletion cost if @p label1 does not equal ged::dummy_label and @p label2 does,
	 * and 0 otherwise.
	 */
	double edge_cost(LabelID label1, LabelID label2) const;

	/*!
	 * @brief Computes an edge label's representation as a real-valued vector.
	 * @param[in] edge_label An edge label.
	 * @param[out] vector_representation The edge label's vector representation.
	 * @note If the selected edit costs do not override ged::EditCosts::vectorize_edge_label(), @p vector_representation is empty.
	 */
	void vectorize_edge_label(LabelID edge_label, std::vector<double> & vector_representation) const;

	/*!
	 * @brief Saves a node map.
	 * @param[in] filename Name of the file where the node map should be saved.
	 * @param[in] g_id ID of a graph.
	 * @param[in] h_id ID of a graph.
	 * @param[in] node_map Node map between the graphs with IDs @p g_id and @p h_id.
	 * @param[in] append Boolean @p true if the node map should be appended to the file and @p false if the file should be overridden.
	 */
	void save_node_map(const std::string & filename, GEDGraph::NodeID g_id, GEDGraph::NodeID h_id, const NodeMap & node_map, bool append = true) const;

	/*!
	 * @brief Loads a node map from a file.
	 * @param[in] filename Name of a file which contains a node map between the graphs with IDs @p g_id and @p h_id.
	 * @param[in] g_id ID of a graph.
	 * @param[in] h_id ID of a graph.
	 * @param[out] node_map Node map between the graphs with IDs @p g_id and @p h_id populated by the method.
	 * @note If the file contains more than one node map between the graphs with IDs @p g_id and @p h_id, the first node map is loaded.
	 */
	void load_node_map(const std::string & filename, GEDGraph::NodeID g_id, GEDGraph::NodeID h_id, NodeMap & node_map) const;

	/*!
	 * @brief Computes the edit cost between two graphs induced by a node map.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in,out] node_map Node map whose induced edit cost is to be computed.
	 */
	void compute_induced_cost(const GEDGraph & g, const GEDGraph & h, NodeMap & node_map) const;

	/*!
	 * @brief Computes the cost of swapping two assignments in a node map while leaving the node map unchanged.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in] assignment_1 First assignment.
	 * @param[in] assignment_2 Second assignment.
	 * @param[in] node_map The node map. Must contain the @p assignment_1 and @p assignment_2 in its relational representation.
	 * @return The cost of swapping @p assignment_1 and @p assignment_2 in @p node_map.
	 */
	double swap_cost(const GEDGraph & g, const GEDGraph & h, const NodeMap::Assignment & assignment_1, const NodeMap::Assignment & assignment_2, NodeMap & node_map) const;

	/*!
	 * @brief Swaps two assignments in a node map.
	 * @param assignment_1 First assignment.
	 * @param assignment_2 Second assignment.
	 * @param swap_cost The cost of swapping @p assignment_1 and @p assignment_2 in @p node_map. Can be computed by call to swap_cost().
	 * @param node_map The node map. Must contain the @p assignment_1 and @p assignment_2 in its relational representation.
	 */
	void swap_assignments(const NodeMap::Assignment & assignment_1, const NodeMap::Assignment & assignment_2, double swap_cost, NodeMap & node_map) const;

	/*!
	 * @brief Checks if the edit costs are quasimetric.
	 * @return Boolean @p true if the edit costs are quasimetric and @p false otherwise.
	 */
	bool quasimetric_costs() const;

	/*!
	 * @brief Checks if the edit costs between two graphs are quasimetric.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Boolean @p true if the edit costs between the graphs @p g and @p h quasimetric and @p false otherwise.
	 */
	bool quasimetric_costs(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the maximal node edit cost between any two graphs contained in the instance.
	 * @return Maximal node edit cost.
	 */
	double max_node_edit_cost() const;

	/*!
	 * @brief Returns the maximal edge edit cost between any two graphs contained in the instance.
	 * @return Maximal edge edit cost.
	 */
	double max_edge_edit_cost() const;

	/*!
	 * @brief Returns the maximal edit cost between two graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Maximal edit cost between the graphs @p g and @p h.
	 */
	double max_edit_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the minimal edit cost between two graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Minimal edit cost between the graphs @p g and @p h.
	 */
	double min_edit_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the maximal node edit cost between two graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Maximal node edit cost between the graphs @p g and @p h.
	 */
	double max_node_edit_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the minimal node edit cost between two graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Minimal node edit cost between the graphs @p g and @p h.
	 */
	double min_node_edit_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the maximal edge edit cost between two graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Maximal node edit cost between the graphs @p g and @p h.
	 */
	double max_edge_edit_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the minimal edge edit cost between two graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Minimal edge edit cost between the graphs @p g and @p h.
	 */
	double min_edge_edit_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the maximal cost of deleting a node contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Maximal cost of deleting a node contained in @p graph.
	 */
	double max_node_del_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the minimal cost of deleting a node contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Minimal cost of deleting a node contained in @p graph.
	 */
	double min_node_del_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the mean cost of deleting a node contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Mean cost of deleting a node contained in @p graph.
	 */
	double mean_node_del_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the maximal cost of inserting a node contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Maximal cost of inserting a node contained in @p graph.
	 */
	double max_node_ins_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the minimal cost of inserting a node contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Minimal cost of inserting a node contained in @p graph.
	 */
	double min_node_ins_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the mean cost of inserting a node contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Mean cost of inserting a node contained in @p graph.
	 */
	double mean_node_ins_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the maximal cost of substituting a node contained in a graph by a node contained in another graph.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Maximal cost of substituting a node contained in @p g by a node contained in @p h.
	 */
	double max_node_subs_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the minimal cost of substituting a node contained in a graph by a node contained in another graph.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Minimal cost of substituting a node contained in @p g by a node contained in @p h.
	 */
	double min_node_subs_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the mean cost of substituting a node contained in a graph by a node contained in another graph.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Mean cost of substituting a node contained in @p g by a node contained in @p h.
	 */
	double mean_node_subs_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the maximal cost of deleting a edge contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Maximal cost of deleting a edge contained in @p graph.
	 */
	double max_edge_del_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the minimal cost of deleting a edge contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Minimal cost of deleting a edge contained in @p graph.
	 */
	double min_edge_del_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the mean cost of deleting a edge contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Mean cost of deleting a edge contained in @p graph.
	 */
	double mean_edge_del_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the maximal cost of inserting a edge contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Maximal cost of inserting a edge contained in @p graph.
	 */
	double max_edge_ins_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the minimal cost of inserting a edge contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Minimal cost of inserting a edge contained in @p graph.
	 */
	double min_edge_ins_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the mean cost of inserting a edge contained in a graph.
	 * @param[in] graph Input graph.
	 * @return Mean cost of inserting a edge contained in @p graph.
	 */
	double mean_edge_ins_cost(const GEDGraph & graph) const;

	/*!
	 * @brief Returns the maximal cost of substituting a edge contained in a graph by a edge contained in another graph.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Maximal cost of substituting a edge contained in @p g by a edge contained in @p h.
	 */
	double max_edge_subs_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the minimal cost of substituting a edge contained in a graph by a edge contained in another graph.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Minimal cost of substituting a edge contained in @p g by a edge contained in @p h.
	 */
	double min_edge_subs_cost(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns the mean cost of substituting a edge contained in a graph by a edge contained in another graph.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @return Mean cost of substituting a edge contained in @p g by a edge contained in @p h.
	 */
	double mean_edge_subs_cost(const GEDGraph & g, const GEDGraph & h) const;

private:

	std::vector<GEDGraph> graphs_;

	std::vector<std::string> graph_names_;

	std::vector<std::string> graph_classes_;

	std::size_t num_graphs_without_shuffled_copies_;

	std::vector<std::map<std::string, GEDGraph::NodeID>> strings_to_internal_node_ids_;

	std::vector<std::map<GEDGraph::NodeID, std::string>> internal_node_ids_to_strings_;

	EditCosts<UserNodeLabel, UserEdgeLabel> * edit_costs_;

	DMatrix node_costs_;

	DMatrix edge_costs_;

    std::vector<UserNodeLabel> node_labels_;

    std::map<UserNodeLabel, LabelID> node_label_ids_;

    std::vector<UserEdgeLabel> edge_labels_;

    std::map<UserEdgeLabel, LabelID> edge_label_ids_;

	Options::InitType init_type_;

	bool delete_edit_costs_;

	std::size_t max_num_nodes_;

	std::size_t max_num_edges_;

	GEDData();

	void set_edit_costs_(Options::EditCosts edit_costs, const std::vector<double> & edit_cost_constants);

	void set_edit_costs_(EditCosts<UserNodeLabel, UserEdgeLabel> * edit_costs);

	void init_cost_matrices_(bool print_to_stdout = false);

	bool eager_init_() const;

	LabelID node_label_to_id_(const UserNodeLabel & node_label);

	UserNodeLabel id_to_node_label(LabelID label_id) const;

	LabelID edge_label_to_id_(const UserEdgeLabel & edge_label);

	UserEdgeLabel id_to_edge_label(LabelID label_id) const;

	GEDGraph::NodeID string_to_node_id_(GEDGraph::GraphID graph_id, const std::string & string) const;

	std::string node_id_to_string_(GEDGraph::GraphID graph_id, GEDGraph::NodeID node_id) const;

};

}

#include "ged_data.ipp"

#endif /* SRC_ENV_GED_DATA_HPP_ */
