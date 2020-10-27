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
 * @file  ged_env.ipp
 * @brief ged::GEDEnv class definition.
 */

#ifndef SRC_ENV_GED_ENV_IPP_
#define SRC_ENV_GED_ENV_IPP_

namespace ged {

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
~GEDEnv() {
	delete ged_method_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
GEDEnv():
initialized_{false},
new_graph_ids_(),
ged_data_(),
lower_bounds_(),
upper_bounds_(),
runtimes_(),
node_maps_(),
original_to_internal_node_ids_(),
internal_to_original_node_ids_(),
ged_method_{nullptr} {}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_edit_costs(Options::EditCosts edit_costs, std::vector<double> edit_cost_constants) {
	ged_data_.set_edit_costs_(edit_costs, edit_cost_constants);
	if (ged_data_.eager_init_()) {
		initialized_ = false;
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_edit_costs(EditCosts<UserNodeLabel, UserEdgeLabel> * edit_costs) {
	ged_data_.set_edit_costs_(edit_costs);
	if (ged_data_.eager_init_()) {
		initialized_ = false;
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDGraph::GraphID
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_graph(const std::string & name, const std::string & graph_class) {
	for (auto & graph : ged_data_.graphs_) {
		graph.un_init();
	}
	initialized_ = false;
	GEDGraph::GraphID graph_id{ged_data_.num_graphs_without_shuffled_copies_++};
	new_graph_ids_.push_back(graph_id);
	ged_data_.graphs_.emplace(ged_data_.graphs_.begin() + graph_id, graph_id);
	ged_data_.graph_names_.insert(ged_data_.graph_names_.begin() + graph_id, name);
	ged_data_.graph_classes_.insert(ged_data_.graph_classes_.begin() + graph_id, graph_class);
	original_to_internal_node_ids_.insert(original_to_internal_node_ids_.begin() + graph_id, std::map<UserNodeID, GEDGraph::NodeID>());
	internal_to_original_node_ids_.insert(internal_to_original_node_ids_.begin() + graph_id, std::map<GEDGraph::NodeID, UserNodeID>());
	ged_data_.strings_to_internal_node_ids_.insert(ged_data_.strings_to_internal_node_ids_.begin() + graph_id, std::map<std::string, GEDGraph::NodeID>());
	ged_data_.internal_node_ids_to_strings_.insert(ged_data_.internal_node_ids_to_strings_.begin() + graph_id, std::map<GEDGraph::NodeID, std::string>());
	return graph_id;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
clear_graph(GEDGraph::GraphID graph_id) {
	if (graph_id >= ged_data_.num_graphs_without_shuffled_copies()) {
		throw Error("The graph " + get_graph_name(graph_id) + " has not been added to the environment.");
	}
	ged_data_.graphs_[graph_id].clear();
	original_to_internal_node_ids_[graph_id].clear();
	internal_to_original_node_ids_[graph_id].clear();
	ged_data_.strings_to_internal_node_ids_[graph_id].clear();
	ged_data_.internal_node_ids_to_strings_[graph_id].clear();
	initialized_ = false;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDGraph::GraphID
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
load_exchange_graph(const ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & exchange_graph, GEDGraph::GraphID graph_id, Options::ExchangeGraphType exchange_graph_type, const std::string & graph_name, const std::string & graph_class) {
	if (graph_id == ged::undefined()) {
		graph_id = add_graph(graph_name, graph_class);
	}
	else {
		clear_graph(graph_id);
	}
	for (GEDGraph::NodeID node_id{0}; node_id < exchange_graph.num_nodes; node_id++) {
		add_node(graph_id, exchange_graph.original_node_ids.at(node_id), exchange_graph.node_labels.at(node_id));
	}
	if (exchange_graph_type == Options::ExchangeGraphType::ADJ_MATRIX) {
		for (GEDGraph::NodeID i{0}; i < exchange_graph.num_nodes; i++) {
			for (GEDGraph::NodeID j{i + 1}; j < exchange_graph.num_nodes; j++) {
				if (exchange_graph.adj_matrix.at(i).at(j) == 1) {
					add_edge(graph_id, exchange_graph.original_node_ids.at(i), exchange_graph.original_node_ids.at(j), exchange_graph.edge_labels.at(std::make_pair(i, j)));
				}
			}
		}
	}
	else if (exchange_graph_type == Options::ExchangeGraphType::ADJ_LISTS) {
		for (GEDGraph::NodeID i{0}; i < exchange_graph.num_nodes; i++) {
			for (const auto & adj : exchange_graph.adj_lists.at(i)) {
				add_edge(graph_id, exchange_graph.original_node_ids.at(i), exchange_graph.original_node_ids.at(adj.first), adj.second);
			}
		}
	}
	else {
		for (const auto & edge : exchange_graph.edge_list) {
			add_edge(graph_id, exchange_graph.original_node_ids.at(edge.first.first), exchange_graph.original_node_ids.at(edge.first.second), edge.second);
		}
	}
	return graph_id;
}

template<>
GEDGraph::GraphID
GEDEnv<GXLNodeID, GXLLabel, GXLLabel>::
load_gxl_graph(const std::string & graph_dir, const std::string & gxl_file_name, Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
		const std::unordered_set<std::string> & irrelevant_node_attributes, const std::unordered_set<std::string> & irrelevant_edge_attributes, GEDGraph::GraphID graph_id, const std::string & graph_class) {

	// read the file into a property tree
	std::string filename("");
	if (graph_dir.back() == '/') {
		filename = graph_dir + gxl_file_name;
	}
	else {
		filename = graph_dir + "/" + gxl_file_name;
	}

	boost::property_tree::ptree root;
	try {
		read_xml(filename, root);
	}
	catch (const boost::property_tree::xml_parser_error & error) {
		throw Error(std::string("Error reading file ") + filename + ": " + error.message() + ".");
	}

	// first sanity checks
	if (root.count("gxl") == 0) {
		throw Error("The file " + gxl_file_name + " has the wrong format: no xml-element <gxl>.");
	}
	if (root.count("gxl") >= 2) {
		throw Error("The file " + gxl_file_name + " has the wrong format: more than one xml-element <gxl>.");
	}
	root = root.get_child("gxl");
	if (root.count("graph") == 0) {
		throw Error("The file " + gxl_file_name + " has the wrong format: no xml-element <gxl>.<graph>");
	}
	if (root.count("graph") >= 2) {
		throw Error("The file " + gxl_file_name + " has the wrong format: more than one xml-element <gxl>.<graph>");
	}
	root = root.get_child("graph");

	// add new graph to the environment
	if (graph_id == ged::undefined()) {
		graph_id = add_graph(gxl_file_name, graph_class);
	}
	else {
		clear_graph(graph_id);
	}

	// initialize local variables needed for construction of the graph
	GXLLabel label;
	std::string attr_name;
	std::string attr_val;
	GXLNodeID v_id, from, to;

	// iterate through the property tree to construct the graph
	for (const boost::property_tree::ptree::value_type & node_or_edge : root) {

		// encountered a new vertex that has to be added to the graph
		if (node_or_edge.first == "node") {
			// determine the vertex ID
			try {
				v_id = node_or_edge.second.get<std::string>("<xmlattr>.id");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + filename + " has the wrong format: missing xml-attribute \"id\" in element <gxl>.<graph>.<node>.");
			}
			// read the node label and add the vertex to the graph

			label.clear();
			if (node_type == Options::GXLNodeEdgeType::LABELED) {
				read_gxl_label_from_ptree_(node_or_edge, irrelevant_node_attributes, filename, label);
			}
			add_node(graph_id, v_id, label);
		}

		// encountered a new edge that has to be added to the graph
		else if (node_or_edge.first == "edge") {
			// determine the edge's tail and head
			try {
				from = node_or_edge.second.get<std::string>("<xmlattr>.from");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + filename + " has the wrong format: missing xml-attribute \"from\" in element <gxl>.<graph>.<edge>.");
			}
			try {
				to = node_or_edge.second.get<std::string>("<xmlattr>.to");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + filename + " has the wrong format: missing xml-attribute \"to\" in element <gxl>.<graph>.<edge>.");
			}
			// read the edge label and add the edge to the graph
			label.clear();
			if (edge_type == Options::GXLNodeEdgeType::LABELED) {
				read_gxl_label_from_ptree_(node_or_edge, irrelevant_edge_attributes, filename, label);
			}
			add_edge(graph_id, from, to, label);
		}

		// sanity check
		else if (node_or_edge.first != "<xmlattr>"){
			throw Error("The file " + filename + " has the wrong format: unexpected element <gxl>.<graph>.<" + node_or_edge.first + ">.");
		}
	}
	// return ID of newly constructed graph
	return graph_id;
}

template<>
std::vector<GEDGraph::GraphID>
GEDEnv<GXLNodeID, GXLLabel, GXLLabel>::
load_gxl_graphs(const std::string & graph_dir, const std::string & file, Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
		const std::unordered_set<std::string> & irrelevant_node_attributes, const std::unordered_set<std::string> & irrelevant_edge_attributes) {
	// read the file into a property tree
	boost::property_tree::ptree root;
	try {
		read_xml(file, root);
	}
	catch (const boost::property_tree::xml_parser_error & error) {
		throw Error(std::string("Error reading file ") + file + ": " + error.message() + ".");
	}
	// first sanity checks
	if (root.count("GraphCollection") == 0) {
		throw Error("The file " + file + " has the wrong format: no xml-element <GraphCollection>.");
	}
	if (root.count("GraphCollection") >= 2) {
		throw Error("The file " + file + " has the wrong format: more than one xml-element <GraphCollection>.");
	}
	root = root.get_child("GraphCollection");


	// Read the listed .gxl-files into the environment.
	std::vector<GEDGraph::GraphID> graph_ids;
	std::string gxl_file("");
	std::string graph_class("");
	for (const boost::property_tree::ptree::value_type & val : root) {
		if (val.first == "graph") {
			try {
				gxl_file = val.second.get<std::string>("<xmlattr>.file");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + file + " has the wrong format: missing xml-attribute \"file\" in element <GraphCollection>.<graph>");
			}
			catch (const boost::property_tree::ptree_bad_data & error) {
				throw Error("The file " + file + " has the wrong format: corrupted content in xml-attribute \"file\" of element <GraphCollection>.<graph>");
			}
			try {
				graph_class = val.second.get<std::string>("<xmlattr>.class");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + file + " has the wrong format: missing xml-attribute \"class\" in element <GraphCollection>.<graph>");
			}
			catch (const boost::property_tree::ptree_bad_data & error) {
				throw Error("The file " + file + " has the wrong format: corrupted content in xml-attribute \"class\" of element <GraphCollection>.<graph>");
			}
			graph_ids.push_back(load_gxl_graph(graph_dir, gxl_file, node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, ged::undefined(), graph_class));
		}
		else if (val.first != "<xmlattr>") {
			throw Error("The file " + file + " has the wrong format: unexpected element <GraphCollection>.<" + val.first + ">.");
		}
	}
	return graph_ids;
}

template<>
void
GEDEnv<GXLNodeID, GXLLabel, GXLLabel>::
save_as_gxl_graph(GEDGraph::GraphID graph_id, const std::string & gxl_file_name) const {
	std::ofstream gxl_file(gxl_file_name.c_str());
	gxl_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	gxl_file << "<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n";
	gxl_file << "<gxl xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
	gxl_file << "<graph id=\"" << get_graph_name(graph_id) << "\" edgeids=\"false\" edgemode=\"undirected\">\n";
	const GEDGraph & graph{ged_data_.graphs_.at(graph_id)};

	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		gxl_file << "<node id=\"_" << *node << "\">";
		gxl_file << gxl_label_to_string_(ged_data_.node_labels_.at(graph.get_node_label(*node) - 1));
		gxl_file << "</node>\n";
	}
	for (auto eitr = graph.edges(); eitr.first != eitr.second; eitr.first++) {
		GEDGraph::EdgeID edge(*eitr.first);
		gxl_file << "<edge from=\"_" << graph.tail(edge) << "\" to=\"_" << graph.head(edge) << "\">";
		gxl_file << gxl_label_to_string_(ged_data_.edge_labels_.at(graph.get_edge_label(edge) - 1));
		gxl_file << "</edge>\n";
	}
	gxl_file << "</graph>\n</gxl>\n";
	gxl_file.close();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
save_graph_collection(const std::string & xml_file_name, const std::vector<std::string> & gxl_file_names, const std::vector<std::string> & graph_classes) const {
	if ((not graph_classes.empty()) and (graph_classes.size() != gxl_file_names.size())) {
		throw Error("The number of GXL files does not match the number of graph classes.");
	}
	std::ofstream xml_file(xml_file_name.c_str());
	xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	xml_file << "<!DOCTYPE GraphCollection SYSTEM \"http://www.inf.unibz.it/~blumenthal/dtd/GraphCollection.dtd\">\n";
	xml_file << "<GraphCollection>\n";
	for (std::size_t graph_id{0}; graph_id < gxl_file_names.size(); graph_id++) {
		xml_file << "\t<graph file=\"" << gxl_file_names.at(graph_id) << "\" class=\"" << (graph_classes.empty() ? "no_class" : graph_classes.at(graph_id)) << "\"/>\n";
	}
	xml_file << "</GraphCollection>\n";
	xml_file.close();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
read_gxl_label_from_ptree_(const boost::property_tree::ptree::value_type & node_or_edge, const std::unordered_set<std::string> & irrelevant_attributes, const std::string & file, GXLLabel & label) {
	std::string attr_name, attr_val;
	std::string info(node_or_edge.first);
	// read the attributes into the label
	for (const boost::property_tree::ptree::value_type & attr : node_or_edge.second) {

		if (attr.first == "attr") {
			// determine the name of the attribute
			try {
				attr_name = attr.second.get<std::string>("<xmlattr>.name");
			}
			catch (const boost::property_tree::ptree_bad_path & error) {
				throw Error("The file " + file + " has the wrong format: missing xml-attribute \"name\" in element <gxl>.<graph>.<" + info + ">.<attr>");
			}
			if (irrelevant_attributes.find(attr_name) != irrelevant_attributes.end()) {
				continue;
			}
			// determine the value of the attribute
			LabelID num_attr {0};
			for (const boost::property_tree::ptree::value_type & val : attr.second) {
				if (val.first != "<xmlattr>") {
					// sanity check
					if (++num_attr > 1) {
						throw Error("The file " + file + " has the wrong format: each element <gxl>.<graph>.<" + info + ">.<attr> has to have exactly one child. Actual number of children: " + std::to_string(num_attr));
					}
					try {
						attr_val = val.second.get<std::string>("");
					}
					catch (const boost::property_tree::ptree_bad_path & error) {
						throw Error("The file " + file + " has the wrong format: missing content <gxl>.<graph>.<" + info + ">.<attr>.<[TYPENAME]>.___.");
					}
					catch (const boost::property_tree::ptree_bad_data & error) {
						throw Error("The file " + file + " has the wrong format: corrupted content <gxl>.<graph>.<" + info + ">.<attr>.<[TYPENAME]>.___.");
					}
				}
			}
			// add the attribute to the label
			label[attr_name] = attr_val;
		}
		// sanity check
		else if (attr.first != "<xmlattr>"){
			throw Error("The file " + file + " has the wrong format: unexpected element <gxl>.<graph>.<" + info + ">.<" + attr.first + ">.");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::string
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
to_string_(UserNodeID node_id) const {
	std::stringstream ss;
	ss << node_id;
	return ss.str();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::string
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
gxl_label_to_string_(GXLLabel gxl_label) const {
	std::stringstream ss;
	for (auto const & attr : gxl_label) {
		ss << "<attr name=\"" << attr.first << "\"><TYPENAME>" << attr.second << "</TYPENAME></attr>";
	}
	return ss.str();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_node(GEDGraph::GraphID graph_id, const UserNodeID & node_id, const UserNodeLabel & node_label) {
	if (graph_id >= ged_data_.num_graphs()) {
		throw Error("The graph " + get_graph_name(graph_id) + " with ID " + std::to_string(graph_id) + " has not been added to the environment. The environment contains " + std::to_string(ged_data_.num_graphs_without_shuffled_copies()) + " graphs.");
	}
	if (original_to_internal_node_ids_[graph_id].find(node_id) != original_to_internal_node_ids_[graph_id].end()) {
		throw Error("The node " + to_string_(node_id) + " has already been added to the graph " + std::to_string(graph_id) + ": " + get_graph_name(graph_id) + ".");
	}
	initialized_ = false;
	GEDGraph::NodeID internal_node_id{ged_data_.graphs_[graph_id].add_node()};
	original_to_internal_node_ids_[graph_id][node_id] = internal_node_id;
	internal_to_original_node_ids_[graph_id][internal_node_id] = node_id;
	ged_data_.strings_to_internal_node_ids_[graph_id][to_string_(node_id)] = internal_node_id;
	ged_data_.internal_node_ids_to_strings_[graph_id][internal_node_id] = to_string_(node_id);
	LabelID label_id{ged_data_.node_label_to_id_(node_label)};
	ged_data_.graphs_[graph_id].set_label(original_to_internal_node_ids_[graph_id][node_id], label_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_edge(GEDGraph::GraphID graph_id, const UserNodeID & from, const UserNodeID & to, const UserEdgeLabel & edge_label, bool ignore_duplicates) {
	if (graph_id >= ged_data_.num_graphs()) {
		throw Error("The graph " + get_graph_name(graph_id) + " with ID " + std::to_string(graph_id) + " has not been added to the environment. The environment contains " + std::to_string(ged_data_.num_graphs_without_shuffled_copies()) + " graphs.");
	}
	if (original_to_internal_node_ids_[graph_id].find(from) == original_to_internal_node_ids_[graph_id].end()) {
		throw Error("The node " + to_string_(from) + " does not exist in the graph " + get_graph_name(graph_id) + ".");
	}
	if (original_to_internal_node_ids_[graph_id].find(to) == original_to_internal_node_ids_[graph_id].end()) {
		throw Error("The node " + to_string_(to) + " does not exist in the graph " + get_graph_name(graph_id) + ".");
	}
	initialized_ = false;
	if (ged_data_.graphs_[graph_id].safe_is_edge(original_to_internal_node_ids_[graph_id][from], original_to_internal_node_ids_[graph_id][to])) {
		if (ignore_duplicates) {
			return;
		}
		throw Error("An edge between " + to_string_(from) + " and " + to_string_(to) + " has already been added to the graph " + get_graph_name(graph_id) + ".");
	}
	GEDGraph::EdgeID edge_id{ged_data_.graphs_[graph_id].add_edge(original_to_internal_node_ids_[graph_id][from], original_to_internal_node_ids_[graph_id][to])};
	LabelID label_id{ged_data_.edge_label_to_id_(edge_label)};
	ged_data_.graphs_[graph_id].set_label(edge_id, label_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_method(Options::GEDMethod method, const std::string & options) {
	delete ged_method_;
	switch (method) {
	case Options::GEDMethod::BRANCH:
		ged_method_ = new Branch<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BRANCH_FAST:
		ged_method_ = new BranchFast<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BRANCH_TIGHT:
		ged_method_ = new BranchTight<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BRANCH_UNIFORM:
		ged_method_ = new BranchUniform<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BRANCH_COMPACT:
		ged_method_ = new BranchCompact<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::PARTITION:
		ged_method_ = new Partition<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::HYBRID:
		ged_method_ = new Hybrid<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::RING:
		ged_method_ = new Ring<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::ANCHOR_AWARE_GED:
		ged_method_ = new AnchorAwareGED<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::WALKS:
		ged_method_ = new Walks<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::IPFP:
		ged_method_ = new IPFP<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BIPARTITE:
		ged_method_ = new Bipartite<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::SUBGRAPH:
		ged_method_ = new Subgraph<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::NODE:
		ged_method_ = new Node<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::RING_ML:
		ged_method_ = new RingML<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BIPARTITE_ML:
		ged_method_ = new BipartiteML<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::REFINE:
		ged_method_ = new Refine<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BP_BEAM:
		ged_method_ = new BPBeam<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::SIMULATED_ANNEALING:
		ged_method_ = new SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::HED:
		ged_method_ = new HED<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::STAR:
		ged_method_ = new Star<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
#ifdef GUROBI
	case Options::GEDMethod::F1:
		ged_method_ = new F1<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::F2:
		ged_method_ = new F2<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::COMPACT_MIP:
		ged_method_ = new CompactMIP<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
	case Options::GEDMethod::BLP_NO_EDGE_LABELS:
		ged_method_ = new BLPNoEdgeLabels<UserNodeLabel, UserEdgeLabel>(ged_data_);
		break;
#endif
	}
	ged_method_->set_options(options);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
run_method(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id, bool use_shuffled_graphs) {
	if (g_id >= ged_data_.num_graphs()) {
		throw Error("The graph with ID " + std::to_string(g_id) + " has not been added to the environment.");
	}
	if (h_id >= ged_data_.num_graphs()) {
		throw Error("The graph with ID " + std::to_string(h_id) + " has not been added to the environment.");
	}
	if (not initialized_) {
		throw Error("The environment is uninitialized. Call init() after adding all graphs to the environment.");
	}
	if (not ged_method_) {
		throw Error("No method has been set. Call set_method() before calling run().");
	}
	// call selected GEDMethod and store results
	if (ged_data_.shuffled_graph_copies_available() and (g_id == h_id)) {
		ged_method_->run(g_id, ged_data_.id_shuffled_graph_copy(h_id));
	}
	else if (ged_data_.shuffled_graph_copies_available() and use_shuffled_graphs) {
        ged_method_->run(ged_data_.id_shuffled_graph_copy(g_id), ged_data_.id_shuffled_graph_copy(h_id));
	}
	else {
		ged_method_->run(g_id, h_id);
	}
	std::pair<GEDGraph::GraphID, GEDGraph::GraphID> key(g_id, h_id);
	lower_bounds_[key] = ged_method_->get_lower_bound();
	upper_bounds_[key] = ged_method_->get_upper_bound();
	runtimes_[key] = ged_method_->get_runtime();
	auto it = node_maps_.find(key);
	if (it == node_maps_.end()) {
		node_maps_.emplace(key, ged_method_->get_node_map());
	}
	else {
		it->second = ged_method_->get_node_map();
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::pair<GEDGraph::GraphID, GEDGraph::GraphID>
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
graph_ids() const {
	return std::make_pair(0, static_cast<GEDGraph::GraphID>(ged_data_.num_graphs_without_shuffled_copies()));
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
num_graphs() const {
	return ged_data_.num_graphs_without_shuffled_copies();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
num_node_labels() const {
	return ged_data_.node_labels_.size();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
UserNodeLabel
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_node_label(LabelID label_id) const {
	if (label_id < 1 or label_id > num_node_labels()) {
		throw Error("The environment does not contain a node label with ID " + std::to_string(label_id) + ".");
	}
	return ged_data_.node_labels_.at(label_id - 1);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
num_edge_labels() const {
	return ged_data_.edge_labels_.size();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
UserEdgeLabel
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_edge_label(LabelID label_id) const {
	if (label_id < 1 or label_id > num_edge_labels()) {
		throw Error("The environment does not contain an edge label with ID " + std::to_string(label_id) + ".");
	}
	return ged_data_.edge_labels_.at(label_id - 1);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
init_method() {
	if (not initialized_) {
		throw Error("The environment is uninitialized. Call init() before calling init_method().");
	}
	if (not ged_method_) {
		throw Error("No method has been set. Call set_method() before calling init_method().");
	}
	ged_method_->init();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_lower_bound(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const {
	std::pair<GEDGraph::GraphID, GEDGraph::GraphID> key(g_id, h_id);
	if (lower_bounds_.find(key) == lower_bounds_.end()) {
		throw Error("Call run(" + std::to_string(g_id) + "," + std::to_string(h_id) + ") before calling get_lower_bound(" + std::to_string(g_id) + "," + std::to_string(h_id) + ").");
	}
	return lower_bounds_.at(key);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_upper_bound(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const {
	std::pair<GEDGraph::GraphID, GEDGraph::GraphID> key(g_id, h_id);
	if (upper_bounds_.find(key) == upper_bounds_.end()) {
		throw Error("Call run(" + std::to_string(g_id) + "," + std::to_string(h_id) + ") before calling get_upper_bound(" + std::to_string(g_id) + "," + std::to_string(h_id) + ").");
	}
	return upper_bounds_.at(key);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const NodeMap &
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_node_map(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const {
	std::pair<GEDGraph::GraphID, GEDGraph::GraphID> key(g_id, h_id);
	if (node_maps_.find(key) == node_maps_.end()) {
		throw Error("Call run(" + std::to_string(g_id) + "," + std::to_string(h_id) + ") before calling get_node_map(" + std::to_string(g_id) + "," + std::to_string(h_id) + ").");
	}
	return node_maps_.at(key);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_runtime(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) const {
	std::pair<GEDGraph::GraphID, GEDGraph::GraphID> key(g_id, h_id);
	if (runtimes_.find(key) == runtimes_.end()) {
		throw Error("Call run(" + std::to_string(g_id) + "," + std::to_string(h_id) + ") before calling get_runtime(" + std::to_string(g_id) + "," + std::to_string(h_id) + ").");
	}
	return runtimes_.at(key).count();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const std::string &
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_graph_class(GEDGraph::GraphID graph_id) const {
	return ged_data_.graph_classes_.at(graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_nodes(GEDGraph::GraphID graph_id) const {
	return ged_data_.graph(graph_id).num_nodes();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_avg_num_nodes() const {
	std::size_t sum_num_nodes{0};
	for (std::size_t graph_id{0}; graph_id < ged_data_.num_graphs_without_shuffled_copies(); graph_id++) {
		sum_num_nodes += ged_data_.graph(graph_id).num_nodes();
	}
	return static_cast<double>(sum_num_nodes) / static_cast<double>(ged_data_.num_graphs_without_shuffled_copies());
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
node_rel_cost(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const {
	return ged_data_.edit_costs_->node_rel_cost_fun(node_label_1, node_label_2);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
node_del_cost(const UserNodeLabel & node_label) const {
	return ged_data_.edit_costs_->node_del_cost_fun(node_label);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
node_ins_cost(const UserNodeLabel & node_label) const {
	return ged_data_.edit_costs_->node_ins_cost_fun(node_label);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
UserNodeLabel
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
median_node_label(const std::vector<UserNodeLabel> & node_labels) const {
	return ged_data_.edit_costs_->median_node_label(node_labels);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
edge_rel_cost(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const {
	return ged_data_.edit_costs_->edge_rel_cost_fun(edge_label_1, edge_label_2);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
edge_del_cost(const UserEdgeLabel & edge_label) const {
	return ged_data_.edit_costs_->edge_del_cost_fun(edge_label);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
edge_ins_cost(const UserEdgeLabel & edge_label) const {
	return ged_data_.edit_costs_->edge_ins_cost_fun(edge_label);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
UserEdgeLabel
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
median_edge_label(const std::vector<UserEdgeLabel> & edge_labels) const {
	return ged_data_.edit_costs_->median_edge_label(edge_labels);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const std::string &
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_graph_name(GEDGraph::GraphID graph_id) const {
	return ged_data_.graph_names_.at(graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_graph(GEDGraph::GraphID graph_id, bool adj_matrix, bool adj_lists, bool edge_list) const {
	ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> exchange_graph;
	const GEDGraph & graph{ged_data_.graphs_.at(graph_id)};
	exchange_graph.id = graph.id();
	exchange_graph.num_nodes = graph.num_nodes();
	exchange_graph.num_edges = graph.num_edges();
	if (adj_matrix) {
		exchange_graph.adj_matrix = std::vector<std::vector<std::size_t>>(exchange_graph.num_nodes, std::vector<std::size_t>(exchange_graph.num_nodes, 0));
	}
	for (GEDGraph::NodeID node_id{0}; node_id < exchange_graph.num_nodes; node_id++) {
		exchange_graph.original_node_ids.emplace_back(internal_to_original_node_ids_.at(graph_id).at(node_id));
		exchange_graph.node_labels.emplace_back(ged_data_.node_labels_.at(graph.get_node_label(node_id) - 1));
		if (adj_lists) {
			exchange_graph.adj_lists.emplace_back();
		}
	}
	for (auto eitr = graph.edges(); eitr.first != eitr.second; eitr.first++) {
		GEDGraph::EdgeID edge(*eitr.first);
		UserEdgeLabel edge_label(ged_data_.edge_labels_.at(graph.get_edge_label(edge) - 1));
		std::size_t i(graph.tail(edge));
		std::size_t j(graph.head(edge));
		std::pair<std::size_t, std::size_t> edge_as_pair(i,j);
		std::pair<std::size_t, std::size_t> inversed_edge_as_pair(j,i);
		if (adj_matrix) {
			exchange_graph.adj_matrix[i][j] = 1;
			exchange_graph.adj_matrix[j][i] = 1;
			exchange_graph.edge_labels[edge_as_pair] = edge_label;
			exchange_graph.edge_labels[inversed_edge_as_pair] = edge_label;
		}
		if (adj_lists) {
			exchange_graph.adj_lists.at(i).emplace_back(j, edge_label);
			exchange_graph.adj_lists.at(j).emplace_back(i, edge_label);
		}
		if (edge_list) {
			exchange_graph.edge_list.emplace_back(edge_as_pair, edge_label);
		}
	}
	return exchange_graph;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_init_time() const {
	return ged_method_->get_init_time().count();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_induced_cost(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id, NodeMap & node_map) const {
	ged_data_.compute_induced_cost(ged_data_.graphs_.at(g_id), ged_data_.graphs_.at(h_id), node_map);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
quasimetric_costs() const {
	return ged_data_.quasimetric_costs();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
init(Options::InitType init_type, bool print_to_stdout) {

	// Throw an exception if no edit costs have been selected.
	if (not ged_data_.edit_costs_) {
		throw Error("No edit costs have been selected. Call set_edit_costs() before calling init().");
	}

	// Return if the environment is initialized.
	if (initialized_) {
		return;
	}

	// Set initialization type.
	ged_data_.init_type_ = init_type;

	// Construct shuffled graph copies if necessary.
	if (ged_data_.shuffled_graph_copies_available()) {
		for (auto graph_id : new_graph_ids_) {
			construct_shuffled_graph_copy_(graph_id);
		}
	}

	// Re-initialize adjacency matrices (also previously initialized graphs must be re-initialized because of possible re-allocation).
	ProgressBar progress(ged_data_.graphs_.size());
	if (print_to_stdout) {
		std::cout << "\rInitializing graphs: " << progress << std::flush;
	}
	for (auto & graph : ged_data_.graphs_) {
		if (not graph.initialized()) {
			graph.setup_adjacency_matrix();
			ged_data_.max_num_nodes_ = std::max(ged_data_.max_num_nodes_, graph.num_nodes());
			ged_data_.max_num_edges_ = std::max(ged_data_.max_num_edges_, graph.num_edges());
		}
		if (print_to_stdout) {
			progress.increment();
			std::cout << "\rInitializing graphs: " << progress << std::flush;
		}
	}
	if (print_to_stdout) {
		std::cout << "\n";
	}

	// Initialize cost matrices if necessary.
	if (ged_data_.eager_init_()) {
		ged_data_.init_cost_matrices_(print_to_stdout);
	}

	// Mark environment as initialized.
	initialized_ = true;
	new_graph_ids_.clear();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
initialized() const {
	return initialized_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
Options::InitType
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_init_type() const {
	return ged_data_.init_type_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDGraph::GraphID
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_or_clear_shuffled_graph_copy_(GEDGraph::GraphID graph_id) {
	GEDGraph::GraphID copied_graph_id{ged_data_.id_shuffled_graph_copy(graph_id)};
	if (copied_graph_id < ged_data_.num_graphs()) {
		ged_data_.graphs_[copied_graph_id].clear();
		original_to_internal_node_ids_[copied_graph_id].clear();
		ged_data_.strings_to_internal_node_ids_[copied_graph_id].clear();
		internal_to_original_node_ids_[copied_graph_id].clear();
		ged_data_.internal_node_ids_to_strings_[copied_graph_id].clear();
	}
	else if (copied_graph_id == ged_data_.num_graphs()) {
		ged_data_.graphs_.emplace_back(copied_graph_id);
		ged_data_.graph_names_.push_back(get_graph_name(graph_id));
		ged_data_.graph_classes_.push_back(get_graph_class(graph_id));
		original_to_internal_node_ids_.push_back(std::map<UserNodeID, GEDGraph::NodeID>());
		ged_data_.strings_to_internal_node_ids_.push_back(std::map<std::string, GEDGraph::NodeID>());
		internal_to_original_node_ids_.push_back(std::map<GEDGraph::NodeID, UserNodeID>());
		ged_data_.internal_node_ids_to_strings_.push_back(std::map<GEDGraph::NodeID, std::string>());
	}
	else {
		throw Error("Unexpected copied graph ID " + std::to_string(copied_graph_id) + " for graph with ID " + std::to_string(graph_id) + ". Number of graphs in environment: " + std::to_string(ged_data_.num_graphs()) + ".");
	}
	return copied_graph_id;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
construct_shuffled_graph_copy_(GEDGraph::GraphID graph_id) {
	GEDGraph::GraphID copied_graph_id{add_or_clear_shuffled_graph_copy_(graph_id)};
	const GEDGraph & graph{ged_data_.graph(graph_id)};
	std::vector<std::pair<UserNodeID, UserNodeLabel>> nodes;
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		nodes.emplace_back(internal_to_original_node_ids_.at(graph_id).at(*node), ged_data_.id_to_node_label(graph.get_node_label(*node)));
	}
	std::random_device rng_g;
	std::mt19937 urng_g(rng_g());
	std::shuffle(nodes.begin(), nodes.end(), urng_g);
	for (const auto & node : nodes) {
		add_node(copied_graph_id, node.first, node.second);
	}
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		add_edge(copied_graph_id, internal_to_original_node_ids_.at(graph_id).at(graph.tail(*edge)), internal_to_original_node_ids_.at(graph_id).at(graph.head(*edge)), ged_data_.id_to_edge_label(graph.get_edge_label(*edge)));
	}
}

}

#endif /* SRC_ENV_GED_ENV_IPP_ */
