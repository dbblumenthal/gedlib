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
set_edit_costs(Options::EditCosts edit_costs, std::initializer_list<double> edit_cost_constants) {
	ged_data_.set_edit_costs_(edit_costs, std::vector<double>(edit_cost_constants));
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_edit_costs(EditCosts<UserNodeLabel, UserEdgeLabel> * edit_costs) {
	ged_data_.set_edit_costs_(edit_costs);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDGraph::GraphID
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_graph(const std::string & name, const std::string & graph_class) {
	if (initialized_) {
		throw Error("The environment is has already been initialized. Don't call add_graph() after calling init().");
	}
	GEDGraph::GraphID graph_id{ged_data_.num_graphs()};
	ged_data_.graphs_.emplace_back(graph_id);
	ged_data_.graph_names_.push_back(name);
	ged_data_.graph_classes_.push_back(graph_class);
	original_to_internal_node_ids_.push_back(std::map<UserNodeID, GEDGraph::NodeID>());
	ged_data_.original_to_internal_node_ids_.push_back(std::map<std::string, GEDGraph::NodeID>());
	internal_to_original_node_ids_.push_back(std::map<GEDGraph::NodeID, UserNodeID>());
	ged_data_.internal_to_original_node_ids_.push_back(std::map<GEDGraph::NodeID, std::string>());
	return graph_id;
}

template<>
GEDGraph::GraphID
GEDEnv<GXLNodeID, GXLLabel, GXLLabel>::
read_graph_from_gxl_(const std::string & dir, const std::string & filename, const std::string & graph_class, Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
		const std::unordered_set<std::string> & irrelevant_node_attributes, const std::unordered_set<std::string> & irrelevant_edge_attributes) {

	// read the file into a property tree
	boost::property_tree::ptree root;
	try {
		read_xml(dir + filename, root);
	}
	catch (const boost::property_tree::xml_parser_error & error) {
		throw Error(std::string("Error reading file ") + filename + ": " + error.message() + ".");
	}

	// first sanity checks
	if (root.count("gxl") == 0) {
		throw Error("The file " + filename + " has the wrong format: no xml-element <gxl>.");
	}
	if (root.count("gxl") >= 2) {
		throw Error("The file " + filename + " has the wrong format: more than one xml-element <gxl>.");
	}
	root = root.get_child("gxl");
	if (root.count("graph") == 0) {
		throw Error("The file " + filename + " has the wrong format: no xml-element <gxl>.<graph>");
	}
	if (root.count("graph") >= 2) {
		throw Error("The file " + filename + " has the wrong format: more than one xml-element <gxl>.<graph>");
	}
	root = root.get_child("graph");

	// add new graph to the environment
	GEDGraph::GraphID graph_id{add_graph(filename, graph_class)};

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
			graph_ids.push_back(read_graph_from_gxl_(graph_dir, gxl_file, graph_class, node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes));
		}
		else if (val.first != "<xmlattr>") {
			throw Error("The file " + file + " has the wrong format: unexpected element <GraphCollection>.<" + val.first + ">.");
		}
	}
	return graph_ids;
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
to_string_(UserNodeID node_id) {
	std::stringstream ss;
	ss << node_id;
	return ss.str();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_node(GEDGraph::GraphID graph_id, const UserNodeID & node_id, const UserNodeLabel & node_label) {
	if (graph_id >= ged_data_.num_graphs()) {
		throw Error("The graph " + get_graph_name(graph_id) + " has not been added to the environment.");
	}
	if (initialized_) {
		throw Error("The environment is has already been initialized. Don't call add_node() after calling init().");
	}
	if (original_to_internal_node_ids_[graph_id].find(node_id) != original_to_internal_node_ids_[graph_id].end()) {
		throw Error("The node " + to_string_(node_id) + " has already been added to the graph " + std::to_string(graph_id) + ": " + get_graph_name(graph_id) + ".");
	}
	LabelID label_id{ged_data_.node_label_to_id_(node_label)};
	GEDGraph::NodeID internal_node_id{ged_data_.graphs_[graph_id].add_node()};
	original_to_internal_node_ids_[graph_id][node_id] = internal_node_id;
	internal_to_original_node_ids_[graph_id][internal_node_id] = node_id;
	ged_data_.original_to_internal_node_ids_[graph_id][to_string_(node_id)] = internal_node_id;
	ged_data_.internal_to_original_node_ids_[graph_id][internal_node_id] = to_string_(node_id);
	ged_data_.graphs_[graph_id].set_label(original_to_internal_node_ids_[graph_id][node_id], label_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_edge(GEDGraph::GraphID graph_id, const UserNodeID & from, const UserNodeID & to, const UserEdgeLabel & edge_label, bool ignore_duplicates) {
	if (graph_id >= ged_data_.num_graphs()) {
		throw Error("The graph " + get_graph_name(graph_id) + " has not been added to the environment.");
	}
	if (initialized_) {
		throw Error("The environment is has already been initialized. Don't call add_edge() after calling init().");
	}
	if (original_to_internal_node_ids_[graph_id].find(from) == original_to_internal_node_ids_[graph_id].end()) {
		throw Error("The node " + to_string_(from) + " does not exist in the graph " + get_graph_name(graph_id) + ".");
	}
	if (original_to_internal_node_ids_[graph_id].find(to) == original_to_internal_node_ids_[graph_id].end()) {
		throw Error("The node " + to_string_(to) + " does not exist in the graph " + get_graph_name(graph_id) + ".");
	}
	if (ged_data_.graphs_[graph_id].safe_is_edge(original_to_internal_node_ids_[graph_id][from], original_to_internal_node_ids_[graph_id][to])) {
		if (ignore_duplicates) {
			return;
		}
		throw Error("An edge between " + to_string_(from) + " and " + to_string_(to) + " has already been added to the graph " + get_graph_name(graph_id) + ".");
	}
	LabelID label_id{ged_data_.edge_label_to_id_(edge_label)};
	GEDGraph::EdgeID edge_id{ged_data_.graphs_[graph_id].add_edge(original_to_internal_node_ids_[graph_id][from], original_to_internal_node_ids_[graph_id][to])};
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
run_method(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) {
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
	// call selected GEDApproximator and store results
	if (ged_data_.shuffled_graph_copies_available() and (g_id == h_id)) {
		ged_method_->run(g_id, ged_data_.id_shuffled_graph_copy(h_id));
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
const std::string &
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_graph_name(GEDGraph::GraphID graph_id) const {
	return ged_data_.graph_names_.at(graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_init_time() const {
	return ged_method_->get_init_time().count();
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
init(Options::InitType init_type) {
	if (initialized_) {
		throw Error("The environment is has already been initialized. Don't call init() twice.");
	}
	if (not ged_data_.edit_costs_) {
		throw Error("No edit costs have been selected. Call set_edit_costs() before calling init().");
	}
	ged_data_.init_type_ = init_type;
	if (ged_data_.shuffled_graph_copies_available()) {
		std::size_t num_graphs{ged_data_.num_graphs()};
		for (std::size_t graph_id{0}; graph_id < num_graphs; graph_id++) {
			construct_shuffled_graph_copy_(graph_id);
		}
	}
	ged_data_.init_(init_type);
	original_to_internal_node_ids_.clear();
	internal_to_original_node_ids_.clear();
	initialized_ = true;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel>::
construct_shuffled_graph_copy_(GEDGraph::GraphID graph_id) {
	GEDGraph::GraphID copied_graph_id{add_graph(get_graph_name(graph_id), get_graph_class(graph_id))};
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
