/*!
 * @file  branch_uniform.ipp
 * @brief BranchUniform class definition.
 */

#ifndef SRC_METHODS_BRANCH_UNIFORM_IPP_
#define SRC_METHODS_BRANCH_UNIFORM_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
~BranchUniform() {}

template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
BranchUniform(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
sort_method_{COUNTING},
wildcard_option_{false},
sorted_edge_labels_() {}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	sorted_edge_labels_[graph.id()] = SortedUserEdgeLabels_(graph, sort_method_);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

	const GEDGraph::SizeTNodeMap & g_ids_to_nodes = this->ids_to_nodes_.at(g.id());
	const GEDGraph::SizeTNodeMap & h_ids_to_nodes = this->ids_to_nodes_.at(h.id());
	const SortedUserEdgeLabels_ & sorted_edge_labels_g = sorted_edge_labels_.at(g.id());
	const SortedUserEdgeLabels_ & sorted_edge_labels_h = sorted_edge_labels_.at(h.id());

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				if (wildcard_option_) {
					master_problem(row_in_master, col_in_master) = compute_wildcard_substitution_cost_(g, h, g_ids_to_nodes.at(row_in_master), h_ids_to_nodes.at(col_in_master));
				}
				else {
					master_problem(row_in_master, col_in_master) = compute_substitution_cost_(g, h, g_ids_to_nodes.at(row_in_master), h_ids_to_nodes.at(col_in_master), sorted_edge_labels_g, sorted_edge_labels_h);
				}
			}
			else if (row_in_master < g.num_nodes()) {
				master_problem(row_in_master, h.num_nodes()) = compute_deletion_cost_(g, g_ids_to_nodes.at(row_in_master));
			}
			else if (col_in_master < h.num_nodes()) {
				if (wildcard_option_) {
					master_problem(g.num_nodes(), col_in_master) = compute_wildcard_insertion_cost_(h, h_ids_to_nodes.at(col_in_master));
				}
				else {
					master_problem(g.num_nodes(), col_in_master) = compute_insertion_cost_(h, h_ids_to_nodes.at(col_in_master));
				}
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	sort_method_ = COUNTING;
	wildcard_option_ = false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "[--sort-method <arg>] [--wildcards <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "sort-method") {
		if (arg == "STD") {
			sort_method_ = STD;
		}
		else if (arg == "COUNTING") {
			sort_method_ = COUNTING;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option upper-bound. Usage: options = \"[--sort-method STD|COUNTING] [...]\"");
		}
		return true;
	}
	else if (option == "wildcards") {
		if (arg == "NO") {
			wildcard_option_ = false;
		}
		else if (arg == "YES") {
			wildcard_option_ = true;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option wildcards. Usage: options = \"[--wildcards YES|NO] [...]\"");
		}
		return true;
	}
	return false;
}

// === Definition of private class SortedUserEdgeLabels_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedUserEdgeLabels_ ::
SortedUserEdgeLabels_(const GEDGraph & g, SortMethod_ sort_method):
sorted_edge_labels_(){
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		sorted_edge_labels_[*node] = std::vector<LabelID>();
		for (auto edge = g.incident_edges(*node).first; edge != g.incident_edges(*node).second; edge++) {
			sorted_edge_labels_[*node].push_back(g.get_edge_label(*edge));
		}
		switch (sort_method) {
		case STD:
			std::sort(sorted_edge_labels_[*node].begin(), sorted_edge_labels_[*node].end());
			break;
		default:
			util::counting_sort(sorted_edge_labels_[*node].begin(), sorted_edge_labels_[*node].end());
			break;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedUserEdgeLabels_ ::
SortedUserEdgeLabels_():
sorted_edge_labels_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedUserEdgeLabels_ ::
operator=(const SortedUserEdgeLabels_ & sorted_edge_labels) {
	sorted_edge_labels_ = sorted_edge_labels.sorted_edge_labels_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const std::vector<LabelID> &
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedUserEdgeLabels_ ::
get_incident_labels(GEDGraph::NodeID node) const {
	return sorted_edge_labels_.at(node);
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k,
		const SortedUserEdgeLabels_ & sorted_edge_labels_g, const SortedUserEdgeLabels_ & sorted_edge_labels_h) const {

	// Collect node substitution cost.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

	// Compute the minimal edge edit costs.

	// Compute the size of the multiset intersection.
	std::size_t intersection_size{0};
	auto label_g = sorted_edge_labels_g.get_incident_labels(i).begin();
	auto label_h = sorted_edge_labels_h.get_incident_labels(k).begin();
	while ((label_g != sorted_edge_labels_g.get_incident_labels(i).end()) and (label_h != sorted_edge_labels_h.get_incident_labels(k).end())) {
		if (*label_g == *label_h) {
			intersection_size++;
			label_g++;
			label_h++;
		}
		else if (*label_g < *label_h) {
			label_g++;
		}
		else {
			label_h++;
		}
	}

	// Add edge edit cost to substitution cost.
	cost += static_cast<double>(std::max(g.degree(i), h.degree(k)) - intersection_size) * 0.5;

	// Return the overall substitution cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_wildcard_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const {

	double cost{0.0};

	if (h.get_node_label(k) != dummy_label()) {
		cost += this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k));
	}

	// Initialize subproblem.
	DMatrix subproblem(g.degree(i) + 1, h.degree(k) + 1, 0.0);

	// Collect edge deletion costs.
	std::size_t j{0};
	for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++, j++) {
		subproblem(j, h.degree(k)) = this->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label()) * 0.5;
	}

	// Collect edge insertion costs.
	std::size_t l{0};
	for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++, l++) {
		if (h.get_edge_label(*kl) != dummy_label()) {
			//std::cout << "label of edge " << *kl << " = " << h.get_edge_label(*kl) << "\n";
			subproblem(g.degree(i), l) = this->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl)) * 0.5;
		}
	}

	// Collect edge relabelling costs.
	j = 0;
	for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++, j++) {
		l = 0;
		for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++, l++) {
			if (h.get_edge_label(*kl) != dummy_label()) {
				subproblem(j, l) = this->ged_data_.edge_cost(g.get_edge_label(*ij), h.get_edge_label(*kl)) * 0.5;
			}
		}
	}

	// Solve subproblem.
	LSAPESolver subproblem_solver(subproblem);
	subproblem_solver.set_model(this->lsape_model_);
	subproblem_solver.solve();

	// Add edge edit cost to substitution cost.
	cost += subproblem_solver.minimal_cost();

	// Return the overall substitution cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const {
	// Collect node deletion cost.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), dummy_label())};

	// Collect edge deletion cost.
	cost += static_cast<double>(g.degree(i)) * 0.5;

	// Return overall deletion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const {
	// Collect node insertion cost.
	double cost{this->ged_data_.node_cost(dummy_label(), h.get_node_label(k))};

	// Collect edge insertion cost.
	cost += static_cast<double>(h.degree(k)) * 0.5;

	// Return overall insertion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_wildcard_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const {
	// Collect node insertion cost.
	double cost{this->ged_data_.node_cost(dummy_label(), h.get_node_label(k))};

	// Collect edge insertion cost.
	std::size_t num_incident_wildard_edges{0};
	for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++) {
		if (h.get_edge_label(*kl) == dummy_label()) {
			num_incident_wildard_edges++;
		}
	}
	cost += static_cast<double>(h.degree(k) - num_incident_wildard_edges) * 0.5;

	// Return overall insertion cost.
	return cost;
}

}

#endif /* SRC_METHODS_BRANCH_UNIFORM_IPP_ */
