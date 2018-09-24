/*!
 * @file  node.ipp
 * @brief Node class definition.
 */

#ifndef SRC_METHODS_NODE_IPP_
#define SRC_METHODS_NODE_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
Node<UserNodeLabel, UserEdgeLabel>::
~Node() {}

template<class UserNodeLabel, class UserEdgeLabel>
Node<UserNodeLabel, UserEdgeLabel>::
Node(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data) {}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Node<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

	const GEDGraph::SizeTNodeMap & g_ids_to_nodes = this->ids_to_nodes_.at(g.id());
	const GEDGraph::SizeTNodeMap & h_ids_to_nodes = this->ids_to_nodes_.at(h.id());

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				master_problem(row_in_master, col_in_master) = this->ged_data_.node_cost(g.get_node_label(g_ids_to_nodes.at(row_in_master)), h.get_node_label(h_ids_to_nodes.at(col_in_master)));
			}
			else if (row_in_master < g.num_nodes()) {
				master_problem(row_in_master, h.num_nodes()) = this->ged_data_.node_cost(g.get_node_label(g_ids_to_nodes.at(row_in_master)), dummy_label());
			}
			else if (col_in_master < h.num_nodes()) {
				master_problem(g.num_nodes(), col_in_master) = this->ged_data_.node_cost(dummy_label(), h.get_node_label(h_ids_to_nodes.at(col_in_master)));
			}
		}
	}
}

}

#endif /* SRC_METHODS_NODE_IPP_ */
