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
 * @file  branch.ipp
 * @brief Branch class definition.
 */

#ifndef SRC_METHODS_BRANCH_IPP_
#define SRC_METHODS_BRANCH_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
Branch<UserNodeLabel, UserEdgeLabel>::
~Branch() {}

template<class UserNodeLabel, class UserEdgeLabel>
Branch<UserNodeLabel, UserEdgeLabel>::
Branch(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data) {}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Branch<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				master_problem(row_in_master, col_in_master) = compute_substitution_cost_(g, h, row_in_master, col_in_master);
			}
			else if (row_in_master < g.num_nodes()) {
				master_problem(row_in_master, h.num_nodes()) = compute_deletion_cost_(g, row_in_master);
			}
			else if (col_in_master < h.num_nodes()) {
				master_problem(g.num_nodes(), col_in_master) = compute_insertion_cost_(h, col_in_master);
			}
		}
	}
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
double
Branch<UserNodeLabel, UserEdgeLabel>::
compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const {
	// Collect node substitution costs.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

	// Initialize subproblem.
	DMatrix subproblem(g.degree(i) + 1, h.degree(k) + 1);

	// Collect edge deletion costs.
	std::size_t j{0};
	for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++, j++) {
		subproblem(j, h.degree(k)) = this->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label()) / 2.0;
	}

	// Collect edge insertion costs.
	std::size_t l{0};
	for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++, l++) {
		subproblem(g.degree(i), l) = this->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl)) / 2.0;
	}
	j = 0;

	// Collect edge relabelling costs.
	for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++, j++) {
		l = 0;
		for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++, l++) {
			subproblem(j, l) = this->ged_data_.edge_cost(g.get_edge_label(*ij), h.get_edge_label(*kl)) / 2.0;
		}
	}

	// Solve subproblem.
	LSAPESolver subproblem_solver(&subproblem);
	subproblem_solver.set_model(this->lsape_model_);
	subproblem_solver.solve();

	// Update and return overall substitution cost.
	cost += subproblem_solver.minimal_cost();
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Branch<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const {
	// Collect node deletion cost.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), ged::dummy_label())};

	// Collect edge deletion costs.
	auto incident_edges_i = g.incident_edges(i);
	for (auto ij = incident_edges_i.first; ij != incident_edges_i.second; ij++) {
		cost += this->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label()) / 2.0;
	}

	// Return overall deletion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Branch<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const {
	// Collect node insertion cost.
	double cost{this->ged_data_.node_cost(ged::dummy_label(), h.get_node_label(k))};

	// Collect edge insertion costs.
	auto incident_edges_k = h.incident_edges(k);
	for (auto kl = incident_edges_k.first; kl != incident_edges_k.second; kl++) {
		cost += this->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl)) / 2.0;
	}

	// Return overall insertion cost.
	return cost;
}

}

#endif /* SRC_METHODS_BRANCH_IPP_ */
