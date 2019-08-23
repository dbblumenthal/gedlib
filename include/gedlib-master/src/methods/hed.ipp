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
 * @file hed.ipp
 * @brief ged::HED class definition.
 */

#ifndef SRC_METHODS_HED_IPP_
#define SRC_METHODS_HED_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
HED<UserNodeLabel, UserEdgeLabel>::
~HED() {}

template<class UserNodeLabel, class UserEdgeLabel>
HED<UserNodeLabel, UserEdgeLabel>::
HED(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
lsape_model_{LSAPESolver::Model::ECBP},
num_threads_{1},
use_hed_for_edge_set_distances_{true} {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
HED<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {
	DMatrix lsape_instance(g.num_nodes() + 1, h.num_nodes() + 1);
	populate_instance_(g, h, lsape_instance);
	double hed{0};
	hed += lsape_instance.matrix().block(0, 0, g.num_nodes(), h.num_nodes() + 1).rowwise().minCoeff().sum();
	hed += lsape_instance.matrix().block(0, 0, g.num_nodes() + 1, h.num_nodes()).colwise().minCoeff().sum();
	result.set_lower_bound(hed);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
HED<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "threads") {
		try {
			num_threads_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		return true;
	}
	else if (option == "lsape-model") {
		if (arg == "EBP") {
			lsape_model_ = LSAPESolver::EBP;
		}
		else if (arg  == "FLWC") {
			lsape_model_ = LSAPESolver::FLWC;
		}
		else if (arg  == "FLCC") {
			lsape_model_ = LSAPESolver::FLCC;
		}
		else if (arg  == "FBP") {
			lsape_model_ = LSAPESolver::FBP;
		}
		else if (arg == "SFBP") {
			lsape_model_ = LSAPESolver::SFBP;
		}
		else if (arg  == "FBP0") {
			lsape_model_ = LSAPESolver::FBP0;
		}
		else if (arg  != "ECBP") {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option lsape-model. Usage: options = \"[--lsape-model ECBP|EBP|FLWC|FLCC|FBP|SFBP|FBP0] [...]\"");
		}
		return true;
	}
	else if (option == "edge-set-distances") {
		if (arg == "OPTIMAL") {
			use_hed_for_edge_set_distances_ = false;
		}
		else if (arg  != "HED") {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option edge-set-distances. Usage: options = \"[--edge-set-distances OPTIMAL|HED] [...]\"");
		}
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
HED<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "[--lsape-model <arg>] [--threads <arg>] [--edge-set-distances <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
HED<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	lsape_model_ = LSAPESolver::ECBP;
	num_threads_ = 1;
	use_hed_for_edge_set_distances_ = true;
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
HED<UserNodeLabel, UserEdgeLabel>::
populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance) const {

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < lsape_instance.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < lsape_instance.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				lsape_instance(row_in_master, col_in_master) = compute_substitution_cost_(g, h, row_in_master, col_in_master);
			}
			else if (row_in_master < g.num_nodes()) {
				lsape_instance(row_in_master, h.num_nodes()) = compute_deletion_cost_(g, row_in_master);
			}
			else if (col_in_master < h.num_nodes()) {
				lsape_instance(g.num_nodes(), col_in_master) = compute_insertion_cost_(h, col_in_master);
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
HED<UserNodeLabel, UserEdgeLabel>::
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
		if (use_hed_for_edge_set_distances_) {
			subproblem(j, l) /= 2.0;
		}
	}

	// Solve subproblem and update overall substitution costs.
	if (use_hed_for_edge_set_distances_) {
		cost += subproblem.matrix().block(0, 0, g.degree(i), h.degree(k) + 1).rowwise().minCoeff().sum();
		cost += subproblem.matrix().block(0, 0, g.degree(i) + 1, h.degree(k)).colwise().minCoeff().sum();
	}
	else {
		LSAPESolver subproblem_solver(&subproblem);
		subproblem_solver.set_model(this->lsape_model_);
		subproblem_solver.solve();
		cost += subproblem_solver.minimal_cost();
	}

	// Return the substitution costs divided by 2.
	return cost / 2.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
HED<UserNodeLabel, UserEdgeLabel>::
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
HED<UserNodeLabel, UserEdgeLabel>::
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



#endif /* SRC_METHODS_HED_IPP_ */
