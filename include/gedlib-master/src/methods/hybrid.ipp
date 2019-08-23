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
 * @file  hybrid.ipp
 * @brief Hybrid class definition.
 */

#ifndef SRC_METHODS_HYBRID_IPP_
#define SRC_METHODS_HYBRID_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
Hybrid<UserNodeLabel, UserEdgeLabel>::
~Hybrid() {}

template<class UserNodeLabel, class UserEdgeLabel>
Hybrid<UserNodeLabel, UserEdgeLabel>::
Hybrid(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
lsape_model_{LSAPESolver::FLWC},
num_threads_{1},
time_limit_in_sec_{0.0} {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Hybrid<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & gg, const GEDGraph & hh, Result & result) {
	Timer timer(time_limit_in_sec_);
	GEDGraph g(gg);
	GEDGraph h(hh);

	// Run Partition.
	Partition<UserNodeLabel, UserEdgeLabel> partition(this->ged_data_);
	Result partition_result;
	partition.run_as_util(g, h, partition_result);

	// Initialize and run BranchUniform without wildcards.
	std::string options(std::string("--lsape-model ") + to_string_(lsape_model_) + " --threads " + std::to_string(num_threads_));
	BranchUniform<UserNodeLabel, UserEdgeLabel> branch_uniform(this->ged_data_);
	branch_uniform.set_options(options);
	Result bu_result;
	branch_uniform.run_as_util(gg, hh, bu_result);

	// Run BranchUniform DFS over mismatching substructures to find minimal wildcard Branch matching cost.
	SubstructItr_ current_substruct{partition.get_unmatched_substructs_().cbegin()};
	SubstructItr_ end_substructs{partition.get_unmatched_substructs_().cend()};

	if (current_substruct == end_substructs) {
		result.set_lower_bound(bu_result.lower_bound());
	}
	else {
		options += " --wildcards YES";
		double lower_bound{std::numeric_limits<double>::infinity()};
		if (branch_uniform_dfs_(timer, options, g, h, current_substruct, end_substructs, lower_bound)) {
			result.set_lower_bound(std::max(partition_result.lower_bound() + lower_bound, bu_result.lower_bound()));
		}
		else {
			result.set_lower_bound(std::max(partition_result.lower_bound(), bu_result.lower_bound()));
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Hybrid<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	lsape_model_ = LSAPESolver::FLWC;
	num_threads_ = 1;
	time_limit_in_sec_ = 0.;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Hybrid<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "[--lsape-model <arg>] [--time-limit <arg>] [--threads <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Hybrid<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "lsape-model") {
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
		else if (arg  == "ECBP") {
			lsape_model_ = LSAPESolver::ECBP;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option lsape-model. Usage: options = \"[--lsape-model ECBP|EBP|FLWC|FLCC|FBP|SFBP|FBP0] [...]\"");
		}
		return true;
	}
	else if (option == "time-limit") {
		try {
			time_limit_in_sec_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option time-limit. Usage: options = \"[--time-limit <convertible to double>] [...]");
		}
		return true;
	}
	else if (option == "threads") {
		try {
			num_threads_ = std::stoi(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		return true;
	}
	return false;
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
bool
Hybrid<UserNodeLabel, UserEdgeLabel>::
branch_uniform_dfs_(const Timer & timer, const std::string & options, GEDGraph & g, GEDGraph & h, SubstructItr_ current_substruct, SubstructItr_ end_substructs, double & lower_bound) {
	if (timer.expired()) {
		return false;
	}

	// Base case: run BranchUniform on on g and h with wildcard labels and update lower bound.
	if (current_substruct == end_substructs) {
		Result bu_result;
		BranchUniform<UserNodeLabel, UserEdgeLabel> branch_uniform(this->ged_data_);
		branch_uniform.set_options(options);
		branch_uniform.run_as_util(g, h, bu_result);
		lower_bound = std::min(lower_bound, bu_result.lower_bound());
		return true;
	}

	// Recursively add wildcard labels to h.
	bool success{true};
	switch (current_substruct->type) {
	case Partition<UserNodeLabel, UserEdgeLabel>::NODE:
	h.set_label(current_substruct->node_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	break;

	case Partition<UserNodeLabel, UserEdgeLabel>::EDGE:
	h.set_label(current_substruct->edge_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	break;

	case Partition<UserNodeLabel, UserEdgeLabel>::NODE_EDGE:
	h.set_label(current_substruct->node_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	h.set_label(current_substruct->node_1, current_substruct->node_label_1);

	h.set_label(current_substruct->edge_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	break;

	case Partition<UserNodeLabel, UserEdgeLabel>::NODE_EDGE_NODE:
	h.set_label(current_substruct->node_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	h.set_label(current_substruct->node_1, current_substruct->node_label_1);

	h.set_label(current_substruct->edge_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	h.set_label(current_substruct->edge_1, current_substruct->edge_label_1);

	h.set_label(current_substruct->node_2, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	break;

	case Partition<UserNodeLabel, UserEdgeLabel>::NODE_EDGE_EDGE:
	h.set_label(current_substruct->node_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	h.set_label(current_substruct->node_1, current_substruct->node_label_1);

	h.set_label(current_substruct->edge_1, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	h.set_label(current_substruct->edge_1, current_substruct->edge_label_1);

	h.set_label(current_substruct->edge_2, dummy_label());
	success = success and branch_uniform_dfs_(timer, options, g, h, current_substruct + 1, end_substructs, lower_bound);
	if (not success) {
		return false;
	}
	break;
	}
	return success;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Hybrid<UserNodeLabel, UserEdgeLabel>::
to_string_(LSAPESolver::Model lsape_model) {
	std::string method_string("");
	switch (lsape_model) {
	case LSAPESolver::ECBP:
		method_string = "ECBP";
		break;
	case LSAPESolver::EBP:
		method_string = "EBP";
		break;
	case LSAPESolver::FLWC:
		method_string = "FLWC";
		break;
	case LSAPESolver::FLCC:
		method_string = "FLCC";
		break;
	case LSAPESolver::FBP:
		method_string = "FBP";
		break;
	case LSAPESolver::FBP0:
		method_string = "FBP0";
		break;
	case LSAPESolver::SFBP:
		method_string = "SFBP";
		break;
	}
	return method_string;
}

}

#endif /* SRC_METHODS_HYBRID_IPP_ */
