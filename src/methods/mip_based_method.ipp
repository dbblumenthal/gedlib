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
 * @file mip_based_method.ipp
 * @brief ged::MIPBasedMethod class definition.
 */

#ifndef SRC_METHODS_MIP_BASED_METHOD_IPP_
#define SRC_METHODS_MIP_BASED_METHOD_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
~MIPBasedMethod() {}

template<class UserNodeLabel, class UserEdgeLabel>
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
MIPBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
relax_{false},
map_root_to_root_{false},
num_threads_{1},
time_limit_in_sec_{0},
tune_{false},
tune_time_limit_in_sec_{0} {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	mip_init_();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {
	try {
		// Initialize optimization environment.
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Threads, num_threads_);
		if (time_limit_in_sec_ > 0) {
			env.set(GRB_DoubleParam_TimeLimit, time_limit_in_sec_);
		}
		if (tune_ and (tune_time_limit_in_sec_ > 0)) {
			env.set(GRB_DoubleParam_TuneTimeLimit, tune_time_limit_in_sec_);
		}

		// Populate and verify the model.
		GRBModel model = GRBModel(env);
		mip_populate_model_(g, h, model);
		if (model.get(GRB_IntAttr_IsMIP) and relax_) {
			throw Error("Relaxed model contains integral variables.");
		}

		// Tune and solve the model.
		if (tune_) {
			model.tune();
		}
		model.optimize();

		// Translate the solved model into a result.
		int model_status{model.get(GRB_IntAttr_Status)};
		if (model_status == GRB_INF_OR_UNBD) {
			throw Error("The model is infeasible or unbounded.");
		}
		else if (model_status == GRB_UNBOUNDED) {
			throw Error("The model is unbounded.");
		}
		else if (model_status == GRB_INFEASIBLE) {
			throw Error("The model is infeasible.");
		}
		else if ((model_status == GRB_OPTIMAL) or (model_status == GRB_SUBOPTIMAL) or (model_status == GRB_TIME_LIMIT)) {
			double obj_val{model.get(GRB_DoubleAttr_ObjVal) + mip_objective_constant_(g, h)};
			if (model_status == GRB_OPTIMAL) {
				result.set_lower_bound(obj_val);
			}
			if (not relax_) {
				std::size_t node_map_id{result.add_node_map(g.num_nodes(), h.num_nodes())};
				result.node_map(node_map_id).set_induced_cost(obj_val);
				mip_model_to_node_map_(g, h, model, result.node_map(node_map_id));
				if (map_root_to_root_ and (result.node_map(node_map_id).image(0) != 0)) {
					throw Error("Root constrained model does not map root to root.");
				}
			}
		}
		else {
			throw Error("Unexpected model status " + std::to_string(model_status) + ".");
		}
	}
	catch (const GRBException & e) {
		throw Error("Optimization failed with GRBException. Error code: " + std::to_string(e.getErrorCode()) + ". Error message: " + e.getMessage() + ".");
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	bool is_valid_option{false};
	if (option == "relax") {
		if (arg == "TRUE") {
			relax_ = true;
		}
		else if (arg != "FALSE") {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option relax. Usage: options = \"[--relax TRUE|FALSE] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "map-root-to-root") {
		if (arg == "TRUE") {
			map_root_to_root_ = true;
		}
		else if (arg != "FALSE") {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option map-root-to-root. Usage: options = \"[--map-root-to-root TRUE|FALSE] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "tune") {
		if (arg == "TRUE") {
			tune_ = true;
		}
		else if (arg != "FALSE") {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option tune. Usage: options = \"[--tune TRUE|FALSE] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "threads") {
		try {
			num_threads_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>]\"");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>]\"");
		}
		is_valid_option = true;
	}
	else if (option == "time-limit") {
		try {
			time_limit_in_sec_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option time-limit.  Usage: options = \"[--time-limit <convertible to double>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "tune-time-limit") {
		try {
			tune_time_limit_in_sec_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option tune-time-limit.  Usage: options = \"[--tune-time-limit <convertible to double>] [...]");
		}
		is_valid_option = true;
	}
	is_valid_option = is_valid_option or mip_parse_option_(option, arg);
	return is_valid_option;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	if (mip_valid_options_string_() == "") {
		return "[--relax <arg>] [--map-root-to-root <arg>] [--tune <arg>] [--threads <arg>] [--time-limit <arg>] [--tune-time-limit <arg>]";
	}
	return mip_valid_options_string_() + " [--relax <arg>] [--map-root-to-root <arg>] [--tune <arg>] [--threads <arg>] [--time-limit <arg>] [--tune-time-limit <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	relax_ = false;
	map_root_to_root_ = false;
	tune_ = false;
	num_threads_ = 1;
	time_limit_in_sec_ = 0;
	tune_time_limit_in_sec_ = 0;
	mip_set_default_options_();
}

// Default definitions of virtual member functions.

template<class UserNodeLabel, class UserEdgeLabel>
void
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model) {}

template<class UserNodeLabel, class UserEdgeLabel>
double
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
mip_objective_constant_(const GEDGraph & g, const GEDGraph & h) {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map) {}

template<class UserNodeLabel, class UserEdgeLabel>
void
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
mip_init_() {}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
mip_parse_option_(const std::string & option, const std::string & arg) {
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
mip_valid_options_string_() const {
	return "";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>::
mip_set_default_options_() {}

}



#endif /* SRC_METHODS_MIP_BASED_METHOD_IPP_ */
