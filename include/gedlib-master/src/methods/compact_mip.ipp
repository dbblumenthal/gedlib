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
 * @file compact_mip.ipp
 * @brief ged::CompactMIP class definition.
 */

#ifndef SRC_METHODS_COMPACT_MIP_IPP_
#define SRC_METHODS_COMPACT_MIP_IPP_


namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
CompactMIP<UserNodeLabel, UserEdgeLabel>::
~CompactMIP() {}

template<class UserNodeLabel, class UserEdgeLabel>
CompactMIP<UserNodeLabel, UserEdgeLabel>::
CompactMIP(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
x_(),
z_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
CompactMIP<UserNodeLabel, UserEdgeLabel>::
mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model) {

	// Clear the variables and the constant for the objective.
	x_.clear();
	z_.clear();

	// Add substitution variables to the model.
	std::pair<GEDGraph::NodeID, GEDGraph::NodeID> key;
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key.first = *i;
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key.second = *k;
			if ((*i == 0) and (*k == 0) and this->map_root_to_root_) {
				x_[key] = model.addVar(1, 1, 0, variable_type_());
			}
			else {
				x_[key] = model.addVar(0, 1, 0, variable_type_());
			}
			z_[key] = model.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS);
		}
	}

	// Add deletion variables to the model.
	key.second = GEDGraph::dummy_node();
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key.first = *i;
		x_[key] = model.addVar(0, 1, 0, variable_type_());
		z_[key] = model.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS);
	}

	// Add insertion variables to the model.
	key.first = GEDGraph::dummy_node();
	for (auto k = h.nodes().first; k != h.nodes().second; k++) {
		key.second = *k;
		x_[key] = model.addVar(0, 1, 0, variable_type_());
		z_[key] = model.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS);
	}


	// Add node constraints to the model.
	GRBLinExpr lhs;
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key.first = *i;
		key.second = GEDGraph::dummy_node();
		lhs = x_.at(key);
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key.second = *k;
			lhs += x_.at(key);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}
	for (auto k = h.nodes().first; k != h.nodes().second; k++) {
		key.second = *k;
		key.first = GEDGraph::dummy_node();
		lhs = x_.at(key);
		for (auto i = g.nodes().first; i != g.nodes().second; i++) {
			key.first = *i;
			lhs += x_.at(key);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}

	// Add substitution constraints to the model.
	double u{0};
	double edit_cost{0};
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key.first = *i;
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key.second = *k;
			edit_cost = this->ged_data_.node_cost(g.get_node_label(*i), h.get_node_label(*k));
			lhs = edit_cost - z_.at(key);
			u = edit_cost;
			for (auto j = g.nodes().first; j != g.nodes().second; j++) {
				key.first = *j;
				for (auto l = h.nodes().first; l != h.nodes().second; l++) {
					key.second = *l;
					edit_cost = this->ged_data_.edge_cost(g.get_edge_label(*i, *j),  h.get_edge_label(*k, *l)) / 2;
					u += edit_cost;
					lhs += edit_cost * x_.at(key);
				}
			}
			key.second = GEDGraph::dummy_node();
			for (auto j = g.nodes().first; j != g.nodes().second; j++) {
				key.first = *j;
				edit_cost = this->ged_data_.edge_cost(g.get_edge_label(*i, *j), dummy_label()) / 2;
				u += edit_cost;
				lhs += edit_cost * x_.at(key);
			}
			key.first = GEDGraph::dummy_node();
			for (auto l = h.nodes().first; l != h.nodes().second; l++) {
				key.second = *l;
				edit_cost = this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(*k, *l)) / 2;
				u += edit_cost;
				lhs += edit_cost * x_.at(key);
			}
			key.first = *i;
			key.second = *k;
			lhs -= (1 - x_.at(key)) * u;
			model.addConstr(lhs, GRB_LESS_EQUAL, 0);
		}
	}

	// Add deletion constraints to the model.
	key.second = GEDGraph::dummy_node();
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key.first = *i;
		edit_cost = this->ged_data_.node_cost(g.get_node_label(*i), dummy_label());
		lhs = edit_cost - z_.at(key);
		u = edit_cost;
		for (auto j = g.nodes().first; j != g.nodes().second; j++) {
			key.first = *j;
			for (auto l = h.nodes().first; l != h.nodes().second; l++) {
				key.second = *l;
				edit_cost = this->ged_data_.edge_cost(g.get_edge_label(*i, *j), dummy_label()) / 2;
				u += edit_cost;
				lhs += edit_cost * x_.at(key);
			}
		}
		key.second = GEDGraph::dummy_node();
		for (auto j = g.nodes().first; j != g.nodes().second; j++) {
			key.first = *j;
			edit_cost = this->ged_data_.edge_cost(g.get_edge_label(*i, *j), dummy_label()) / 2;
			u += edit_cost;
			lhs += edit_cost * x_.at(key);
		}
		key.first = *i;
		lhs -= (1 - x_.at(key)) * u;
		model.addConstr(lhs, GRB_LESS_EQUAL, 0);
	}

	// Add insertion constraints to the model.
	key.first = GEDGraph::dummy_node();
	for (auto k = h.nodes().first; k != h.nodes().second; k++) {
		key.second = *k;
		edit_cost = this->ged_data_.node_cost(dummy_label(), h.get_node_label(*k));
		lhs = edit_cost - z_.at(key);
		u = edit_cost;
		for (auto j = g.nodes().first; j != g.nodes().second; j++) {
			key.first = *j;
			for (auto l = h.nodes().first; l != h.nodes().second; l++) {
				key.second = *l;
				edit_cost = this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(*k, *l)) / 2;
				u += edit_cost;
				lhs += edit_cost * x_.at(key);
			}
		}
		key.first = GEDGraph::dummy_node();
		for (auto l = h.nodes().first; l != h.nodes().second; l++) {
			key.second = *l;
			edit_cost = this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(*k, *l)) / 2;
			u += edit_cost;
			lhs += edit_cost * x_.at(key);
		}
		key.second = *k;
		lhs -= (1 - x_.at(key)) * u;
		model.addConstr(lhs, GRB_LESS_EQUAL, 0);
	}


}

template<class UserNodeLabel, class UserEdgeLabel>
void
CompactMIP<UserNodeLabel, UserEdgeLabel>::
mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map) {

	GEDGraph::NodeID i, k;
	for (auto it = x_.begin(); it != x_.end(); it++) {
		if (it->second.get(GRB_DoubleAttr_X) > 0) {
			i = it->first.first;
			k = it->first.second;
			node_map.add_assignment(i, k);
		}
	}

	// Set induced cost.
	node_map.set_induced_cost(model.get(GRB_DoubleAttr_ObjVal));

}

template<class UserNodeLabel, class UserEdgeLabel>
bool
CompactMIP<UserNodeLabel, UserEdgeLabel>::
mip_model_to_lsape_projection_problem_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, DMatrix & lsape_instance) {
	GEDGraph::NodeID i, k;
	for (auto it = x_.begin(); it != x_.end(); it++) {
		i = std::min(it->first.first, g.num_nodes());
		k = std::min(it->first.second, h.num_nodes());
		lsape_instance(i, k) -= it->second.get(GRB_DoubleAttr_X);
	}
	return true;
}

template<class UserNodeLabel, class UserEdgeLabel>
char
CompactMIP<UserNodeLabel, UserEdgeLabel>::
variable_type_() const {
	if (this->relax_) {
		return GRB_CONTINUOUS;
	}
	return GRB_BINARY;
}

}


#endif /* SRC_METHODS_COMPACT_MIP_IPP_ */
