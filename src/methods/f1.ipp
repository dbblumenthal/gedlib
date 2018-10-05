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
 * @file f1.ipp
 * @brief ged::F1 class definition.
 */

#ifndef SRC_METHODS_F1_IPP_
#define SRC_METHODS_F1_IPP_


namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
F1<UserNodeLabel, UserEdgeLabel>::
~F1() {}

template<class UserNodeLabel, class UserEdgeLabel>
F1<UserNodeLabel, UserEdgeLabel>::
F1(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
x_(),
y_(),
u_(),
v_(),
e_(),
f_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
F1<UserNodeLabel, UserEdgeLabel>::
mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model) {

	// Clear the variables.
	x_.clear();
	y_.clear();
	u_.clear();
	v_.clear();
	e_.clear();
	f_.clear();

	// Add u variables to the model (node deletions).
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		u_.push_back(model.addVar(0, 1, this->ged_data_.node_cost(g.get_node_label(*i), dummy_label()), variable_type_()));
	}

	// Add v variables to the model (node insertions).
	for (auto k = h.nodes().first; k != h.nodes().second; k++) {
		v_.push_back(model.addVar(0, 1, this->ged_data_.node_cost(dummy_label(), h.get_node_label(*k)), variable_type_()));
	}

	// Add x variables to the model (node substitutions).
	std::pair<GEDGraph::NodeID, GEDGraph::NodeID> key_x;
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key_x.first = *i;
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key_x.second = *k;
			if ((*i == 0) and (*k == 0) and this->map_root_to_root_) {
				x_[key_x] = model.addVar(1, 1, this->ged_data_.node_cost(g.get_node_label(*i), h.get_node_label(*k)), variable_type_());
			}
			else {
				x_[key_x] = model.addVar(0, 1, this->ged_data_.node_cost(g.get_node_label(*i), h.get_node_label(*k)), variable_type_());
			}
		}
	}

	// Add e variables to the model (edge deletions).
	for (auto e = g.edges().first; e != g.edges().second; e++) {
		e_[*e] = model.addVar(0, 1, this->ged_data_.edge_cost(g.get_edge_label(*e), dummy_label()), variable_type_());
	}

	// Add f variables to the model (edge insertions).
	for (auto f = h.edges().first; f != h.edges().second; f++) {
		f_[*f] = model.addVar(0, 1, this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(*f)), variable_type_());
	}

	// Add y variables to the model (edge substitutions).
	std::pair<GEDGraph::EdgeID, GEDGraph::EdgeID> key_y;
	for (auto e = g.edges().first; e != g.edges().second; e++) {
		key_y.first = *e;
		for (auto f = h.edges().first; f != h.edges().second; f++) {
			key_y.second = *f;
			y_[key_y] = model.addVar(0, 1, this->ged_data_.edge_cost(g.get_edge_label(*e), h.get_edge_label(*f)), variable_type_());
		}
	}

	// Add node constraints to the model.
	GRBLinExpr lhs;
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key_x.first = *i;
		lhs = u_.at(*i);
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key_x.second = *k;
			lhs += x_.at(key_x);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}
	for (auto k = h.nodes().first; k != h.nodes().second; k++) {
		key_x.second = *k;
		lhs = v_.at(*k);
		for (auto i = g.nodes().first; i != g.nodes().second; i++) {
			key_x.first = *i;
			lhs += x_.at(key_x);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}

	// Add edge constraints to the model.
	for (auto e = g.edges().first; e != g.edges().second; e++) {
		key_y.first = *e;
		lhs = e_.at(*e);
		for (auto f = h.edges().first; f != h.edges().second; f++) {
			key_y.second = *f;
			lhs += y_.at(key_y);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}
	for (auto f = h.edges().first; f != h.edges().second; f++) {
		key_y.second = *f;
		lhs = f_.at(*f);
		for (auto e = g.edges().first; e != g.edges().second; e++) {
			key_y.first = *e;
			lhs += y_.at(key_y);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}

	// Add topology constraints to the model.
	for (auto e = g.edges().first; e != g.edges().second; e++) {
		key_y.first = *e;
		GEDGraph::NodeID i{g.tail(*e)};
		GEDGraph::NodeID j{g.head(*e)};
		for (auto f = h.edges().first; f != h.edges().second; f++) {
			key_y.second = *f;
			GEDGraph::NodeID k{h.tail(*f)};
			GEDGraph::NodeID l{h.head(*f)};
			key_x.first = i;
			key_x.second = k;
			lhs = x_.at(key_x);
			key_x.second = l;
			lhs += x_.at(key_x);
			model.addConstr(lhs, GRB_GREATER_EQUAL, y_.at(key_y));
			key_x.first = j;
			lhs = x_.at(key_x);
			key_x.second = k;
			lhs += x_.at(key_x);
			model.addConstr(lhs, GRB_GREATER_EQUAL, y_.at(key_y));
		}
	}

}

template<class UserNodeLabel, class UserEdgeLabel>
void
F1<UserNodeLabel, UserEdgeLabel>::
mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map) {

	// Add node substitutions.
	GEDGraph::NodeID i, k;
	for (auto it = x_.begin(); it != x_.end(); it++) {
		if (it->second.get(GRB_DoubleAttr_X) > 0) {
			i = it->first.first;
			k = it->first.second;
			node_map.add_assignment(i, k);
		}
	}

	// Add node deletions.
	for (i = 0; i < g.num_nodes(); i++) {
		if (u_.at(i).get(GRB_DoubleAttr_X) > 0) {
			node_map.add_assignment(i, GEDGraph::dummy_node());
		}
	}

	// Add node insertions.
	for (k = 0; k < h.num_nodes(); k++) {
		if (v_.at(k).get(GRB_DoubleAttr_X) > 0) {
			node_map.add_assignment(GEDGraph::dummy_node(), k);
		}
	}

	// Set induced cost.
	node_map.set_induced_cost(model.get(GRB_DoubleAttr_ObjVal));

}

template<class UserNodeLabel, class UserEdgeLabel>
bool
F1<UserNodeLabel, UserEdgeLabel>::
mip_model_to_lsape_projection_problem_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, DMatrix & lsape_instance) {
	GEDGraph::NodeID i, k;
	for (auto it = x_.begin(); it != x_.end(); it++) {
		i = it->first.first;
		k = it->first.second;
		lsape_instance(i, k) -= it->second.get(GRB_DoubleAttr_X);
	}
	for (i = 0; i < g.num_nodes(); i++) {
		lsape_instance(i, h.num_nodes()) -= u_.at(i).get(GRB_DoubleAttr_X);
	}
	for (k = 0; k < h.num_nodes(); k++) {
		lsape_instance(g.num_nodes(), k) -= v_.at(k).get(GRB_DoubleAttr_X);
	}
	return true;
}

template<class UserNodeLabel, class UserEdgeLabel>
char
F1<UserNodeLabel, UserEdgeLabel>::
variable_type_() const {
	if (this->relax_) {
		return GRB_CONTINUOUS;
	}
	return GRB_BINARY;
}

}


#endif /* SRC_METHODS_F1_IPP_ */
