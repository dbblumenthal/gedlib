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
 * @file f2.ipp
 * @brief ged::F2 class definition.
 */

#ifndef SRC_METHODS_F2_IPP_
#define SRC_METHODS_F2_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
F2<UserNodeLabel, UserEdgeLabel>::
~F2() {}

template<class UserNodeLabel, class UserEdgeLabel>
F2<UserNodeLabel, UserEdgeLabel>::
F2(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
x_(),
y_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
F2<UserNodeLabel, UserEdgeLabel>::
mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model) {

	// Clear the variables and the constant for the objective.
	x_.clear();
	y_.clear();
	GRBLinExpr obj = 0.0;

	// Collect node deletion costs.
	std::vector<double> del_cost;
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		del_cost.push_back(this->ged_data_.node_cost(g.get_node_label(*i), dummy_label()));
		obj += del_cost.back();
	}

	// Collect node insertion costs.
	std::vector<double> ins_cost;
	for (auto k = h.nodes().first; k != h.nodes().second; k++) {
		ins_cost.push_back(this->ged_data_.node_cost(dummy_label(), h.get_node_label(*k)));
		obj += ins_cost.back();
	}

	// Add x variables to the model.
	std::pair<GEDGraph::NodeID, GEDGraph::NodeID> key_x;
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key_x.first = *i;
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key_x.second = *k;
			if ((*i == 0) and (*k == 0) and this->map_root_to_root_) {
				x_[key_x] = model.addVar(1, 1, 0, variable_type_());
			}
			else {
				x_[key_x] = model.addVar(0, 1, 0, variable_type_());
			}
			obj += (this->ged_data_.node_cost(g.get_node_label(*i), h.get_node_label(*k)) - del_cost.at(*i) - ins_cost.at(*k)) * x_.at(key_x);
		}
	}

	// Collect edge deletion costs.
	del_cost.clear();
	for (auto e = g.edges().first; e != g.edges().second; e++) {
		del_cost.push_back(this->ged_data_.edge_cost(g.get_edge_label(*e), dummy_label()));
		obj += del_cost.back();
	}


	// Collect edge insertion costs.
	ins_cost.clear();
	for (auto f = h.edges().first; f != h.edges().second; f++) {
		ins_cost.push_back(this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(*f)));
		obj += ins_cost.back();
	}

	// Add y variables to the model.
	std::pair<GEDGraph::EdgeID, GEDGraph::EdgeID> key_y;
	std::size_t counter_e{0};
	std::size_t counter_f{0};
	for (auto e = g.edges().first; e != g.edges().second; e++, counter_e++) {
		key_y.first = *e;
		counter_f = 0;
		for (auto f = h.edges().first; f != h.edges().second; f++, counter_f++) {
			key_y.second = *f;
			y_[key_y] = model.addVar(0, 1, 0, variable_type_());
			obj += (this->ged_data_.edge_cost(g.get_edge_label(*e), h.get_edge_label(*f)) - del_cost.at(counter_e) - ins_cost.at(counter_f)) * y_.at(key_y);
		}
	}

	// Set the objective.
	model.setObjective(obj, GRB_MINIMIZE);

	// Add node constraints to the model.
	GRBLinExpr lhs;
	for (auto i = g.nodes().first; i != g.nodes().second; i++) {
		key_x.first = *i;
		lhs = 0;
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key_x.second = *k;
			lhs += x_.at(key_x);
		}
		model.addConstr(lhs, GRB_LESS_EQUAL, 1);
	}
	for (auto k = h.nodes().first; k != h.nodes().second; k++) {
		key_x.second = *k;
		lhs = 0;
		for (auto i = g.nodes().first; i != g.nodes().second; i++) {
			key_x.first = *i;
			lhs += x_.at(key_x);
		}
		model.addConstr(lhs, GRB_LESS_EQUAL, 1);
	}

	// Add node-edge constraints to the model.
	for (auto e = g.edges().first; e != g.edges().second; e++, counter_e++) {
		key_y.first = *e;
		for (auto k = h.nodes().first; k != h.nodes().second; k++) {
			key_x.second = *k;
			lhs = 0;
			for (auto f = h.incident_edges(*k).first; f != h.incident_edges(*k).second; f++) {
				key_y.second = *f;
				lhs += y_.at(key_y);
			}
			key_x.first = g.tail(*e);
			lhs -= x_.at(key_x);
			key_x.first = g.head(*e);
			lhs -= x_.at(key_x);
			model.addConstr(lhs, GRB_LESS_EQUAL, 0);
		}
	}

}

template<class UserNodeLabel, class UserEdgeLabel>
void
F2<UserNodeLabel, UserEdgeLabel>::
mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map) {

	// Initialize local variables.
	std::vector<bool> delete_node(g.num_nodes(), true);
	std::vector<bool> insert_node(h.num_nodes(), true);
	GEDGraph::NodeID i, k;

	// Add node substitutions.
	for (auto it = x_.begin(); it != x_.end(); it++) {
		if (it->second.get(GRB_DoubleAttr_X) > 0) {
			i = it->first.first;
			k = it->first.second;
			node_map.add_assignment(i, k);
			delete_node[i] = false;
			insert_node[k] = false;
		}
	}

	// Add node deletions.
	for (i = 0; i < g.num_nodes(); i++) {
		if (delete_node.at(i)) {
			node_map.add_assignment(i, GEDGraph::dummy_node());
		}
	}

	// Add node insertions.
	for (k = 0; k < h.num_nodes(); k++) {
		if (insert_node.at(k)) {
			node_map.add_assignment(GEDGraph::dummy_node(), k);
		}
	}

	// Set induced cost.
	node_map.set_induced_cost(model.get(GRB_DoubleAttr_ObjVal));

}

template<class UserNodeLabel, class UserEdgeLabel>
bool
F2<UserNodeLabel, UserEdgeLabel>::
mip_model_to_lsape_projection_problem_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, DMatrix & lsape_instance) {

	std::vector<double> delete_node(g.num_nodes(), 0);
	std::vector<double> insert_node(h.num_nodes(), 0);

	GEDGraph::NodeID i, k;
	for (auto it = x_.begin(); it != x_.end(); it++) {
		i = it->first.first;
		k = it->first.second;
		lsape_instance(i, k) -= it->second.get(GRB_DoubleAttr_X);
		delete_node[i] += it->second.get(GRB_DoubleAttr_X);
		insert_node[k] += it->second.get(GRB_DoubleAttr_X);
	}
	for (i = 0; i < g.num_nodes(); i++) {
		lsape_instance(i, h.num_nodes()) = delete_node.at(i);
	}
	for (k = 0; k < h.num_nodes(); k++) {
		lsape_instance(g.num_nodes(), k) = insert_node.at(k);
	}
	return true;
}

template<class UserNodeLabel, class UserEdgeLabel>
char
F2<UserNodeLabel, UserEdgeLabel>::
variable_type_() const {
	if (this->relax_) {
		return GRB_CONTINUOUS;
	}
	return GRB_BINARY;
}

}

#endif /* SRC_METHODS_F2_IPP_ */
