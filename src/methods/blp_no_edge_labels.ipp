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
 * @file blp_no_edge_labels.ipp
 * @brief ged::BLPNoEdgeLabels class definition.
 */

#ifndef SRC_METHODS_BLP_NO_EDGE_LABELS_IPP_
#define SRC_METHODS_BLP_NO_EDGE_LABELS_IPP_


namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
BLPNoEdgeLabels<UserNodeLabel, UserEdgeLabel>::
~BLPNoEdgeLabels() {}

template<class UserNodeLabel, class UserEdgeLabel>
BLPNoEdgeLabels<UserNodeLabel, UserEdgeLabel>::
BLPNoEdgeLabels(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
MIPBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
x_(),
s_(),
t_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
BLPNoEdgeLabels<UserNodeLabel, UserEdgeLabel>::
mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model) {

	// Clear the variables and the constant for the objective.
	x_.clear();
	s_.clear();
	t_.clear();

	// Determine minimal edge deletion or insertion cost.
	double min_edge_ins_del_cost{std::min(this->ged_data_.min_edge_del_cost(g), this->ged_data_.min_edge_ins_cost(h))};

	// Determine N.
	std::size_t N{g.num_nodes() + h.num_nodes()};

	// Add substitution variables to the model.
	NodeMap::Assignment key;
	for (GEDGraph::GraphID i{0}; i < N; i++) {
		for (GEDGraph::GraphID k{0}; k < N; k++) {
			key.first = i;
			key.second = k;
			if ((i == 0) and (k == 0) and this->map_root_to_root_) {
				x_[key] = model.addVar(1, 1, this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k)), variable_type_());
			}
			else if ((i < g.num_nodes()) and (k < h.num_nodes())) {
				x_[key] = model.addVar(0, 1, this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k)), variable_type_());
			}
			else if (i < g.num_nodes()) {
				x_[key] = model.addVar(0, 1, this->ged_data_.node_cost(g.get_node_label(i), dummy_label()), variable_type_());
			}
			else if (k < h.num_nodes()) {
				x_[key] = model.addVar(0, 1, this->ged_data_.node_cost(dummy_label(), h.get_node_label(k)), variable_type_());
			}
			else {
				x_[key] = model.addVar(0, 1, 0, variable_type_());
			}
			s_[key] = model.addVar(0, 1, min_edge_ins_del_cost / 2, variable_type_());
			t_[key] = model.addVar(0, 1, min_edge_ins_del_cost / 2, variable_type_());
		}
	}

	// Add node constraints to the model.
	GRBLinExpr lhs;
	for (GEDGraph::GraphID i{0}; i < N; i++) {
		key.first = i;
		lhs = 0;
		for (GEDGraph::GraphID k{0}; k < N; k++) {
			key.second = k;
			lhs += x_.at(key);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}
	for (GEDGraph::GraphID k{0}; k < N; k++) {
		key.second = k;
		lhs = 0;
		for (GEDGraph::GraphID i{0}; i < N; i++) {
			key.first = i;
			lhs += x_.at(key);
		}
		model.addConstr(lhs, GRB_EQUAL, 1);
	}

	// Add topology constraints to the model.
	for (GEDGraph::GraphID i{0}; i < N; i++) {
		for (GEDGraph::GraphID k{0}; k < N; k++) {
			key.first = i;
			key.second = k;
			lhs = s_.at(key) - t_.at(key);
			if (i < g.num_nodes()) {
				for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++) {
					key.first = g.head(*ij);
					lhs += x_.at(key);
				}
			}
			if (k < h.num_nodes()) {
				key.first = i;
				for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++) {
					key.second = h.head(*kl);
					lhs -= x_.at(key);
				}
			}
			model.addConstr(lhs, GRB_EQUAL, 0);
		}
	}


}

template<class UserNodeLabel, class UserEdgeLabel>
void
BLPNoEdgeLabels<UserNodeLabel, UserEdgeLabel>::
mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map) {

	// Initialize local variables.
	std::vector<bool> delete_node(g.num_nodes(), true);
	std::vector<bool> insert_node(h.num_nodes(), true);
	GEDGraph::NodeID i, k;

	// Add node substitutions.
	for (auto it = x_.begin(); it != x_.end(); it++) {
		i = it->first.first;
		k = it->first.second;
		if ((i < g.num_nodes()) and (k < h.num_nodes()) and (it->second.get(GRB_DoubleAttr_X) > 0)) {
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
	this->ged_data_.compute_induced_cost(g, h, node_map);

}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BLPNoEdgeLabels<UserNodeLabel, UserEdgeLabel>::
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
BLPNoEdgeLabels<UserNodeLabel, UserEdgeLabel>::
variable_type_() const {
	if (this->relax_) {
		return GRB_CONTINUOUS;
	}
	return GRB_BINARY;
}

}


#endif /* SRC_METHODS_BLP_NO_EDGE_LABELS_IPP_ */
