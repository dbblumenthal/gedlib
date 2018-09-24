/*!
 * @file edit_costs.ipp
 * @brief ged::EditCosts class definition.
 */

#ifndef SRC_EDIT_COSTS_EDIT_COSTS_IPP_
#define SRC_EDIT_COSTS_EDIT_COSTS_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
EditCosts<UserNodeLabel, UserEdgeLabel>::
~EditCosts() {}

template<class UserNodeLabel, class UserEdgeLabel>
EditCosts<UserNodeLabel, UserEdgeLabel>::
EditCosts() {}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
node_del_cost_fun(const UserNodeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
node_ins_cost_fun(const UserNodeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
EditCosts<UserNodeLabel, UserEdgeLabel>::
vectorize_node_label(const UserNodeLabel & node_label, std::vector<double> & vector_representation) const {
	vector_representation.clear();
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
edge_rel_cost_fun(const UserEdgeLabel & node_label_1, const UserEdgeLabel & node_label_2) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
edge_del_cost_fun(const UserEdgeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
edge_ins_cost_fun(const UserEdgeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
EditCosts<UserNodeLabel, UserEdgeLabel>::
vectorize_edge_label(const UserEdgeLabel & edge_label, std::vector<double> & vector_representation) const {
	vector_representation.clear();
}

}

#endif /* SRC_EDIT_COSTS_EDIT_COSTS_IPP_ */
