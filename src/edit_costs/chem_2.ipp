/*!
 * @file chem_2.ipp
 * @brief ged::CHEM2 class definition.
 */

#ifndef SRC_EDIT_COSTS_CHEM_2_IPP_
#define SRC_EDIT_COSTS_CHEM_2_IPP_

namespace ged {

template<>
CHEM2<GXLLabel, GXLLabel>::
~CHEM2() {}

template<>
CHEM2<GXLLabel, GXLLabel>::
CHEM2(double node_ins_del_cost, double edge_ins_del_cost, double alpha) :
node_ins_del_cost_{node_ins_del_cost},
edge_ins_del_cost_{edge_ins_del_cost},
alpha_{alpha} {}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	if (node_label_1.at("chem") != node_label_2.at("chem")) {
		return alpha_ * 2 * node_ins_del_cost_;
	}
	return 0.0;
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_;
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_;
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	if (edge_label_1.at("valence") != edge_label_2.at("valence")) {
		return (1 - alpha_) * edge_ins_del_cost_;
	}
	return 0.0;
}

}




#endif /* SRC_EDIT_COSTS_CHEM_2_IPP_ */
