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
 * @file  ring.ipp
 * @brief ged::Ring class definition.
 */

#ifndef SRC_METHODS_RING_IPP_
#define SRC_METHODS_RING_IPP_

namespace ged {

// === Definitions of destructors and constructors. ===
template<class UserNodeLabel, class UserEdgeLabel>
Ring<UserNodeLabel, UserEdgeLabel>::
~Ring() {}

template<class UserNodeLabel, class UserEdgeLabel>
Ring<UserNodeLabel, UserEdgeLabel>::
Ring(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
rings_(),
led_method_{LSAPE_OPTIMAL},
sort_method_{COUNTING},
num_layers_{undefined()},
alpha_(),
lambda_(),
mu_{1.0},
num_evals_{50},
num_x0s_{100},
infile_(""),
outfile_("") {
	this->compute_lower_bound_ = false;
}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	build_rings_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	led_method_ = LSAPE_OPTIMAL;
	sort_method_ = COUNTING;
	mu_ = 1;
	num_evals_ = 50;
	num_x0s_ = 100;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_pre_graph_init_(bool called_at_runtime) {
	if (load_config_file_()) {
		read_params_from_file_();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_default_post_graph_init_() {
	if (not load_config_file_())  {
		set_num_layers_();
		for (std::size_t level{0}; level < num_layers_; level++) {
			lambda_.push_back(1.0 / static_cast<double>(num_layers_));
		}
		for (std::size_t i{0}; i < 3; i++) {
			alpha_.push_back(1.0 / 3.0);
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_init_() {

	// Return if the method is initialized from a file.
	if (load_config_file_()) {
		return;
	}

	// Set the number of layers to the maximum diameter of the graphs in the instance.
	set_num_layers_();

	// Initialise the starting points employed by NOMAD.
	std::vector<NOMAD::Point> x0s;
	init_x0s_(x0s);

	// Run NOMAD.
	std::vector<NOMAD::Point> solutions;
	std::vector<NOMAD::Double> objectives;
	ProgressBar progress_bar(num_x0s_);
	std::cout << "\rNOMAD: " << progress_bar << std::flush;
	// Run NOMAD for all the starting points.

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t i = 0; i < num_x0s_; i++) {
		NOMAD::Display out (std::cout);
		out.precision (NOMAD::DISPLAY_PRECISION_STD);
		try {
			// Setup NOMAD.
			NOMAD::begin(0, nullptr);
			NOMAD::Parameters p(out);
			std::vector<NOMAD::bb_output_type> bbot(3);
			bbot[0] = NOMAD::OBJ; // objective function
			bbot[1] = NOMAD::PB; // sum over alpha > 0
			bbot[2] = NOMAD::PB; // sum over lambda > 0
			p.set_BB_OUTPUT_TYPE(bbot);
			p.set_DISPLAY_STATS ("bbe ( sol ) obj");
			p.set_DIMENSION(3 + num_layers_);
			p.set_LOWER_BOUND(NOMAD::Point(3 + num_layers_, 0.0));
			p.set_UPPER_BOUND(NOMAD::Point(3 + num_layers_, 1.0));
			p.set_MAX_BB_EVAL(num_evals_);
			p.set_DISPLAY_DEGREE("NO_DISPLAY");

			// Setup evaluator and Mads.
			p.set_X0(x0s.at(i));
			p.check();
			Evaluator_ ev(p, this);
			NOMAD::Mads mads(p, &ev);

			// Run Mads.
			mads.run();
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				solutions.push_back(*mads.get_best_feasible());
				objectives.push_back(mads.get_best_feasible()->get_f());
			}
		}
		catch (NOMAD::Exception & e) {
			throw Error(std::string("NOMAD error: ") + e.what());
		}
		NOMAD::Slave::stop_slaves(out);
		NOMAD::end();
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			progress_bar.increment();
			std::cout << "\rNOMAD: " << progress_bar << std::flush;
		}
	}
	std::cout << "\n";

	NOMAD::Double best_objective(objectives.at(0));
	std::size_t pos_best_solution{0};
	for (std::size_t pos{1}; pos < objectives.size(); pos++) {
		if (objectives.at(pos) < best_objective) {
			best_objective = objectives.at(pos);
			pos_best_solution = pos;
		}
	}
	NOMAD::Point best_solution(solutions.at(pos_best_solution));

	// Write the best parameters in alpha_ and lambda_ and write them to output file.
	nomad_point_to_params_(best_solution, alpha_, lambda_);
	normalize_params_();

	if (outfile_ != "") {
		write_params_to_file_();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "[--led-method <arg>] [--sort-method <arg>] [--load <arg>] [--save <arg>] [--init-evaluations <arg>] [--init-initial-solutions <arg>] [--init-mu <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "led-method") {
		if (arg == "LSAPE_OPTIMAL") {
			led_method_ = LSAPE_OPTIMAL;
		}
		else if (arg  == "LSAPE_GREEDY") {
			led_method_ = LSAPE_GREEDY;
		}
		else if (arg == "GAMMA") {
			led_method_ = GAMMA;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option led-method. Usage: options = \"[--led-method LSAPE_OPTIMAL|LSAPE_GREEDY|GAMMA] [...]\"");
		}
		return true;
	}
	else if (option == "sort-method") {
		if (arg == "COUNTING") {
			sort_method_ = COUNTING;
		}
		else if (arg  == "STD") {
			sort_method_ = STD;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option sort-method. Usage: options = \"[--sort-method COUNTING|STD] [...]\"");
		}
		return true;
	}
	else if (option == "load") {
		infile_ = arg;
		return true;
	}
	else if (option == "save") {
		outfile_ = arg;
		return true;
	}
	else if (option == "init-evaluations") {
		try {
			num_evals_ = std::stoi(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option init-evaluations. Usage: options = \"[--init-evaluations <convertible to int greater 0>] [...]");
		}
		if (num_evals_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option init-evaluations. Usage: options = \"[--init-evaluations <convertible to int greater 0>] [...]");
		}
		return true;
	}
	else if (option == "init-initial-solutions") {
		try {
			num_x0s_ = std::stoi(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option init-initial-solutions. Usage: options = \"[--init-initial-solutions <convertible to int greater 0>] [...]");
		}
		if (num_x0s_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option init-initial-solutions. Usage: options = \"[--init-initial-solutions <convertible to int greater 0>] [...]");
		}
		return true;
	}
	else if (option == "init-mu") {
		try {
			mu_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option init-mu. Usage: options = \"[--init-mu <convertible to double between 0 and 1>] [...]");
		}
		if (mu_ < 0 or mu_ > 1) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option init-mu.  Usage: options = \"[--init-mu <convertible to double between 0 and 1>] [...]");
		}
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {
	const NodeRingMap_ & rings_g = rings_.at(g.id());
	const NodeRingMap_ & rings_h = rings_.at(h.id());
#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			master_problem(row_in_master, col_in_master) = compute_ring_distance_(g, h, rings_g, rings_h, alpha_, lambda_, row_in_master, col_in_master);
		}
	}
}

// === Defintion of private class Ring :: Layer_. ===
template<class UserNodeLabel, class UserEdgeLabel>
Ring<UserNodeLabel, UserEdgeLabel>::
Layer_ ::
Layer_(std::size_t level) :
level{level},
node_labels(),
inner_edge_labels(),
outer_edge_labels() {}

// === Definition of private class Ring<UserNodeLabel, UserEdgeLabel>:: Ring_. ===
template<class UserNodeLabel, class UserEdgeLabel>
Ring<UserNodeLabel, UserEdgeLabel>::
Ring_ ::
Ring_() :
layers() {}

// === Definition of private class Ring<UserNodeLabel, UserEdgeLabel>:: Evaluator_. ===
template<class UserNodeLabel, class UserEdgeLabel>
Ring<UserNodeLabel, UserEdgeLabel>::
Evaluator_ ::
Evaluator_(const NOMAD::Parameters & param, Ring<UserNodeLabel, UserEdgeLabel> * ring) :
NOMAD::Evaluator(param),
ring_{ring}{}

template<class UserNodeLabel, class UserEdgeLabel>
Ring<UserNodeLabel, UserEdgeLabel>::
Evaluator_ ::
~Evaluator_(){}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Ring<UserNodeLabel, UserEdgeLabel>::
Evaluator_ ::
eval_x(NOMAD::Eval_Point & x, const NOMAD::Double & h_max, bool & count_eval) const {
	count_eval = true;
	ring_->eval_x_(x);
	return true;
}

// === Definition of helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
nomad_point_to_params_(const NOMAD::Point & x, std::vector<double> & alpha, std::vector<double> & lambda) const {
	alpha.clear();
	for (std::size_t i{0} ; i < 3; i++) {
		alpha.push_back(x.get_coord(i).value());
	}
	std::size_t max_non_zero_level{0};
	for (std::size_t i{0} ; i < num_layers_; i++) {
		if (x.get_coord(i + 3).value() > 0) {
			max_non_zero_level = i;
		}
	}
	lambda.clear();
	for (std::size_t i{0} ; i <= max_non_zero_level; i++) {
		lambda.push_back(x.get_coord(i + 3).value());
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
normalize_params_() {
	double sum_alpha {0.0};
	for (auto alpha = alpha_.begin(); alpha != alpha_.end(); alpha++) {
		sum_alpha += *alpha;
	}
	for (auto alpha = alpha_.begin(); alpha != alpha_.end(); alpha++) {
		*alpha /= sum_alpha;
	}
	double sum_lambda {0.0};
	for (auto lambda = lambda_.begin(); lambda != lambda_.end(); lambda++) {
		sum_lambda += *lambda;
	}
	for (auto lambda = lambda_.begin(); lambda != lambda_.end(); lambda++) {
		*lambda /= sum_lambda;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
eval_x_(NOMAD::Eval_Point & x) const {
	std::vector<double> alpha;
	std::vector<double> lambda;
	nomad_point_to_params_(x, alpha, lambda);
	double val {0.0};
	LSAPESolver lsape_solver;
	lsape_solver.set_model(this->lsape_model_);
	for (auto g = this->ged_data_.begin(); g != this->ged_data_.end(); g++) {
		if (this->ged_data_.is_shuffled_graph_copy(g->id())) {
			continue;
		}
		for (auto h = this->ged_data_.begin(); h != this->ged_data_.end(); h++) {
			if (this->ged_data_.is_shuffled_graph_copy(h->id())) {
				continue;
			}
			NodeMap matching(g->num_nodes(), h->num_nodes());
			DMatrix lsape_instance(g->num_nodes() + 1, h->num_nodes() + 1, 0.0);
			lsape_solver.set_problem(&lsape_instance);
			if (this->ged_data_.shuffled_graph_copies_available() and (g->id() == h->id())) {
				GEDGraph::GraphID id_shuffled_graph_copy{this->ged_data_.id_shuffled_graph_copy(h->id())};
				populate_instance_with_params_(*g, this->ged_data_.graph(id_shuffled_graph_copy), alpha, lambda, lsape_instance);
				lsape_solver.solve();
				util::construct_node_map_from_solver(lsape_solver, matching);
				this->ged_data_.compute_induced_cost(*g, this->ged_data_.graph(id_shuffled_graph_copy), matching);
			}
			else {
				populate_instance_with_params_(*g, *h, alpha, lambda, lsape_instance);
				lsape_solver.solve();
				util::construct_node_map_from_solver(lsape_solver, matching);
				this->ged_data_.compute_induced_cost(*g, *h, matching);
			}
			val += matching.induced_cost();
		}
	}
	val /= static_cast<double>(this->ged_data_.num_graphs() * this->ged_data_.num_graphs());
	std::size_t supp_lambda{0};
	double sum_lambda{0.0};
	for (auto lambda_l = lambda.begin(); lambda_l != lambda.end(); lambda_l++) {
		if (*lambda_l > 0) {
			supp_lambda++;
			sum_lambda += *lambda_l;
		}
	}
	if (num_layers_ > 1 and mu_ < 1) {
		val *=  (mu_ + ((1 - mu_) * static_cast<double>(supp_lambda - 1)  / static_cast<double>(num_layers_ - 1)));
	}
	double sum_alpha{0.0};
	for (auto alpha_i = alpha.begin(); alpha_i != alpha.end(); alpha_i++) {
		sum_alpha += *alpha_i;
	}
	double epsilon{0.00001};
	x.set_bb_output(0, val);
	x.set_bb_output(1, epsilon - sum_alpha);
	x.set_bb_output(2, epsilon - sum_lambda);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
build_rings_(const GEDGraph & graph) {
	rings_[graph.id()] = NodeRingMap_();
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		build_ring_(graph, *node, rings_.at(graph.id()));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
build_ring_(const GEDGraph & graph, GEDGraph::NodeID root, NodeRingMap_ & rings) {
	std::map<GEDGraph::NodeID, int> distance_to_root;
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		distance_to_root[*node] = -1;
	}
	distance_to_root[root] = 0;

	std::map<GEDGraph::EdgeID, bool> discovered_edge;
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		discovered_edge[*edge] = false;
	}

	Layer_ current_layer(0);
	std::queue<GEDGraph::NodeID> queue;
	queue.push(root);
	while (not queue.empty()) {
		GEDGraph::NodeID current_node{queue.front()};
		queue.pop();
		if (static_cast<int>(current_layer.level) < distance_to_root.at(current_node)) {
			if (sort_method_ == COUNTING) {
				util::counting_sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
				util::counting_sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
				util::counting_sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
			}
			else {
				std::sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
				std::sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
				std::sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
			}
			rings[root].layers.push_back(current_layer);
			current_layer = Layer_(current_layer.level + 1);
		}
		current_layer.node_labels.push_back(graph.get_node_label(current_node));
		for (auto edge = graph.incident_edges(current_node).first; edge != graph.incident_edges(current_node).second; edge++) {
			GEDGraph::NodeID next_node{graph.head(*edge)};
			if (distance_to_root.at(next_node) == -1) {
				distance_to_root[next_node] = current_layer.level + 1;
				if (current_layer.level < num_layers_) {
					queue.push(next_node);
				}
			}
			if (not discovered_edge.at(*edge)) {
				discovered_edge[*edge] = true;
				if (distance_to_root.at(current_node) == distance_to_root.at(next_node)) {
					current_layer.inner_edge_labels.push_back(graph.get_edge_label(*edge));
				}
				else if (distance_to_root.at(current_node) < distance_to_root.at(next_node)) {
					current_layer.outer_edge_labels.push_back(graph.get_edge_label(*edge));
				}
				else {
					throw Error(std::string("Error when building ring rooted at ") + std::to_string(root) +
							" for graph " + std::to_string(graph.id()) + ": dist(" +
							std::to_string(current_node) +") = " + std::to_string(distance_to_root.at(current_node)) +
							" > dist(" + std::to_string(next_node) +") = " + std::to_string(distance_to_root.at(next_node)));
				}
			}
		}
	}
	if (led_method_ == GAMMA) {
		if (sort_method_ == COUNTING) {
			util::counting_sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
			util::counting_sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
			util::counting_sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
		}
		else {
			std::sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
			std::sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
			std::sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
		}
	}
	rings[root].layers.push_back(current_layer);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
set_num_layers_() {
	std::size_t max_num_layers{0};
	for (auto graph_rings_pair = rings_.begin(); graph_rings_pair != rings_.end(); graph_rings_pair++) {
		for (auto ring = graph_rings_pair->second.begin(); ring != graph_rings_pair->second.end(); ring++) {
			max_num_layers = std::max(max_num_layers, ring->second.layers.size());
		}
	}
	num_layers_ = std::min(num_layers_, max_num_layers);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Ring<UserNodeLabel, UserEdgeLabel>::
load_config_file_() const {
	return (infile_ != "");
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
init_x0s_(std::vector<NOMAD::Point> & x0s) const {
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> uni(0.0, 1.0);
	uni(rng);
	NOMAD::Point x0(3 + num_layers_, 0.0);
	bool found_new_x0{true};
	std::vector<double> alpha_gen;
	std::vector<double> lambda_gen;
	while(x0s.size() < num_x0s_) {
		alpha_gen.clear();
		alpha_gen.push_back(0.0);
		alpha_gen.push_back(1.0);
		for (std::size_t i{0}; i < 2; i++) {
			alpha_gen.push_back(uni(rng));
		}
		std::sort(alpha_gen.begin(), alpha_gen.end());
		for (std::size_t i{0}; i < 3; i++) {
			x0[i] = alpha_gen.at(i+1) - alpha_gen.at(i);
		}
		lambda_gen.clear();
		lambda_gen.push_back(0.0);
		lambda_gen.push_back(1.0);
		for (std::size_t level{0}; level < num_layers_ - 1; level++) {
			lambda_gen.push_back(uni(rng));
		}
		std::sort(lambda_gen.begin(), lambda_gen.end());
		for (std::size_t level{0}; level < num_layers_; level++) {
			x0[level + 3] = lambda_gen.at(level + 1) - lambda_gen.at(level);
		}
		for (auto old_x0 = x0s.begin(); old_x0 != x0s.end(); old_x0++) {
			if (*old_x0 == x0) {
				found_new_x0 = false;
				break;
			}
		}
		if (found_new_x0) {
			x0s.push_back(x0);
		}
		found_new_x0 = true;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Ring<UserNodeLabel, UserEdgeLabel>::
compute_ring_distance_(const GEDGraph & g, const GEDGraph & h, const NodeRingMap_ & rings_g, const NodeRingMap_ & rings_h,
		const std::vector<double> & alpha, const std::vector<double> & lambda, std::size_t row_in_master, std::size_t col_in_master) const {
	double red{0.0};
	if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) { // compute substitution cost
		const Ring_ & ring_i = rings_g.at(row_in_master);
		const Ring_ & ring_k = rings_h.at(col_in_master);
		for (std::size_t level{0}; level < lambda.size(); level++) {
			red += lambda.at(level) * compute_substitution_cost_(ring_i, ring_k, alpha, level);
		}
	}
	else if (row_in_master < g.num_nodes()) { // compute deletion cost
		const Ring_ & ring_i = rings_g.at(row_in_master);
		for (std::size_t level{0}; level < lambda.size(); level++) {
			red += lambda.at(level) * compute_deletion_cost_(ring_i, alpha, level);
		}
	}
	else if (col_in_master < h.num_nodes()) { // compute insertion cost
		const Ring_ & ring_k = rings_h.at(col_in_master);
		for (std::size_t level{0}; level < lambda.size(); level++) {
			red += lambda.at(level) * compute_insertion_cost_(ring_k, alpha, level);
		}
	}
	return red;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Ring<UserNodeLabel, UserEdgeLabel>::
compute_substitution_cost_(const Ring_ & ring_i, const Ring_ & ring_k, const std::vector<double> & alpha, std::size_t level) const {
	if ((ring_i.layers.size() > level) and (ring_k.layers.size() > level)) {
		return compute_layer_distance_(ring_i.layers.at(level), ring_k.layers.at(level), alpha);
	}
	if (ring_i.layers.size() > level) {
		return compute_layer_distance_(ring_i.layers.at(level), Layer_(0), alpha);
	}
	if (ring_k.layers.size() > level) {
		return compute_layer_distance_(Layer_(0), ring_k.layers.at(level), alpha);
	}
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Ring<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const Ring_ & ring,  const std::vector<double> & alpha, std::size_t level) const {
	if (ring.layers.size() > level) {
		return compute_layer_distance_(ring.layers.at(level), Layer_(0), alpha);
	}
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Ring<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const Ring_ & ring, const std::vector<double> & alpha, std::size_t level) const {
	if (ring.layers.size() > level) {
		return compute_layer_distance_(Layer_(0), ring.layers.at(level), alpha);
	}
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Ring<UserNodeLabel, UserEdgeLabel>::
compute_layer_distance_(const Layer_ & lhs, const Layer_ & rhs, const std::vector<double> & alpha) const {
	double node_cost{0.0};
	double inner_edge_cost{0.0};
	double outer_edge_cost{0.0};

	switch (led_method_) {
	case GAMMA:
		node_cost = gamma_multiset_cost_(lhs.node_labels, rhs.node_labels, true);
		inner_edge_cost = gamma_multiset_cost_(lhs.inner_edge_labels, rhs.inner_edge_labels, false);
		outer_edge_cost = gamma_multiset_cost_(lhs.outer_edge_labels, rhs.outer_edge_labels, false);
		break;
	default:
		node_cost = lsape_multiset_cost_(lhs.node_labels, rhs.node_labels, true);
		inner_edge_cost = lsape_multiset_cost_(lhs.inner_edge_labels, rhs.inner_edge_labels, false);
		outer_edge_cost = lsape_multiset_cost_(lhs.outer_edge_labels, rhs.outer_edge_labels, false);
		break;
	}

	std::size_t max_num_node_labels{std::max(lhs.node_labels.size(), rhs.node_labels.size())};
	if (max_num_node_labels > 0) {
		node_cost /= static_cast<double>(max_num_node_labels);
	}
	std::size_t max_num_inner_edge_labels{std::max(lhs.inner_edge_labels.size(), rhs.inner_edge_labels.size())};
	if (max_num_inner_edge_labels > 0) {
		inner_edge_cost /= static_cast<double>(max_num_inner_edge_labels);
	}
	std::size_t max_num_outer_edge_labels{std::max(lhs.outer_edge_labels.size(), rhs.outer_edge_labels.size())};
	if (max_num_outer_edge_labels > 0) {
		outer_edge_cost /= static_cast<double>(max_num_outer_edge_labels);
	}

	return alpha.at(0) * node_cost + alpha.at(1) * inner_edge_cost + alpha.at(2) * outer_edge_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Ring<UserNodeLabel, UserEdgeLabel>::
lsape_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const {

	if ((lhs.size() == 0) and (rhs.size() == 0)) {
		return 0.0;
	}

	if ((lhs.size() > 0) and (rhs.size() == 0)) {
		double cost{0.0};
		for (std::size_t row{0}; row < lhs.size(); row++) {
			if (node_labels) {
				cost += this->ged_data_.node_cost(lhs.at(row), dummy_label());
			}
			else {
				cost += this->ged_data_.edge_cost(lhs.at(row), dummy_label());
			}
		}
		return cost;
	}

	if ((lhs.size() == 0) and (rhs.size() > 0)) {
		double cost{0.0};
		for (std::size_t col{0}; col < rhs.size(); col++) {
			if (node_labels) {
				cost += this->ged_data_.node_cost(dummy_label(), rhs.at(col));
			}
			else {
				cost += this->ged_data_.edge_cost(dummy_label(), rhs.at(col));
			}
		}

		return cost;
	}

	DMatrix problem(lhs.size() + 1, rhs.size() + 1, 0.0);
	// Collect deletion costs.
	for (std::size_t row{0}; row < lhs.size(); row++) {
		if (node_labels) {
			problem(row, rhs.size()) = this->ged_data_.node_cost(lhs.at(row), dummy_label());
		}
		else {
			problem(row, rhs.size()) = this->ged_data_.edge_cost(lhs.at(row), dummy_label());
		}
	}

	// Collect insertion costs.
	for (std::size_t col{0}; col < rhs.size(); col++) {
		if (node_labels) {
			problem(lhs.size(), col) = this->ged_data_.node_cost(dummy_label(), rhs.at(col));
		}
		else {
			problem(lhs.size(), col) = this->ged_data_.edge_cost(dummy_label(), rhs.at(col));
		}
	}

	// Collect substitution costs.
	for (std::size_t row{0}; row < lhs.size(); row++) {
		for (std::size_t col{0}; col < rhs.size(); col++) {
			if (node_labels) {
				problem(row, col) = this->ged_data_.node_cost(lhs.at(row), rhs.at(col));
			}
			else {
				problem(row, col) = this->ged_data_.edge_cost(lhs.at(row), rhs.at(col));
			}
		}
	}

	LSAPESolver problem_solver(&problem);
	if (led_method_ == LSAPE_OPTIMAL) {
		problem_solver.set_model(this->lsape_model_);
	}
	else {
		problem_solver.set_greedy_method(this->greedy_method_);
	}
	problem_solver.solve();

	return problem_solver.minimal_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Ring<UserNodeLabel, UserEdgeLabel>::
gamma_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const {
	double cost{0.0};

	// Compute and add minimal edge insertion costs.
	if (lhs.size() < rhs.size()) {
		double avg_ins_cost{0.0};
		for (auto label_rhs = rhs.begin(); label_rhs != rhs.end(); label_rhs++) {
			if (node_labels) {
				avg_ins_cost += this->ged_data_.node_cost(dummy_label(), *label_rhs);
			}
			else {
				avg_ins_cost += this->ged_data_.edge_cost(dummy_label(), *label_rhs);
			}
		}
		avg_ins_cost /= static_cast<double>(rhs.size());
		cost += static_cast<double>(rhs.size() - lhs.size()) * avg_ins_cost;
	}

	// Compute and add minimal edge deletion costs.
	if (lhs.size() > rhs.size()) {
		double avg_del_cost{0.0};
		for (auto label_lhs = lhs.begin(); label_lhs != lhs.end(); label_lhs++) {
			if (node_labels) {
				avg_del_cost += this->ged_data_.node_cost(*label_lhs, dummy_label());
			}
			else {
				avg_del_cost += this->ged_data_.edge_cost(*label_lhs, dummy_label());
			}
		}
		avg_del_cost /= static_cast<double>(lhs.size());
		cost += static_cast<double>(lhs.size() - rhs.size()) * avg_del_cost;
	}

	// Compute minimal edge relabelling costs.
	double avg_rel_cost{0.0};
	std::size_t count_diff_labels{0};
	for (auto label_lhs = lhs.begin(); label_lhs != lhs.end(); label_lhs++) {
		for (auto label_rhs = rhs.begin(); label_rhs != rhs.end(); label_rhs++) {
			if (*label_lhs != *label_rhs) {
				count_diff_labels++;
				if (node_labels) {
					avg_rel_cost += this->ged_data_.node_cost(*label_lhs, *label_rhs);
				}
				else {
					avg_rel_cost += this->ged_data_.edge_cost(*label_lhs, *label_rhs);
				}
			}
		}
	}
	avg_rel_cost /= static_cast<double>(count_diff_labels);

	// Compute multiset intersection size.
	std::size_t intersection_size{0};
	auto label_lhs = lhs.begin();
	auto label_rhs = rhs.begin();
	while ((label_lhs != lhs.end()) and (label_rhs != rhs.end())) {
		if (*label_lhs == *label_rhs) {
			intersection_size++;
			label_lhs++;
			label_rhs++;
		}
		else if (*label_lhs < *label_rhs) {
			label_lhs++;
		}
		else {
			label_rhs++;
		}
	}

	std::size_t gamma(std::min(lhs.size(), rhs.size()) - intersection_size);
	if (gamma > 0) {
		cost += static_cast<double>(gamma) * avg_rel_cost;
	}

	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
populate_instance_with_params_(const GEDGraph & g, const GEDGraph & h, const vector<double> & alpha, const vector<double> & lambda, DMatrix & lsape_instance) const {
	const NodeRingMap_ & rings_g = rings_.at(g.id());
	const NodeRingMap_ & rings_h = rings_.at(h.id());

	for (std::size_t row_in_master = 0; row_in_master < lsape_instance.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < lsape_instance.num_cols(); col_in_master++) {
			lsape_instance(row_in_master, col_in_master) = compute_ring_distance_(g, h, rings_g, rings_h, alpha, lambda, row_in_master, col_in_master);
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
write_params_to_file_() const {
	std::map<std::string, std::string> options;
	options["num_layers"] = std::to_string(lambda_.size());
	for (std::size_t i{0}; i < alpha_.size(); i++) {
		double val{alpha_.at(i)};
		options["alpha_" + std::to_string(i)] = std::to_string(val);
	}
	for (std::size_t i{0}; i <lambda_.size(); i++) {
		double val{lambda_.at(i)};
		options["lambda_" + std::to_string(i)] = std::to_string(val);
	}
	util::save_as_config_file(outfile_, options);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Ring<UserNodeLabel, UserEdgeLabel>::
read_params_from_file_() {
	std::map<std::string, std::string> options;
	util::parse_config_file(infile_, options);
	num_layers_ = std::stoul(options.at("num_layers"));
	alpha_.clear();
	for (std::size_t i{0}; i < 3; i++) {
		alpha_.push_back(std::stod(options.at("alpha_" + std::to_string(i))));
	}
	lambda_.clear();
	for (std::size_t i{0}; i < num_layers_; i++) {
		lambda_.push_back(std::stod(options.at("lambda_" + std::to_string(i))));
	}
}

}

#endif /* SRC_METHODS_RING_IPP_ */
