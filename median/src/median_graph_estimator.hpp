/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2019 by David B. Blumenthal                              *
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
 * @file median_graph_estimator.hpp
 * @brief ged::MedianGraphEstimator class declaration.
 */

#ifndef MEDIAN_SRC_MEDIAN_GRAPH_ESTIMATOR_HPP_
#define MEDIAN_SRC_MEDIAN_GRAPH_ESTIMATOR_HPP_

#include "../../src/env/ged_env.hpp"

namespace ged {

/*!
 * @brief Class for estimating generalized median graphs.
 *
 * @details Implements the algorithm suggested in:
 * - N. Boria, S. Bougleux, B. Ga&uuml;z&egrave;re, D. B. Blumenthal, and L. Brun:
 *   &ldquo;Scalable generalized graph median estimation and applications in clustering, classification, and indexing&rdquo;
 *   submitted to VLDB J.,
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--init-type RANDOM\|MEDOID\|MIN\|MAX\|MEAN</tt> | method for computing the initial medians | @p RANDOM | unless @p RANDOM, the option @p \--random-inits has no effect |
 * | <tt>\--update-order TRUE\|FALSE</tt> | update the order of the medians during optimization | @p TRUE | should always be set to TRUE, we included the option only for test purposes |
 * | <tt>\--random-inits @<convertible to int greater 0@></tt> | number of randomly constructed initial medians | @p 50 | n.a. |
 * | <tt>\--randomness REAL\|PSEUDO</tt> | use real randomness or pseudo randomness | @p REAL | if @p REAL, the option @p \--seed has no effect |
 * | <tt>\--seed @<convertible to int greater equal 0@></tt> | seed for generating pseudo random numbers | @p 0 | n.a. |
 * | <tt>\--refine TRUE\|FALSE</tt> | improve node maps and sum of distances for converged median | @p TRUE | n.a. |
 * | <tt>\--max-itrs @<convertible to int@></tt> | maximal number of iterations in main block gradient descent | @p 100 | if negative, no maximal number of iterations is enforced |
 * | <tt>\--max-itrs-without-update @<convertible to int@></tt> | maximal number of consecutive iterations in main block gradient descent where the median is not updated | @p 3 | if negative, no maximal number of iterations without update is enforced |
 * | <tt>\--time-limit @<convertible to double@></tt> | time limit in seconds for main block gradient descent | @p 0 | if less or equal @p 0, no time limit is enforced |
 * | <tt>\--epsilon @<convertible to double greater 0@></tt> | convergence threshold used everywhere | @p 0.0001 | n.a. |
 * | <tt>\--inits-increase-order @<convertible to int greater 0@></tt> | number of initial solutions for generic heuristic to increase the order of the median | @p 5 | used as starting points for (parallel) block gradient descents to determine node label of inserted node |
 * | <tt>\--init-type-increase-order CLUSTERS\|K-MEANS++</tt> | initialization type for generic heuristic to increase the order of the median | @p K-MEANS++ | if @p K-MEANS++, well distributed node labels are used as starting points |
 * | <tt>\--max-itrs-increase-order @<convertible to int@></tt> | maximal number of iterations used in generic heuristic to increase the order of the median | @p 10 |  if negative, no iteration based termination criterion is used |
 * | <tt>\--stdout 0\|1\|2</tt> | print runtime information to standard output stream | @p 2 | @p 0: no output; @p 1: output only before termination; @p 2: output also during optimization |
 */
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
class MedianGraphEstimator {

public:

	/*!
	 * @brief Constructor.
	 * @param[in,out] ged_env Pointer to initialized environment. The edit costs must be set by the user.
	 * @param[in] constant_node_costs Set to @p true if the node relabeling costs are constant.
	 */
	MedianGraphEstimator(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, bool constant_node_costs);

	/*!
	 * @brief Sets the options of the estimator.
	 * @param[in] options String that specifies with which options to run the estimator.
	 */
	void set_options(const std::string & options);

	/*!
	 * @brief Selects method to be used for computing the initial medoid graph.
	 * @param[in] init_method The selected method. Default: ged::Options::GEDMethod::BRANCH_UNIFORM.
	 * @param[in] init_options The options for the selected method. Default: "".
	 * @note Has no effect unless "--init-type MEDOID" is passed to set_options().
	 */
	void set_init_method(Options::GEDMethod init_method, const std::string & init_options = "");

	/*!
	 * @brief Selects method to be used for block gradient descent..
	 * @param[in] descent_method The selected method. Default: ged::Options::GEDMethod::BRANCH_FAST.
	 * @param[in] descent_options The options for the selected method. Default: "".
	 * @note Has no effect unless "--init-type MEDOID" is passed to set_options().
	 */
	void set_descent_method(Options::GEDMethod descent_method, const std::string & descent_options = "");

	/*!
	 * @brief Selects method to be used for improving the sum of distances and the node maps for the converged median.
	 * @param[in] refine_method The selected method. Default: ged::Options::GEDMethod::IPFP.
	 * @param[in] refine_options The options for the selected method. Default: "".
	 * @note Has no effect if "--refine FALSE" is passed to set_options().
	 */
	void set_refine_method(Options::GEDMethod refine_method, const std::string & refine_options = "");

	/*!
	 * @brief Computes a generalized median graph.
	 * @param[in] graph_ids The IDs of the graphs for which the median should be computed.
	 * Must have been added to the environment passed to the constructor.
	 * @param[in] median_id The ID of the computed median.
	 * A dummy graph with this ID must have been added to the environment passed to the constructor.
	 * Upon termination, the computed median can be obtained via ged::GEDEnv::get_graph().
	 */
	void run(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID median_id);

	/*!
	 * @brief Returns the state of the estimator.
	 * @return The state of the estimator at termination of the last call to run().
	 * If run without time limit or maximum number of iterations, this method always returns ged::Options::AlgorithmState::TERMINATED.
	 * Otherwise, the last uninterrupted state is returned.
	 */
	Options::AlgorithmState get_state() const;

	/*!
	 * @brief Returns the sum of distances.
	 * @param[in] state The state of the estimator.
	 * @return The sum of distances of the median when the estimator was in the state @p state during the last call to run().
	 */
	double get_sum_of_distances(Options::AlgorithmState state = Options::AlgorithmState::TERMINATED) const;

	/*!
	 * @brief Returns distance from the median.
	 * @param[in] graph_id ID of the graph whose distance from the median should be returned.
	 * Must have been contained in the collection of IDs passed to run().
	 * @return Distance of the graph with ID @p graph_id from the median upon termination.
	 */
	double get_distance_from_median(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Returns node map from the median.
	 * @param[in] graph_id ID of the graph whose node map from the median should be returned.
	 * Must have been contained in the collection of IDs passed to run().
	 * @return Node map from the median to the graph with ID @p graph_id upon termination.
	 */
	const NodeMap & get_node_map_from_median(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Computes node map from the median by using the GED method employed to obtain the final node maps.
	 * @param graph_id ID of the graph whose node map from the median should be returned.
	 * Must not necessarily have been contained in the collection of IDs passed to run().
	 * @return  Node map from the median to the graph with ID @p graph_id.
	 */
	const NodeMap & compute_node_map_from_median(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Returns the runtime.
	 * @param[in] state The state of the estimator.
	 * @return The runtime up to the point where the estimator entered the state @p state during the last call to run().
	 */
	double get_runtime(Options::AlgorithmState state = Options::AlgorithmState::TERMINATED) const;

	/*!
	 * @brief Returns number of iterations.
	 * @return A vector that contains the number of iterations for each initial median for the last call to run().
	 */
	const std::vector<std::size_t> & get_num_itrs() const;

	/*!
	 * @brief Returns the number of times the order of the median decreased.
	 * @return Overall number of times the order of the median decreased during the last call to run().
	 */
	std::size_t get_num_times_order_decreased() const;

	/*!
	 * @brief Returns the number of times the order of the median increased.
	 * @return Overall number of times the order of the median increased during the last call to run().
	 */
	std::size_t get_num_times_order_increased() const;

	/*!
	 * @brief Returns the number of converged descents.
	 * @return Overall number of converged descents during the last call to run().
	 */
	std::size_t get_num_converged_descents() const;

	/*!
	 * @brief Returns pointer to the environment employed by the estimator.
	 * @return Pointer to the environment employed by the estimator.
	 */
	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * get_ged_env();

private:

	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env_;

	Options::GEDMethod init_method_;

	std::string init_options_;

	Options::GEDMethod descent_method_;

	std::string descent_options_;

	Options::GEDMethod refine_method_;

	std::string refine_options_;

	bool constant_node_costs_;

	bool labeled_nodes_;

	double node_del_cost_;

	double node_ins_cost_;

	bool labeled_edges_;

	double edge_del_cost_;

	double edge_ins_cost_;

	std::string init_type_;

	bool update_order_;

	std::size_t num_random_inits_;

	std::size_t desired_num_random_inits_;

	bool use_real_randomness_;

	std::size_t seed_;

	bool refine_;

	double time_limit_in_sec_;

	double epsilon_;

	std::size_t max_itrs_;

	std::size_t max_itrs_without_update_;

	std::size_t num_inits_increase_order_;

	std::string init_type_increase_order_;

	std::size_t max_itrs_increase_order_;

	std::size_t print_to_stdout_;

	GEDGraph::GraphID median_id_;

	std::map<GEDGraph::GraphID, NodeMap> node_maps_from_median_;

	double sum_of_distances_;

	double best_init_sum_of_distances_;

	double converged_sum_of_distances_;

	Seconds runtime_;

	Seconds runtime_initialized_;

	Seconds runtime_converged_;

	std::vector<std::size_t> itrs_;

	std::size_t num_decrease_order_;

	std::size_t num_increase_order_;

	std::size_t num_converged_descents_;

	Options::AlgorithmState state_;

	void set_default_options_();

	void construct_initial_medians_(const std::vector<GEDGraph::GraphID> & graph_ids, const Timer & timer, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians);

	void compute_medoid_(const std::vector<GEDGraph::GraphID> & graph_ids, const Timer & timer, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians);

	void compute_max_order_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) const;

	void compute_min_order_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) const;

	void compute_mean_order_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) const;

	void sample_initial_medians_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians);

	bool termination_criterion_met_(bool converged, const Timer & timer, std::size_t itr, std::size_t itrs_without_update);

	bool update_median_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) const;

	void update_node_labels_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) const;

	void update_edges_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) const;

	bool update_node_maps_();

	bool decrease_order_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median);

	double compute_best_deletion_delta_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median, std::size_t & id_deleted_node) const;

	void delete_node_from_median_(std::size_t id_deleted_node, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median);

	bool increase_order_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median);

	double compute_best_insertion_delta_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, 
			std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const;

	double compute_insertion_delta_unlabeled_(const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes,
			std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const;

	double compute_insertion_delta_constant_(const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes,
			std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const;

	double compute_insertion_delta_generic_(const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes,
			std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const;

	void compute_initial_node_labels_(const std::vector<UserNodeLabel> & node_labels, std::vector<UserNodeLabel> & median_labels_) const;

	bool insertion_termination_criterion_met_(bool converged, std::size_t itr) const;

	bool update_config_(const UserNodeLabel & node_label, const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes, std::map<GEDGraph::GraphID, std::pair<std::size_t, UserNodeLabel>> & config, std::vector<UserNodeLabel> & node_labels) const;

	bool update_node_label_(const std::vector<UserNodeLabel> & node_labels, UserNodeLabel & node_label) const;

	bool update_clusters_(const std::vector<UserNodeLabel> & node_labels, const std::vector<UserNodeLabel> & median_labels, std::vector<std::size_t> & closest_median_ids) const;

	void add_node_to_median_(const std::map<GEDGraph::GraphID, std::size_t> & best_config, const UserNodeLabel & best_label, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median);

	void improve_sum_of_distances_(const Timer & timer);

	bool median_available_() const;

};

}

#include "median_graph_estimator.ipp"

#endif /* MEDIAN_SRC_MEDIAN_GRAPH_ESTIMATOR_HPP_ */
