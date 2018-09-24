/*!
 * 	@file  misc.hpp
 *  @brief Declaration of miscellaneous utility functions.
 */

#ifndef SRC_UTIL_MISC_HPP_
#define SRC_UTIL_MISC_HPP_

#include "../env/ged_graph.hpp"
#include "../env/node_map.hpp"
#include "../util/lsap_solver.hpp"
#include "../util/lsape_solver.hpp"

namespace ged {

/*!
 * @namespace ged::util
 * @brief Contains miscellaneous utility functions.
 */
namespace util {

/*!
 * @brief Constructs a node map from a solution to LSAPE or LSAPE stored in a ged::LSAPESolver or a ged::LSAPSolver object.
 * @tparam Solver Must be either ged::LSAPESolver or ged::LSAPSolver.
 * @param[in] solver Solver object that has been solved by call to ged::LSAPESolver::solve() or ged::LSAPSolver::solve(), respectively.
 * @param[in] g_ids_to_nodes All keys must be rows in @p solver.
 * @param[in] h_ids_to_nodes All keys must be columns in @p solver.
 * @param[out] node_map The constructed node map.
 * @param[in] solution_id The ID of the solution stored in @p solver from which the matching should be constructed.
 */
template<class Solver>
void construct_node_map_from_solver(const Solver & solver, const GEDGraph::SizeTNodeMap & g_ids_to_nodes,
		const GEDGraph::SizeTNodeMap & h_ids_to_nodes, NodeMap & node_map, std::size_t solution_id = 0);

/*!
 * @brief Initializes an index from node IDs to integer IDs and stores it in a map with the graph's ID as its key.
 * @param[in] graph Input graph whose nodes should be indexed.
 * @param[out] nodes_to_ids Map where the constructed index is stored.
 */
void init_node_to_id_indices(const GEDGraph & graph, std::map<GEDGraph::GraphID, GEDGraph::NodeSizeTMap> & nodes_to_ids);

/*!
 * @brief Initializes an index from node IDs to integer IDs.
 * @param[in] graph Input graph whose nodes should be indexed.
 * @param[out] nodes_to_ids The constructed index.
 */
void init_node_to_id_indices(const GEDGraph & graph, GEDGraph::NodeSizeTMap & nodes_to_ids);

/*!
 * @brief Initializes an index from integer IDs to node IDs and stores it in a map with the graph's ID as its key.
 * @param[in] graph Input graph whose nodes should be indexed.
 * @param[out] ids_to_nodes Map where the constructed index is stored.
 */
void init_id_to_node_indices(const GEDGraph & graph, std::map<GEDGraph::GraphID, GEDGraph::SizeTNodeMap> & ids_to_nodes);

/*!
 * @brief Initializes an index from integer IDs to node IDs.
 * @param[in] graph Input graph whose nodes should be indexed.
 * @param[out] ids_to_nodes The constructed index.
 */
void init_id_to_node_indices(const GEDGraph & graph, GEDGraph::SizeTNodeMap & ids_to_nodes);

/*!
 * @brief Implementation of counting sort.
 * @param[in] first Iterator to first element in subvector that should be sorted.
 * @param[in] last Iterator to last element in subvector that should be sorted.
 * @see https://en.wikipedia.org/wiki/Counting_sort
 */
void counting_sort(std::vector<LabelID>::iterator first, std::vector<LabelID>::iterator last);

/*!
 * @brief Initalizes the adjacency matrix of a graph.
 * @param[in] graph The graph whose adjacency matrix should be constructed.
 * @param[in] ids_to_nodes An index that maps integer IDs to the IDs of the nodes in @p graph.
 * @param[out] adj_matrix The constructed adjacency matrix.
 */
void init_adj_matrix(const GEDGraph & graph, const GEDGraph::SizeTNodeMap & ids_to_nodes, DMatrix & adj_matrix);

/*!
 * @brief Parses a configuration file.
 * @param[in] filename Name of a file whose lines are of the form <tt>@<key@>=@<value@></tt>.
 * @param[out] options String map that contains the file's keys as keys and the file's values as values.
 */
void parse_config_file(const std::string & filename, std::map<std::string, std::string> & options);

/*!
 * @brief Saves a string map as a configuration file as expected by parse_config_file().
 * @param[in] filename Name of the configuration file that should be created.
 * @param[in] options String map that should be saved as a configuration file.
 */
void save_as_config_file(const std::string & filename, const std::map<std::string, std::string> & options);

}

}

#include "misc.ipp"

#endif /* SRC_UTIL_MISC_HPP_ */
