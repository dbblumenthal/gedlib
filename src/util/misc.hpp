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
 * @param[out] node_map The constructed node map.
 * @param[in] solution_id The ID of the solution stored in @p solver from which the matching should be constructed.
 */
template<class Solver>
void construct_node_map_from_solver(const Solver & solver, NodeMap & node_map, std::size_t solution_id = 0);

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
 * @param[out] adj_matrix The constructed adjacency matrix.
 */
void init_adj_matrix(const GEDGraph & graph, DMatrix & adj_matrix);

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

/*!
 * @brief Transforms an options string into an options map.
 * @param[in] options_string Options string of the form "[--<option> <arg>] [...]".
 * @param[out] options_map Map with one key-value pair (<option>,<arg>) for each option contained in the string.
 */
void options_string_to_options_map(const std::string & options_string, std::map<std::string, std::string> & options_map);

/*!
 * @brief Separates a sentence into words separated by @p sep (unless contained in single quotes).
 * @param[in] sentence The sentence that should be tokenized.
 * @param[in] sep The separator. Must be different from "'".
 * @param[out] words The obtained words.
 */
void tokenize(const std::string & sentence, char sep, std::vector<std::string> & words);

/*!
 * @brief Checks whether a word is an option name and, if so, removes the leading dashes.
 * @param[in,out] word
 * @return True if @p word is of the form "--<option>".
 */
bool is_option_name(std::string & word);

/*!
 * @brief Returns a string or numeric value not contained in a given vector of values.
 * @tparam StringOrNumeric Either std::string or a numeric type.
 * @param[in] values Vector of values that should be different from returned value.
 * @return A value which is different from all elements of @p values.
 */
template<class StringOrNumeric>
StringOrNumeric new_string_or_numeric_(const std::vector<StringOrNumeric> & values);

}

}

#include "misc.ipp"

#endif /* SRC_UTIL_MISC_HPP_ */
