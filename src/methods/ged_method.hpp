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
 * @file  ged_method.hpp
 * @brief ged::GEDMethod class declaration.
 */

#ifndef SRC_METHODS_GED_METHOD_HPP_
#define SRC_METHODS_GED_METHOD_HPP_

#include "../env/common_types.hpp"
#include "../env/ged_data.hpp"
#include "../util/misc.hpp"

namespace ged {

/*!
 * \brief Abstract class for the (suboptimal) computation of the graph edit distance.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class GEDMethod {

public:

	/*!
	 * @brief Pure virtual destructor.
	 * @note Must be implemented by derived classes.
	 */
	virtual ~GEDMethod() = 0;

	/*!
	 * @brief Constructor.
	 * @param[in] ged_data The instance on which the method should be run.
	 */
	GEDMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

	/*!
	 * @brief Sets the options of the method.
	 * @param[in] options String of the form <tt>[--@<option@> @<arg@>] [...]</tt>, where @p option contains neither spaces nor single quotes,
	 * and @p arg contains neither spaces nor single quotes or is of the form <tt>'[--@<sub-option@> @<sub-arg@>] [...]'</tt>,
	 * where both @p sub-option and @p sub-arg contain neither spaces nor single quotes.
	 */
	void set_options(const std::string & options);

	/*!
	 * @brief Runs the method with options specified by set_options().
	 * @param[in] g_id ID of input graph.
	 * @param[in] h_id ID of input graph.
	 */
	void run(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id);

	/*!
	 * @brief Runs the method with options specified by set_options().
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[out] result Result variable.
	 */
	void run_as_util(const GEDGraph & g, const GEDGraph & h, Result & result);

	/*!
	 * @brief Initializes the method with options specified by set_options().
	 */
	void init();

	/*!
	 * @brief Returns an upper bound.
	 * @return Upper bound for graph edit distance provided by last call to run() or -1 if the method does not yield an upper bound.
	 */
	double get_upper_bound() const;

	/*!
	 * @brief Returns a lower bound.
	 * @return Lower bound for graph edit distance provided by last call to run() or -1 if the method does not yield a lower bound.
	 */
	double get_lower_bound() const;

	/*!
	 * @brief Returns the runtime.
	 * @return Runtime of last call to run() in seconds.
	 */
	Seconds get_runtime() const;

	/*!
	 * @brief Returns the initialization time.
	 * @return Runtime of last call to init() in seconds.
	 */
	Seconds get_init_time() const;

	/*!
	 * @brief Returns a graph matching.
	 * @return Constant reference to graph matching provided by last call to run() or to an empty matching if the method does not yield a matching.
	 */
	const NodeMap & get_node_map() const;

protected:

	/*!
	 * @brief A flag that equals true if init() has been called and false otherwise.
	 */
	bool initialized_;

	/*!
	 * @brief The data on which the method is run.
	 */
	const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data_;

private:

	std::map<std::string, std::string> options_;

	double lower_bound_;

	double upper_bound_;

	NodeMap node_map_;

	Seconds runtime_;

	Seconds init_time_;

	// Virtual member functions to be overridden by derived classes.

	/*!
	 * @brief Initializes the method.
	 * @note Must be overridden by derived classes that require initialization.
	 */
	virtual void ged_init_();

	/*!
	 * @brief Runs the method with options specified by set_options().
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[out] result Result variable.
	 * @note Must be overridden by derived classes.
	 */
	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result);

	/*!
	 * @brief Parses one option.
	 * @param[in] option The name of the option.
	 * @param[in] arg The argument of the option.
	 * @return Boolean @p true if @p option is a valid option name for the method and @p false otherwise.
	 * @note Must be overridden by derived classes that have options.
	 */
	virtual bool ged_parse_option_(const std::string & option, const std::string & arg);

	/*!
	 * @brief Returns string of all valid options.
	 * @return String of the form <tt>[--@<option@> @<arg@>] [...]</tt>.
	 * @note Must be overridden by derived classes that have options.
	 */
	virtual std::string ged_valid_options_string_() const;

	/*!
	 * @brief Sets all options to default values.
	 * @note Must be overridden by derived classes that have options.
	 */
	virtual void ged_set_default_options_();

};

}

#endif /* SRC_METHODS_GED_METHOD_HPP_ */
