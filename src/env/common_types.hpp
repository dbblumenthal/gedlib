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
 * @file common_types.hpp
 * @brief Type declarations used by various classes.
 */

#ifndef SRC_ENV_COMMON_TYPES_HPP_
#define SRC_ENV_COMMON_TYPES_HPP_

// Include standard libraries.
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <string>
#include <unordered_set>
#include <chrono>
#include <cmath>
#include <iostream>
#include <ios>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <sstream>
#include <typeinfo>
#ifdef _OPENMP
#include <omp.h>
#endif

// Include Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

// Include Gurobi.
#ifdef GUROBI
#include <gurobi_c++.h>
#endif

// Include external non-standard libraries.
#ifndef LSAPE_IndexType
#define LSAPE_IndexType std::size_t
#endif
#include <lsap.hh>
#include <lsape.hh>
#include <svm.h>
#include <doublefann.h>
#include <fann_cpp.h>
#include <Dense>
#include <nomad.hpp>

// Include helper classes.
#include "timer.hpp"
#include "error.hpp"
#include "progress_bar.hpp"


namespace ged {

/*!
 * @brief Internally used type for measurements in seconds.
 */
typedef std::chrono::duration<double> Seconds;

/*!
 * @brief Type of node and edge labels of graphs given in the .gxl file format.
 */
typedef std::map<std::string, std::string> GXLLabel;

/*!
 * @brief Streams std::map.
 * @param[in] os Output stream
 * @param[in] map Map that should be streamed.
 * @return Output stream.
 */
template<class Key, class Value>
std::ostream & operator<<(std::ostream & os, const std::map<Key, Value> & map) {
	os << "{ ";
	for (const auto & key_val : map) {
		os << "(" << key_val.first << "," << key_val.second << ") ";
	}
	os << "}";
	return os;
}

/*!
 * @brief Type of node IDs of graphs given in the .gxl file format.
 */
typedef std::string GXLNodeID;

/*!
 * @brief Internally used type of node and edge labels.
 */
typedef std::size_t LabelID;

/*!
 * @brief Type of dummy labels for unlabeled nodes and edges.
 */
struct NoLabel {};

constexpr bool operator==(NoLabel const &, NoLabel const &) {return true;}

/*!
 * @brief Returns an invalid label.
 * @return Invalid label.
 */
constexpr LabelID invalid_label() {return std::numeric_limits<LabelID>::max();}

/*!
 * @brief Returns a dummy label.
 * @return Dummy label.
 */
constexpr LabelID dummy_label() {return 0;}


/*!
 * @brief Returns undefined size.
 * @return Undefined size.
 */
constexpr std::size_t undefined() {return std::numeric_limits<std::size_t>::max();}

/*!
 * @brief Provides constant @f$\pi@f$ with double precision.
 * @return Returns 3.141592653589793238463.
 */
constexpr double pi() {return 3.141592653589793238463;}

/*!
 * @brief Contains enums for options employed by ged::GEDEnv.
 */
struct Options {


	/*!
	 * @brief Selects the method.
	 */
	enum class GEDMethod {
#ifdef GUROBI
		F1,                  //!< Selects ged::F1.
		F2,                  //!< Selects ged::F2.
		COMPACT_MIP,         //!< Selects ged::CompactMIP.
		BLP_NO_EDGE_LABELS,  //!< Selects ged::BLPNoEdgeLabels.
#endif /* GUROBI */
		BRANCH,              //!< Selects ged::Branch.
		BRANCH_FAST,         //!< Selects ged::BranchFast.
		BRANCH_TIGHT,        //!< Selects ged::BranchTight.
		BRANCH_UNIFORM,      //!< Selects ged::BranchUniform.
		BRANCH_COMPACT,      //!< Selects ged::BranchCompact.
		PARTITION,           //!< Selects ged::Partition.
		HYBRID,              //!< Selects ged::Hybrid.
		RING,                //!< Selects ged::Ring.
		ANCHOR_AWARE_GED,    //!< Selects ged::AnchorAwareGED.
		WALKS,               //!< Selects ged::Walks.
		IPFP,                //!< Selects ged::IPFP
		BIPARTITE,           //!< Selects ged::Bipartite.
		SUBGRAPH,            //!< Selects ged::Subgraph.
		NODE,                //!< Selects ged::Node.
		RING_ML,             //!< Selects ged::RingML.
		BIPARTITE_ML,        //!< Selects ged::BipartiteML.
		REFINE,              //!< Selects ged::Refine.
		BP_BEAM,             //!< Selects ged::BPBeam.
		SIMULATED_ANNEALING, //!< Selects ged::SimulatedAnnealing.
		HED,				 //!< Selects ged::HED.
		STAR				 //!< Selects ged::Star.
	};

	/*!
	 * @brief Selects the edit costs.
	 */
	enum class EditCosts {
		CHEM_1,      //!< Selects ged::CHEM1.
		CHEM_2,      //!< Selects ged::CHEM2.
		CMU,         //!< Selects ged::CMU.
		GREC_1,      //!< Selects ged::GREC1.
		GREC_2,      //!< Selects ged::GREC2.
		PROTEIN,     //!< Selects ged::Protein.
		FINGERPRINT, //!< Selects ged::Fingerprint.
		LETTER,      //!< Selects ged::Letter.
		CONSTANT     //!< Selects ged::Constant.
	};

	/*!
	 * @brief Selects whether nodes or edges of graphs given in GXL file format are labeled or unlabeled.
	 */
	enum class GXLNodeEdgeType {
		LABELED,  //!< Labeled nodes or edges.
		UNLABELED //!< Unlabeled nodes or edges.
	};

	/*!
	 * @brief Selects the initialization type of the environment.
	 * @details If eager initialization is selected, all edit costs are pre-computed when initializing the environment.
	 * Otherwise, they are computed at runtime. If initialization with shuffled copies is selected, shuffled copies of
	 * all graphs are created. These copies are used when calling ged::GEDEnv::run_method() with two identical graph IDs.
	 * In this case, one of the IDs is internally replaced by the ID of the shuffled copy and the graph is hence
	 * compared to an isomorphic but non-identical graph. If initialization without shuffled copies is selected, no shuffled copies
	 * are created and calling ged::GEDEnv::run_method() with two identical graph IDs amounts to comparing a graph to itself.
	 */
	enum class InitType {
		LAZY_WITHOUT_SHUFFLED_COPIES,  //!< Lazy initialization, no shuffled graph copies are constructed.
		EAGER_WITHOUT_SHUFFLED_COPIES, //!< Eager initialization, no shuffled graph copies are constructed.
		LAZY_WITH_SHUFFLED_COPIES,     //!< Lazy initialization, shuffled graph copies are constructed.
		EAGER_WITH_SHUFFLED_COPIES     //!< Eager initialization, shuffled graph copies are constructed.
	};

	/*!
	 * @brief Specifies type of exchange graph.
	 */
	enum class ExchangeGraphType {
		ADJ_MATRIX,//!< Exchange graph is given as adjacency matrix.
		ADJ_LISTS, //!< Exchange graph is given as adjacency lists.
		EDGE_LIST  //!< Exchange graph is given as list of edges.
	};

	/*!
	 * @brief can be used to specify the state of an algorithm.
	 */
	enum class AlgorithmState {
		CALLED,     //!< The algorithm has been called.
		INITIALIZED,//!< The algorithm has been initialized.
		CONVERGED,  //!< The algorithm has converged.
		TERMINATED  //!< The algorithm has terminated.
	};

};

std::ostream & operator<<(std::ostream & os, const Options::AlgorithmState & state) {
	switch (state) {
	case Options::AlgorithmState::CALLED:
		os << 0;
		break;
	case Options::AlgorithmState::INITIALIZED:
		os << 1;
		break;
	case Options::AlgorithmState::CONVERGED:
		os << 2;
		break;
	case Options::AlgorithmState::TERMINATED:
		os << 3;
		break;
	}
	return os;

}

/*!
 * @brief Streams Options::GEDMethod object.
 * @param[in] os Output stream.
 * @param[in] ged_method Method selector.
 * @return Output stream.
 */
std::ostream & operator<<(std::ostream & os, const Options::GEDMethod & ged_method) {
	switch (ged_method) {
	case Options::GEDMethod::BRANCH:
		os << "BRANCH";
		break;
	case Options::GEDMethod::BRANCH_FAST:
		os << "BRANCH_FAST";
		break;
	case Options::GEDMethod::BRANCH_TIGHT:
		os << "BRANCH_TIGHT";
		break;
	case Options::GEDMethod::BRANCH_UNIFORM:
		os << "BRANCH_UNIFORM";
		break;
	case Options::GEDMethod::BRANCH_COMPACT:
		os << "BRANCH_COMPACT";
		break;
	case Options::GEDMethod::PARTITION:
		os << "PARTITION";
		break;
	case Options::GEDMethod::HYBRID:
		os << "HYBRID";
		break;
	case Options::GEDMethod::RING:
		os << "RING";
		break;
	case Options::GEDMethod::ANCHOR_AWARE_GED:
		os << "ANCHOR_AWARE_GED";
		break;
	case Options::GEDMethod::WALKS:
		os << "WALKS";
		break;
	case Options::GEDMethod::IPFP:
		os << "IPFP";
		break;
	case Options::GEDMethod::BIPARTITE:
		os << "BIPARTITE";
		break;
	case Options::GEDMethod::SUBGRAPH:
		os << "SUBGRAPH";
		break;
	case Options::GEDMethod::NODE:
		os << "NODE";
		break;
	case Options::GEDMethod::RING_ML:
		os << "RING_ML";
		break;
	case Options::GEDMethod::BIPARTITE_ML:
		os << "BIPARTITE_ML";
		break;
	case Options::GEDMethod::REFINE:
		os << "REFINE";
		break;
	case Options::GEDMethod::BP_BEAM:
		os << "BP_BEAM";
		break;
	case Options::GEDMethod::SIMULATED_ANNEALING:
		os << "SIMULATED_ANNEALING";
		break;
	case Options::GEDMethod::HED:
		os << "HED";
		break;
	case Options::GEDMethod::STAR:
		os << "STAR";
		break;
#ifdef GUROBI
	case Options::GEDMethod::F1:
		os << "F1";
		break;
	case Options::GEDMethod::F2:
		os << "F2";
		break;
	case Options::GEDMethod::COMPACT_MIP:
		os << "COMPACT_MIP";
		break;
	case Options::GEDMethod::BLP_NO_EDGE_LABELS:
		os << "BLP_NO_EDGE_LABELS";
		break;
#endif
	}
	return os;
}

/*!
 * @brief Simple graph class used for communication with user.
 * @tparam UserNodeID Class of user-specific node IDs.
 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
 */
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel> struct ExchangeGraph {

	std::size_t id;                                                                     //!< Internal ID of the graph.

	std::size_t num_nodes;                                                              //!< The number of nodes. Nodes have IDs between @p 0 and <tt>num_nodes - 1</tt>.

	std::size_t num_edges;                                                              //!< The number of edges.

	std::vector<UserNodeID> original_node_ids;                                          //!< The original IDs of all nodes.

	std::vector<UserNodeLabel> node_labels;                                             //!< The labels of all nodes.

	std::vector<std::vector<std::size_t>> adj_matrix;                                   //!< Adjacency matrix.

	std::map<std::pair<std::size_t, std::size_t>, UserEdgeLabel> edge_labels;           //!< A hash map with a key-value pair <tt>((internal_tail_id, internal_head_id), label)</tt> for each edge.

	std::vector<std::list<std::pair<std::size_t, UserEdgeLabel>>> adj_lists;            //!< Adjacency lists for all nodes.

	std::list<std::pair<std::pair<std::size_t, std::size_t>, UserEdgeLabel>> edge_list; //!< A list of all edges.

	bool operator==(const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & rhs) const {
		return ((original_node_ids == rhs.original_node_ids) and (node_labels == rhs.node_labels) and (adj_matrix == rhs.adj_matrix) and (edge_labels == rhs.edge_labels) and (adj_lists == rhs.adj_lists) and (edge_list == rhs.edge_list));
	};
};

#ifdef ENABLE_GRAPH_STREAMING
/*!
 * @brief Streams ged::ExchangeGraph object in GML format.
 * @tparam UserNodeID Class of user-specific node IDs.
 * @tparam UserNodeLabel Class of user-specific node labels. If nodes are unlabeled, use ged::NoLabel or define your own dummy label class.
 * @tparam UserEdgeLabel Class of user-specific edge labels. If edges are unlabeled, use ged::NoLabel or define your own dummy label class.
 * @param os Output stream.
 * @param[in] graph The graph that should be streamed.
 * @return Output stream.
 * @note Define @p ENABLE_GRAPH_STREAMING to compile this function. Requires the <tt>operator\<\<<\tt> to be implemented for the classes UserNodeID, UserNodeLabel, and UserEdgeLabel.
 */
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::ostream & operator<<(std::ostream & os, const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & graph) {
	os << "graph [\n";
	os << "\tid = " << graph.id << "\n";
	os << "\tdirected = 0\n";
	os << "\tnum_nodes = " << graph.num_nodes << "\n";
	os << "\tnum_edges = " << graph.num_edges << "\n";
	for (std::size_t node_id{0}; node_id < graph.num_nodes; node_id++) {
		os << "\tnode [\n";
		os << "\t\tinternal_id = " << node_id << "\n";
		os << "\t\toriginal_id = \"" << graph.original_node_ids.at(node_id) << "\"\n";
		os << "\t\tlabel = \"" << graph.node_labels.at(node_id) << "\"\n";
		os << "\t]\n";
	}
	for (std::size_t tail_id{0}; tail_id < graph.num_nodes; tail_id++) {
		for (std::size_t head_id{tail_id + 1}; head_id < graph.num_nodes; head_id++) {
			if (graph.adj_matrix[tail_id][head_id] == 1) {
				os << "\tedge [\n";
				os << "\t\tinternal_id_tail = " << tail_id << "\n";
				os << "\t\tinternal_id_head = " << head_id << "\n";
				os << "\t\tlabel = \"" << graph.edge_labels[std::make_pair(tail_id, head_id)] << "\"\n";
				os << "\t]\n";
			}
		}
	}
	os << "]\n";
	return os;
}
#endif /* ENABLE_GRAPH_STREAMING */

}

#endif /* SRC_ENV_COMMON_TYPES_HPP_ */
