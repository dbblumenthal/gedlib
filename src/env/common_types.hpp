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
#include <initializer_list>
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
#include <lsap.h>
#include <lsape.h>
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
#endif
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
		LAZY_WITHOUT_SHUFFLED_COPIES, //!< Lazy initialization, no shuffled graph copies are constructed.
		EAGER_WITHOUT_SHUFFLED_COPIES, //!< Eager initialization, no shuffled graph copies are constructed.
		LAZY_WITH_SHUFFLED_COPIES, //!< Lazy initialization, shuffled graph copies are constructed.
		EAGER_WITH_SHUFFLED_COPIES //!< Eager initialization, shuffled graph copies are constructed.
	};

};

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

}

#endif /* SRC_ENV_COMMON_TYPES_HPP_ */
