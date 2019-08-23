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
 * @file all_methods.hpp
 * @brief Includes all methods to enable mutual usability while avoiding circular dependencies.
 * @warning Be careful to respect the dependencies specified in the comments when changing the order of the include directives.
 */

#ifndef SRC_METHODS_ALL_METHODS_HPP_
#define SRC_METHODS_ALL_METHODS_HPP_

// Declarations of abstract classes.

#include "ged_method.hpp"          // Declares GEDMethod.
#ifdef GUROBI
#include "mip_based_method.hpp"    // Declares MIPBasedMethod. Dependencies: "ged_method.hpp".
#endif
#include "lsape_based_method.hpp"  // Declares LSAPEBasedMethod. Dependencies: "ged_method.hpp".
#include "ls_based_method.hpp"     // Declares LSBasedMethod. Dependencies: "lsape_based_method.hpp".
#include "ml_based_method.hpp"     // Declares MLBasedMethod. Dependencies: "lsape_based_method.hpp".

// Declarations of concrete derived classes of GEDMethod.

#include "branch_tight.hpp"        // Declares BranchTight. Dependencies: "ged_method.hpp".
#include "anchor_aware_ged.hpp"    // Declares Exact. Dependencies: "ged_method.hpp".
#include "partition.hpp"           // Declares Partition. Dependencies: "ged_method.hpp".
#include "hybrid.hpp"              // Declares Hybrid. Dependencies: "ged_method.hpp", "partition.hpp".
#include "branch_compact.hpp"      // Declares BranchCompact. Dependencies: "ged_method.hpp".
#include "hed.hpp"                 // Declares HED. Dependencies: "ged_method.hpp".
#include "simulated_annealing.hpp" // Declares SimulatedAnnealing. Dependencies: "ged_method.hpp", "lsape_based_method.hpp"

// Declarations of concrete derived classes of MIPBasedMethod.

#ifdef GUROBI
#include "f1.hpp"				   // Declares F1. Dependencies: "mip_based_method.hpp".
#include "f2.hpp"				   // Declares F2. Dependencies: "mip_based_method.hpp".
#include "compact_mip.hpp"		   // Declares CompactMIP. Dependencies: "mip_based_method.hpp".
#include "blp_no_edge_labels.hpp"  // Declares BLPNoEdgeLabels. Dependencies "mip_based_method.hpp".
#endif

// Declarations of concrete derived classes of LSAPEBasedMethod.

#include "bipartite.hpp"           // Declares Bipartite. Dependencies: "lsape_based_method.hpp".
#include "branch.hpp"              // Declares Branch. Dependencies: "lsape_based_method.hpp".
#include "branch_fast.hpp"         // Declares BranchFast. Dependencies: "lsape_based_method.hpp".
#include "branch_uniform.hpp"      // Declares BranchUniform. Dependencies: "lsape_based_method.hpp".
#include "node.hpp"                // Declares Node. Dependencies: "lsape_based_method.hpp".
#include "ring.hpp"                // Declares Ring. Dependencies: "lsape_based_method.hpp".
#include "star.hpp"                // Declares Star. Dependencies: "lsape_based_method.hpp".
#include "subgraph.hpp"            // Declares Subgraph. Dependencies: "lsape_based_method.hpp".
#include "walks.hpp"               // Declares Walks. Dependencies: "lsape_based_method.hpp".

// Declarations of concrete derived classes of LSBasedMethod.

#include "ipfp.hpp"                // Declares IPFP. Dependencies: "ls_based_method.hpp".
#include "refine.hpp"              // Declares Refine. Dependencies: "ls_based_method.hpp".
#include "bp_beam.hpp"             // Declares BPBeam. Dependencies: "ls_based_method.hpp".

// Declarations of concrete derived classes of MLBasedMethod.

#include "bipartite_ml.hpp"        // Declares BipartiteML. Dependencies: "ml_based_method.hpp".
#include "ring_ml.hpp"             // Declares RingML. Dependencies: "ml_based_method.hpp".

// Definitions of abstract classes.

#include "ged_method.ipp"          // Defines GEDMethod. Dependencies: "ged_method.hpp".
#include "lsape_based_method.ipp"  // Defines LSAPEBasedMethod. Dependencies: "lsape_based_method.hpp", "ged_method.hpp", "branch_tight.hpp".
#include "ls_based_method.ipp"     // Defines LSBasedMethod. Dependencies: "ls_based_method.hpp", "lsape_based_method.hpp", "bipartite.hpp", "branch_hpp", "branch_fast.hpp", "node.hpp", "ring.hpp", "subgraph.hpp", "walks.hpp", "bipartite_ml.hpp", "ring_ml.hpp".
#include "ml_based_method.ipp"     // Defines MLBasedMethod. Dependencies: "ml_based_method.hpp", "lsape_based_method.hpp", "exact.hpp", "ipfp.hpp".
#ifdef GUROBI
#include "mip_based_method.ipp"    // Defines MIPBasedMethod. Dependencies: "mip_based_method.hpp", "ged_method.hpp".
#endif

// Definitions of concrete derived classes of GEDMethod.

#include "branch_tight.ipp"        // Defines BranchTight. Dependencies: "branch_tight.hpp", "ged_method.hpp".
#include "anchor_aware_ged.ipp"    // Defines Exact. Dependencies: "exact.hpp", "ged_method.hpp", "ipfp.hpp".
#include "partition.ipp"           // Defines Partition. Dependencies: "partition.hpp", "ged_method.hpp".
#include "hybrid.ipp"              // Defines Hybrid. Dependencies: "hybrid.hpp", "ged_method.hpp", "partition.hpp", "branch_uniform.hpp".
#include "branch_compact.ipp"      // Defines BranchCompact. Dependencies: "branch_compact.hpp", "ged_method.hpp".
#include "hed.ipp"                 // Defines HED. Dependencies: "hed.hpp", "ged_method.hpp".
#include "simulated_annealing.ipp" // Defines SimulatedAnnealing. Dependencies: "simulated_annealing.hpp", "ged_method.hpp", "lsape_based_method.hpp", "bipartite.hpp", "branch_hpp", "branch_fast.hpp", "node.hpp", "ring.hpp", "subgraph.hpp", "walks.hpp".

// Definitions of concrete derived classes of MIPBasedMethod.

#ifdef GUROBI
#include "f1.ipp"				   // Defines F1. Dependencies: "f1.hpp", "mip_based_method.hpp".
#include "f2.ipp"				   // Defines F2. Dependencies: "f2.hpp", "mip_based_method.hpp".
#include "compact_mip.ipp"		   // Defines CompactMIP. Dependencies: "compact_mip.hpp", "mip_based_method.hpp".
#include "blp_no_edge_labels.ipp"  // Defines BLPNoEdgeLabels. Dependencies "blp_no_edge_labels.hpp", "mip_based_method.hpp".
#endif

// Definitions of concrete derived classes of LSAPEBasedMethod.

#include "bipartite.ipp"           // Defines Bipartite. Dependencies: "bipartite.hpp", "lsape_based_method.hpp".
#include "branch.ipp"              // Defines Branch. Dependencies: "branch.hpp", "lsape_based_method.hpp".
#include "branch_fast.ipp"         // Defines BranchFast. Dependencies: "branch_fast.hpp", "lsape_based_method.hpp".
#include "branch_uniform.ipp"      // Defines BranchUniform. Dependencies: "branch_uniform.hpp", "lsape_based_method.hpp".
#include "node.ipp"                // Defines Node. Dependencies: "node.hpp", "lsape_based_method.hpp".
#include "ring.ipp"                // Defines Ring. Dependencies: "ring.hpp", "lsape_based_method.hpp".
#include "star.ipp"                // Defines Star. Dependencies: "star.hpp", "lsape_based_method.hpp".
#include "subgraph.ipp"            // Defines Subgraph. Dependencies: "subgraph.hpp", "lsape_based_method.hpp", "exact.hpp".
#include "walks.ipp"               // Defines Walks. Dependencies: "walks.hpp", "lsape_based_method.hpp".

// Definitions of concrete derived classes of LSBasedMethod.

#include "ipfp.ipp"                // Defines IPFP. Dependencies: "ipfp.hpp", "ls_based_method.hpp".
#include "refine.ipp"              // Declares Refine. Dependencies: "refine.hpp", "ls_based_method.hpp".
#include "bp_beam.ipp"             // Declares BPBeam. Dependencies: "bp_beam.hpp", "ls_based_method.hpp", "ml_based_method.hpp", "bipartite_ml.ipp", "ring_ml.ipp".

// Definitions of concrete derived classes of MLBasedMethod.

#include "bipartite_ml.ipp"        // Defines BipartiteML. Dependencies: "bipartite_ml.ipp", "ml_based_method.hpp".
#include "ring_ml.ipp"             // Defines RingML. Dependencies: "ring_ml.ipp", "ml_based_method.hpp".


#endif /* SRC_METHODS_ALL_METHODS_HPP_ */
