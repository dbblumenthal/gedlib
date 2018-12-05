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

#include "catch.hpp"

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"
#include<string>

TEST_CASE("testing on MUTA graphs") {
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(ged::Options::EditCosts::CHEM_1);
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs("../../../data/datasets/Mutagenicity/data/", "../collections/MUTA_10.xml"));

	ged::GEDGraph::GraphID g {graph_ids[1]};
	ged::GEDGraph::GraphID h {graph_ids[2]};
	//env.init();
	env.init(ged::Options::InitType::EAGER_WITH_SHUFFLED_COPIES);
	double lower_bound{0.0};
	double exact{0.0};
	double upper_bound;
	double runtime;
	//std::size_t num_runs{1};
	std::size_t num_runs{graph_ids.size() * graph_ids.size()};
	ged::ProgressBar progress(num_runs);
        std::vector<std::string> exp_numsols{"20","10"};
 	std::vector<std::string> exp_numloops{"0","1"};


	SECTION("RANDPOST") {


		for (std::size_t i=0;i<exp_numsols.size();i++){
			std::cout << "\n=== running REFINE RANDPOST (T1, I"<< exp_numsols[i]<<", L"<<exp_numloops[i] <<", R0, P0) ===\n";
			std::cout << "\r" << progress << std::flush;
			env.set_method(ged::Options::GEDMethod::REFINE, "--threads 1 --initial-solutions "+ exp_numsols[i] +" --num-randpost-loops " + exp_numloops[i]);
			upper_bound = 0;
			runtime = 0;
			progress.reset();
			for (ged::GEDGraph::GraphID g : graph_ids) {
				for (ged::GEDGraph::GraphID h : graph_ids) {
					env.run_method(g, h);
					upper_bound += env.get_upper_bound(g, h);
					runtime += env.get_runtime(g, h);
					progress.increment();
					std::cout << "\r" << progress << std::flush;
				}
      }
    }
		std::cout << "\n=== running Star ===\n";
		env.set_method(ged::Options::GEDMethod::STAR, "--threads 5");
		upper_bound = 0;
		lower_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				lower_bound += env.get_lower_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", lower bound = " << lower_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running BranchUniform ===\n";
		env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads 5");
		upper_bound = 0;
		lower_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				lower_bound += env.get_lower_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", lower bound = " << lower_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running BranchFast ===\n";
		env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads 5");
		upper_bound = 0;
		lower_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				lower_bound += env.get_lower_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", lower bound = " << lower_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running Branch ===\n";
		env.set_method(ged::Options::GEDMethod::BRANCH, "--threads 5");
		upper_bound = 0;
		lower_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				lower_bound += env.get_lower_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", lower bound = " << lower_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running BranchTight ===\n";
		env.set_method(ged::Options::GEDMethod::BRANCH_TIGHT, "--threads 5");
		upper_bound = 0;
		lower_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				lower_bound += env.get_lower_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", lower bound = " << lower_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running SimulatedAnnealing ===\n";
		env.set_method(ged::Options::GEDMethod::SIMULATED_ANNEALING, "--threads 5 --lower-bound-method BRANCH_FAST");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running Partition ===\n";
		env.set_method(ged::Options::GEDMethod::PARTITION);
		lower_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				lower_bound += env.get_lower_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nlower bound = " << lower_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running Hybrid ===\n";
		env.set_method(ged::Options::GEDMethod::HYBRID, "--threads 5");
		lower_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				lower_bound += env.get_lower_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nlower bound = " << lower_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running IPFP (QAPE) ===\n";
		env.set_method(ged::Options::GEDMethod::IPFP, "--threads 5 --initial-solutions 4");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";


		std::cout << "\n=== running REFINE ===\n";
		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 5 --initial-solutions 4");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running BP-BEAM ===\n";
		env.set_method(ged::Options::GEDMethod::BP_BEAM, "--threads 5 --initial-solutions 4");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running IBP-BEAM ===\n";
		env.set_method(ged::Options::GEDMethod::BP_BEAM, "--threads 5 --initial-solutions 4 --num-orderings 4");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		/*

		std::cout << "\n=== running REFINE RANDPOST (T5, I5, L0, R0, P0) ===\n";
		std::cout << "\r" << progress << std::flush;
		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 5 --lower-bound-method BRANCH_FAST --initial-solutions 5 --num-randpost-loops 0");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running REFINE RANDPOST (T5, I5, L5, R5, P0) ===\n";
		std::cout << "\r" << progress << std::flush;
		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 5 --lower-bound-method BRANCH_FAST  --initial-solutions 5 --num-randpost-loops 5 --max-randpost-retrials 5");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running REFINE RANDPOST (T5, I5, L5, R5, P1) ===\n";
		std::cout << "\r" << progress << std::flush;
		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 5 --lower-bound-method BRANCH_FAST --initial-solutions 5 --num-randpost-loops 5 --max-randpost-retrials 5 --randpost-penalty 1");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";
		std::cout << "\n=== running IPFP RANDPOST (T4, I5, L0, R0, P0) ===\n";
		std::cout << "\r" << progress << std::flush;
		env.set_method(ged::Options::GEDMethod::IPFP, "--threads 4 --lower-bound-method BRANCH_FAST --initial-solutions 5 --num-randpost-loops 0");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}
		}
		std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";

		std::cout << "\n=== running IPFP RANDPOST (T4, I5, L5, R5, P0) ===\n";
		std::cout << "\r" << progress << std::flush;
		env.set_method(ged::Options::GEDMethod::IPFP, "--threads 4 --lower-bound-method BRANCH_FAST --initial-solutions 5 --num-randpost-loops 5 --max-randpost-retrials 5");
		upper_bound = 0;
		runtime = 0;
		progress.reset();
		for (ged::GEDGraph::GraphID g : graph_ids) {
			for (ged::GEDGraph::GraphID h : graph_ids) {
				env.run_method(g, h);
				upper_bound += env.get_upper_bound(g, h);
				runtime += env.get_runtime(g, h);
				progress.increment();
				std::cout << "\r" << progress << std::flush;
			}

			std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";
		}
		/*
   		for (std::size_t i=0;i<exp_numsols.size();i++){
			std::cout << "\n=== running REFINE RANDPOST (T8, I"<< exp_numsols[i]<<", L"<<exp_numloops[i] <<", R0, P0,SW3) ===\n";
			std::cout << "\r" << progress << std::flush;
			env.set_method(ged::Options::GEDMethod::REFINE, "--threads 8 --initial-solutions "+ exp_numsols[i] +" --num-randpost-loops " + exp_numloops[i] + " --max-swap-size 3");
			upper_bound = 0;
			runtime = 0;
			progress.reset();
			for (ged::GEDGraph::GraphID g : graph_ids) {
				for (ged::GEDGraph::GraphID h : graph_ids) {
					env.run_method(g, h);
					upper_bound += env.get_upper_bound(g, h);
					runtime += env.get_runtime(g, h);
					progress.increment();
					std::cout << "\r" << progress << std::flush;
				}
			}

			std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";
		}
		*/

		for (std::size_t i=0;i<exp_numsols.size();i++){
			std::cout << "\n=== running IPFP RANDPOST (T1, I"<< exp_numsols[i]<<", L"<<exp_numloops[i] <<", R0, P0) ===\n";
			std::cout << "\r" << progress << std::flush;
			env.set_method(ged::Options::GEDMethod::IPFP, " --threads 1 --initial-solutions "+ exp_numsols[i] +" --num-randpost-loops " + exp_numloops[i]);
			upper_bound = 0;
			runtime = 0;
			progress.reset();
			for (ged::GEDGraph::GraphID g : graph_ids) {
				for (ged::GEDGraph::GraphID h : graph_ids) {
					env.run_method(g, h);
					upper_bound += env.get_upper_bound(g, h);
					runtime += env.get_runtime(g, h);
					progress.increment();
					std::cout << "\r" << progress << std::flush;
				}
			}

			std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";
		}


		 */
	}

	/*
	SECTION("Branch, BranchFast, BranchTight, Node") {
		std::cout << "\n=== running Branch ===\n";
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model FLWC --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --threads 5\"\n";
		lower_bound = env.get_lower_bound(g, h);
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound));
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model FLWC --centrality-method PAGERANK --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --centrality-method PAGERANK --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model FLWC --centrality-method EIGENVECTOR --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --centrality-method EIGENVECTOR --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model ECBP --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model ECBP --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		lower_bound = env.get_lower_bound(g, h);
		CHECK(env.get_lower_bound(g, h) <= env.get_upper_bound(g, h));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model EBP --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model EBP --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model FBP --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FBP --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model SFBP --threads 5 --max-num-solutions 10");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model SFBP --threads 5 --max-num-solutions 10\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH, "--lsape-model FBP0 --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FBP0 --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));

		std::cout << "\n=== running BranchFast ===\n";
		double lower_bound_branch_fast;
		env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--lsape-model FLWC --threads 1");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --threads 1\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		lower_bound_branch_fast = env.get_lower_bound(g, h);
		CHECK(env.get_lower_bound(g, h) <= Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--lsape-model EBP --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model EBP --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound_branch_fast));
		CHECK(env.get_lower_bound(g, h) <= Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--lsape-model FLWC --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound_branch_fast));
		CHECK(env.get_lower_bound(g, h) <= Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--lsape-model FBP --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FBP --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound_branch_fast));
		CHECK(env.get_lower_bound(g, h) <= Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--lsape-model SFBP --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model SFBP --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound_branch_fast));
		CHECK(env.get_lower_bound(g, h) <= Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--lsape-model FBP0 --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FBP0 --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) == Approx(lower_bound_branch_fast));
		CHECK(env.get_lower_bound(g, h) <= Approx(lower_bound));
		CHECK(exact <= env.get_upper_bound(g, h));

		std::cout << "\n=== running BranchTight ===\n";
		env.set_method(ged::Options::GEDMethod::BRANCH_TIGHT, "--upper-bound BEST --iterations 20 --threads 5");
		env.run_method(g, h);
		std::cout << "\noptions = \"--upper-bound BEST --iterations 20 --threads 5\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) >= Approx(lower_bound));
		CHECK(env.get_lower_bound(g, h) <= env.get_upper_bound(g, h));

		std::cout << "\n=== running NODE ===\n";
		env.set_method(ged::Options::GEDMethod::NODE);
		env.run_method(g, h);
		std::cout << "\noptions = \"\"\n";
		std::cout << "lower_bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(env.get_lower_bound(g, h) <= Approx(lower_bound));
		CHECK(env.get_lower_bound(g, h) <= env.get_upper_bound(g, h));

		std::cout << "\n=== running BIPARTITE ===\n";
		env.set_method(ged::Options::GEDMethod::BIPARTITE);
		env.run_method(g, h);
		std::cout << "\noptions = \"\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}

	SECTION("HED") {
		std::cout << "\n===running HED ===\n";
		env.set_method(ged::Options::GEDMethod::HED, "--threads 1");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 1\"\n";
		std::cout << "lower bound = " << env.get_lower_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		env.set_method(ged::Options::GEDMethod::HED, "--threads 4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4\"\n";
		std::cout << "lower bound = " << env.get_lower_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}

#ifdef GUROBI
	SECTION("BLPNoEdgeLabels") {
		std::cout << "\n===running BLPNoEdgeLabels ===\n";
		env.set_method(ged::Options::GEDMethod::BLP_NO_EDGE_LABELS, "--threads 4 --tune TRUE --relax TRUE");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --relax TRUE\"\n";
		std::cout << "lower bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::BLP_NO_EDGE_LABELS, "--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1\"\n";
		std::cout << "lower bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}
	SECTION("F2") {
		std::cout << "\n===running F2 ===\n";
		env.set_method(ged::Options::GEDMethod::F2, "--threads 4 --tune TRUE --relax TRUE");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --relax TRUE\"\n";
		std::cout << "lower bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::F2, "--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}
	SECTION("F1") {
		std::cout << "\n===running F1 ===\n";
		env.set_method(ged::Options::GEDMethod::F1, "--threads 4 --tune TRUE --relax TRUE");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --relax TRUE\"\n";
		std::cout << "lower bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::F1, "--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}
	SECTION("CompactMIP") {
		std::cout << "\n===running CompactMIP ===\n";
		env.set_method(ged::Options::GEDMethod::COMPACT_MIP, "--threads 4 --tune TRUE --relax TRUE");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --relax TRUE\"\n";
		std::cout << "lower bound = " << env.get_lower_bound(g, h) << ", upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::COMPACT_MIP, "--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --tune TRUE --time-limit 10 --tune-time-limit 1\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}
#endif

	SECTION("Refine") {
		std::cout << "\n===running Refine ===\n";
		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 4 --initial-solutions 8 --runs-from-initial-solutions 4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 4 --initial-solutions 8 --runs-from-initial-solutions 4\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 1 --initial-solutions 4 --initialization-method BRANCH_FAST");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 1 --initial-solutions 4 --initialization-method BRANCH_FAST\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 1 --initialization-method BRANCH --initial-solutions 10 --lower-bound-method BRANCH_TIGHT");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 1 --initialization-method BRANCH --initial-solutions 10 --lower-bound-method BRANCH_TIGHT\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}

	SECTION("Simulated Annealing") {
		std::cout << "\n===running Simulated Annealing ===\n";
		env.set_method(ged::Options::GEDMethod::SIMULATED_ANNEALING, "--lsape-method BIPARTITE");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-method BIPARTITE\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::SIMULATED_ANNEALING, "--lsape-method NODE");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-method NODE\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}


	SECTION("AnchorAwareGED") {
		std::cout << "\n===running AnchorAwareGED ===\n";
		env.set_method(ged::Options::GEDMethod::ANCHOR_AWARE_GED, "--threads 5 --time-limit 1.0");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --time-limit 1.0\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
	}

	SECTION("Subgraph") {
		std::cout << "\n=== Subgraph ===\n";
		env.set_method(ged::Options::GEDMethod::SUBGRAPH, "--threads 5 --subproblem-solver-options '--time-limit 0.0001'");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --subproblem-solver-options '--time-limit 0.0001'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
#ifdef GUROBI
		env.set_method(ged::Options::GEDMethod::SUBGRAPH, "--threads 5 --subproblem-solver F2 --subproblem-solver-options '--time-limit 0.01 --relax TRUE'");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --subproblem-solver F2 --subproblem-solver-options '--time-limit 0.01 --relax TRUE'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
#endif

	}

	SECTION("Walks") {
		std::cout << "\n===running Walks ===\n";
		env.set_method(ged::Options::GEDMethod::WALKS, "--lsape-model FLWC --depth-range 4,4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --depth-range 4,4\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::WALKS, "--lsape-model FLWC --depth-range 3,3");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --depth-range 3,3\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::WALKS, "--lsape-model FLWC --depth-range 2,2");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --depth-range 2,2\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::WALKS, "--lsape-model FLWC --depth-range 1,1");
		env.run_method(g, h);
		std::cout << "\noptions = \"--lsape-model FLWC --depth-range 1,1\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
	}

	SECTION("IPFP") {
		std::cout << "\n===running IPFP ===\n";
		env.set_method(ged::Options::GEDMethod::IPFP, "--quadratic-model B-QAP --threads 5 --initial-solutions 4 --initialization-method BRANCH_FAST --initialization-options '--threads 5'");
		env.run_method(g, h);
		std::cout << "\noptions = \"--quadratic-model B-QAP --threads 5 --initial-solutions 4 --initialization-method BRANCH_FAST --initialization-options '--threads 5'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::IPFP, "--quadratic-model B-QAP --threads 5 --initialization-method BRANCH --initial-solutions 10 --lower-bound-method BRANCH_TIGHT");
		env.run_method(g, h);
		std::cout << "\noptions = \"--quadratic-model B-QAP --threads 5 --initialization-method BRANCH --initial-solutions 10 --lower-bound-method BRANCH_TIGHT\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::IPFP, "--quadratic-model QAPE --threads 5 --initial-solutions 4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--quadratic-model QAPE --threads 5 --initial-solutions 4\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::IPFP, "--quadratic-model QAPE --threads 5 --initialization-method BRANCH --initial-solutions 4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--quadratic-model QAPE --threads 5 --initialization-method BRANCH --initial-solutions 4\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::IPFP, "--quadratic-model C-QAP --threads 5 --initial-solutions 4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--quadratic-model C-QAP --threads 5 --initial-solutions 4\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));

		env.set_method(ged::Options::GEDMethod::IPFP, "--quadratic-model C-QAP --threads 5 --initialization-method BRANCH --initial-solutions 4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--quadratic-model C-QAP --threads 5 --initialization-method BRANCH --initial-solutions 4\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
	}

	SECTION("Ring") {
		std::cout << "\n===running Ring ===\n";
		env.set_method(ged::Options::GEDMethod::RING, "--threads 5 --led-method LSAPE_OPTIMAL --save ../output/mao_ring_LSAPE_OPTIMAL.ini");
		env.init_method();
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --led-method LSAPE_OPTIMAL --save ../output/MAO/ring_LSAPE_OPTIMAL.ini\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::RING, "--threads 5 --led-method LSAPE_GREEDY --save ../output/mao_ring_LSAPE_GREEDY.ini");
		env.init_method();
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --led-method LSAPE_GREEDY --save ../output/mao_ring_LSAPE_GREEDY.ini\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::RING, "--threads 5 --led-method GAMMA --save ../output/mao_ring_GAMMA.ini");
		env.init_method();
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --led-method GAMMA --save ../output/mao_ring_GAMMA.ini\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
	}

	SECTION("BipartiteML") {
		std::cout << "\n=== BipartiteML ===\n";
		env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, "--ml-method DNN --save ../output/mao_bipartite_ml_dnn.ini --save-train ../output/mao_bipartite_ml.data --ground-truth-options '--initial-solutions 1'");
		env.init_method();
		env.run_method(g, h);
		std::cout << "\noptions = \"--ml-method DNN --save ../output/mao_bipartite_ml_dnn.ini --save-train ../output/mao_bipartite.train --ground-truth-options '--initial-solutions 1'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, "--ml-method SVM --save ../output/mao_bipartite_ml_svm.ini --load-train ../output/mao_bipartite_ml.data --ground-truth-options '--initial-solutions 1'");
		env.init_method();
		env.run_method(g, h);
		std::cout << "\noptions = \"--ml-method SVM --save ../output/mao_bipartite_ml_svm.ini --load-train ../output/mao_bipartite.train --ground-truth-options '--initial-solutions 1'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
	}

	SECTION("RingML") {
		std::cout << "\n=== RingML ===\n";
		env.set_method(ged::Options::GEDMethod::RING_ML, "--ml-method DNN --save ../output/mao_ring_ml_dnn.ini --save-train ../output/mao_ring_ml.data --ground-truth-options '--initial-solutions 1'");
		env.init_method();
		env.run_method(g, h);
		std::cout << "\noptions = \"--ml-method DNN --save ../output/mao_ring_ml_dnn.ini --save-train ../output/mao_ring_ml.data --ground-truth-options '--initial-solutions 1'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::RING_ML, "--ml-method SVM --save ../output/mao_ring_ml_svm.ini --load-train ../output/mao_ring_ml.data --ground-truth-options '--initial-solutions 1'");
		env.init_method();
		env.run_method(g, h);
		std::cout << "\noptions = \"--ml-method SVM --save ../output/mao_ring_ml_svm.ini --load-train ../output/mao_ring_ml.data --ground-truth-options '--initial-solutions 1'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
	}

	SECTION("BPBeam") {
		std::cout << "\n===running BPBeam ===\n";
		env.set_method(ged::Options::GEDMethod::BP_BEAM, "--initial-solutions 4");
		env.run_method(g, h);
		std::cout << "\noptions = \"--initial-solutions 4\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::BP_BEAM, "--initial-solutions 4 --ordering-method RANDOM");
		env.run_method(g, h);
		std::cout << "\noptions = \"--initial-solutions 4 --ordering-method RANDOM\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::BP_BEAM, "--initial-solutions 4 --ordering-method RING_ML --ordering-options '--ml-method DNN --load ../output/MAO/ring_ml_dnn.ini'");
		env.run_method(g, h);
		std::cout << "\noptions = \"--initial-solutions 4 --ordering-method RING_ML --ordering-options '--ml-method DNN --load ../output/mao_ring_ml_dnn.ini'\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}*/

}




