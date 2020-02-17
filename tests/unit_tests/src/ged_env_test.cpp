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

TEST_CASE("determining number of AIDS node labels") {
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::unordered_set<std::string> irrelevant_attributes{"x", "y", "symbol", "charge"};
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs("../../../data/datasets/AIDS/data/", "../collections/AIDS.xml", ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_attributes));
	env.set_edit_costs(ged::Options::EditCosts::CHEM_2);
	std::size_t num_node_labels{env.num_node_labels()};

	std::cout << "number of AIDS node labels: " << num_node_labels << "\n[";
	for (std::size_t label_id{1}; label_id <= num_node_labels; label_id++) {
		std::cout << env.get_node_label(label_id)["chem"] << ",";
	}
	std::cout << "]\n";
}

TEST_CASE("testing on Letter graphs") {
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs("../../../data/datasets/Mutagenicity/data/", "../collections/MUTA_10.xml"));
	env.set_edit_costs(ged::Options::EditCosts::CHEM_2);

	ged::GEDGraph::GraphID g {graph_ids[1]};
	//env.init();
	env.init(ged::Options::InitType::EAGER_WITH_SHUFFLED_COPIES);
	std::vector<ged::GEDGraph::GraphID> graph_ids_2(env.load_gxl_graphs("../../../data/datasets/Mutagenicity/data/", "../collections/MUTA_20.xml"));
	env.init(ged::Options::InitType::EAGER_WITH_SHUFFLED_COPIES);
	ged::GEDGraph::GraphID h {graph_ids[0]};
	//std::cout << env.get_graph(g);
	double lower_bound{0.0};
	double exact{0.0};
	double upper_bound;
	double runtime;
	//std::size_t num_runs{1};
	std::size_t num_runs{graph_ids.size() * graph_ids.size()};
	ged::ProgressBar progress(num_runs);
	std::vector<std::string> exp_numsols{"20","10"};
	std::vector<std::string> exp_numloops{"0","1"};

	SECTION("HED") {
		std::cout << "\n=== running HED ===\n";
		env.set_method(ged::Options::GEDMethod::HED);
		env.run_method(g, h);
		std::cout << "\nupper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
	}


	SECTION("RANDPOST") {

		std::cout << "\n=== running REFINE 40, 0.5, 1 ===\n";
		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 1 --initial-solutions 40 --randomness PSEUDO");
		env.run_method(g, h);
		std::cout << "\nupper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";

		env.set_method(ged::Options::GEDMethod::REFINE, "--threads 1 --initial-solutions 40 --randomness PSEUDO --naive TRUE");
		env.run_method(g, h);
		std::cout << "\nupper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";


		for (std::size_t i=0;i<exp_numsols.size();i++){
			std::cout << "\n=== running REFINE RANDPOST (T1, I"<< exp_numsols[i]<<", L"<<exp_numloops[i] <<", R0, P0) ===\n";
			std::cout << "\r" << progress << std::flush;
			env.set_method(ged::Options::GEDMethod::REFINE, "--threads 1 --initial-solutions 40 --randomness PSEUDO");
			upper_bound = 0;
			runtime = 0;
			progress.reset();
			for (ged::GEDGraph::GraphID g : graph_ids_2) {
				for (ged::GEDGraph::GraphID h : graph_ids_2) {
					env.run_method(g, h);
					upper_bound += env.get_upper_bound(g, h);
					runtime += env.get_runtime(g, h);
					progress.increment();
					std::cout << "\r" << progress << std::flush;
				}
			}
			std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";
		}

		for (std::size_t i=0;i<exp_numsols.size();i++){
			std::cout << "\n=== running REFINE RANDPOST (T8, I"<< exp_numsols[i]<<", L"<<exp_numloops[i] <<", R0, P0,SW3) ===\n";
			std::cout << "\r" << progress << std::flush;
			env.set_method(ged::Options::GEDMethod::REFINE, "--threads 8 --initial-solutions "+ exp_numsols[i] +" --num-randpost-loops " + exp_numloops[i] + " --max-swap-size 3");
			upper_bound = 0;
			runtime = 0;
			progress.reset();
			for (ged::GEDGraph::GraphID g : graph_ids) {
				for (ged::GEDGraph::GraphID h : graph_ids) {
					std::cout << "\r" << progress << " ID " << g << ": " << env.get_graph_name(g) << " ID " << h << ": "<< env.get_graph_name(h) << std::flush;
					env.run_method(g, h);
					upper_bound += env.get_upper_bound(g, h);
					runtime += env.get_runtime(g, h);
					progress.increment();
				}
			}

			std::cout << "\nupper bound = " << upper_bound / static_cast<double>(num_runs) << ", runtime = " << runtime / static_cast<double>(num_runs) << "\n";
		}


		for (std::size_t i=0;i<exp_numsols.size();i++){
			std::cout << "\n=== running IPFP RANDPOST (T5, I"<< exp_numsols[i]<<", L"<<exp_numloops[i] <<", R0, P0) ===\n";
			std::cout << "\r" << progress << std::flush;
			env.set_method(ged::Options::GEDMethod::IPFP, " --threads 5 --initial-solutions "+ exp_numsols[i] +" --num-randpost-loops " + exp_numloops[i]);
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



	}


}




