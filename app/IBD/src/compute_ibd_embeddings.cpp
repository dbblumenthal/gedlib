/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2020 by David B. Blumenthal                              *
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
 * @file  compute_ibd_embeddings.cpp
 * @brief Computes GED embeddings of IBD graphs.
 */

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"
#include "ibd_costs.hpp"
#include "CLI11.hpp"

int main(int argc, char* argv[]) {

	// Parse command line argument.
	CLI::App app;
	ged::Options::GEDMethod method{ged::Options::GEDMethod::BRANCH};
	std::map<std::string, ged::Options::GEDMethod> map{{"BRANCH", ged::Options::GEDMethod::BRANCH}, {"BRANCH_FAST", ged::Options::GEDMethod::BRANCH_FAST}};
	app.add_option("-m,--method", method, "Employed GED method.")->transform(CLI::CheckedTransformer(map, CLI::ignore_case));
	CLI11_PARSE(app, argc, argv);

	// Load the graphs and initialize the environment.
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::string graph_collection("../data/IBD.xml");
	std::string graph_dir("../data/");
	std::unordered_set<std::string> irrelevant_node_attributes;
	irrelevant_node_attributes.insert("Bacteria");
	IBDCosts<ged::GXLLabel, ged::GXLLabel> ibd_costs("../data/phylogenetic_OTU_distances.csv");
	env.load_gxl_graphs(graph_dir, graph_collection, ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes);
	env.set_edit_costs(&ibd_costs);
	env.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);

	// Set and initialize the GED method.
	env.set_method(method, "--threads 10");
	env.init_method();

	// Compute the lower and upper bounds for GED.
	ged::DMatrix lower_bounds(env.num_graphs(), env.num_graphs(), 0.0);
	ged::DMatrix upper_bounds(env.num_graphs(), env.num_graphs(), 0.0);
	ged::ProgressBar progress((env.num_graphs() * (env.num_graphs() - 1)) / 2);
	std::cout << "\rComputing patient embeddings: " << progress << std::flush;
	for (std::size_t g_id{0}; g_id < env.num_graphs(); g_id++) {
		for (std::size_t h_id{g_id + 1}; h_id < env.num_graphs(); h_id++) {
			env.run_method(g_id, h_id);
			lower_bounds(g_id, h_id) = env.get_lower_bound(g_id, h_id);
			lower_bounds(h_id, g_id) = env.get_lower_bound(g_id, h_id);
			upper_bounds(g_id, h_id) = env.get_upper_bound(g_id, h_id);
			upper_bounds(h_id, g_id) = env.get_upper_bound(g_id, h_id);
			progress.increment();
			std::cout << "\rComputing patient embeddings: " << progress << std::flush;
		}
	}
	std::cout << "\n";

	// Write the embeddings as CSV files.
	std::ofstream lb_embedding("../output/LB_embedding.csv");
	std::ofstream ub_embedding("../output/UB_embedding.csv");
	std::ofstream mean_embedding("../output/MEAN_embedding.csv");
	for (std::size_t g_id{0}; g_id < env.num_graphs(); g_id++) {
		lb_embedding << env.get_graph_name(g_id);
		ub_embedding << env.get_graph_name(g_id);
		mean_embedding << env.get_graph_name(g_id);
		for (std::size_t h_id{0}; h_id < env.num_graphs(); h_id++) {
			lb_embedding << "," << lower_bounds(g_id, h_id);
			ub_embedding << "," << upper_bounds(g_id, h_id);
			mean_embedding << "," << (lower_bounds(g_id, h_id) + upper_bounds(g_id, h_id)) / 2.0;
		}
		lb_embedding << "\n";
		ub_embedding << "\n";
		mean_embedding << "\n";
	}
	lb_embedding.close();
	ub_embedding.close();
	mean_embedding.close();
}
