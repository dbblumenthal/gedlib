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
 * @file test_fgw_ged.cpp
 * @brief
 */

#define GXL_GEDLIB_SHARED
#define FGW_GEDLIB_SHARED

#include "../../../src/env/ged_env.hpp"

void determine_dir_collection_edit_costs(const std::string & dataset, std::string & dir, std::string & collection, ged::Options::EditCosts & edit_costs) {
	if (dataset == "Letter") {
		dir = std::string("../../../data/datasets/Letter/HIGH/");
		collection = std::string("../collections/Letter_100.xml");
		edit_costs = ged::Options::EditCosts::LETTER;
	}
	else if (dataset == "pah") {
		dir = std::string("../../../data/datasets/pah/");
		collection = std::string("../collections/pah.xml");
		edit_costs = ged::Options::EditCosts::CHEM_2;
	}
	else if (dataset == "AIDS") {
		dir = std::string("../../../data/datasets/AIDS/data/");
		collection = std::string("../collections/AIDS_100.xml");
		edit_costs = ged::Options::EditCosts::CHEM_2;
	}
	else if (dataset == "GREC") {
		dir = std::string("../../../data/datasets/GREC/data/");
		collection = std::string("../collections/GREC_100.xml");
		edit_costs = ged::Options::EditCosts::GREC_2;
	}
	else if (dataset == "Mutagenicity") {
		dir = std::string("../../../data/datasets/Mutagenicity/data/");
		collection = std::string("../collections/Mutagenicity_100.xml");
		edit_costs = ged::Options::EditCosts::CHEM_2;
	}
	else if (dataset == "Protein") {
		dir = std::string("../../../data/datasets/Protein/data/");
		collection = std::string("../collections/Protein_100.xml");
		edit_costs = ged::Options::EditCosts::PROTEIN;
	}
	else if (dataset == "Fingerprint") {
		dir = std::string("../../../data/datasets/Fingerprint/data/");
		collection = std::string("../collections/Fingerprint_100.xml");
		edit_costs = ged::Options::EditCosts::FINGERPRINT;
	}
	else {
		throw ged::Error("Invalid dataset.");
	}
}

void test_ged_on_dataset(const std::string & dataset, double & ratio, double & runtime) {
	std::string dir;
	std::string collection;
	ged::Options::EditCosts edit_costs;
	determine_dir_collection_edit_costs(dataset, dir, collection, edit_costs);
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(dir, collection));
	env.set_edit_costs(edit_costs);
	env.init();
	env.set_method(ged::Options::GEDMethod::BIPARTITE, "--threads 10");
	std::size_t num_runs{(graph_ids.size() * graph_ids.size()) - graph_ids.size()};
	ged::ProgressBar progress_bar(num_runs);
	std::cout << "\r\t Standard GED on " << dataset << ": "  << progress_bar << std::flush;
	std::size_t num_correct_classifications{0};
	runtime = 0.0;
	for (ged::GEDGraph::GraphID g_id : graph_ids) {
		double distance_to_closest_graph{std::numeric_limits<double>::infinity()};
		ged::GEDGraph::GraphID closest_graph{std::numeric_limits<ged::GEDGraph::GraphID>::max()};
		for (ged::GEDGraph::GraphID h_id : graph_ids) {
			if (g_id == h_id) {
				continue;
			}
			env.run_method(g_id, h_id);
			runtime += env.get_runtime(g_id, h_id);
			if (env.get_upper_bound(g_id, h_id) < distance_to_closest_graph) {
				distance_to_closest_graph = env.get_upper_bound(g_id, h_id);
				closest_graph = h_id;
			}
			progress_bar.increment();
			std::cout << "\r\t Standard GED on " << dataset << ": "  << progress_bar << std::flush;
		}
		if (env.get_graph_class(g_id) == env.get_graph_class(closest_graph)) {
			num_correct_classifications++;
		}
	}
	std::cout << "\n";
	ratio = static_cast<double>(num_correct_classifications) / static_cast<double>(graph_ids.size());
	runtime /= static_cast<double>(num_runs);
}


void test_fgw_ged_on_dataset(const std::string & dataset, double & ratio, double & runtime) {
	std::string dir;
	std::string collection;
	ged::Options::EditCosts edit_costs;
	determine_dir_collection_edit_costs(dataset, dir, collection, edit_costs);
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, double> env;
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_fgw_gxl_graphs(dir, collection));
	env.set_edit_costs(edit_costs);
	env.init();
	env.set_method(ged::Options::GEDMethod::BIPARTITE, "--threads 10");
	std::size_t num_runs{(graph_ids.size() * graph_ids.size()) - graph_ids.size()};
	ged::ProgressBar progress_bar(num_runs);
	std::cout << "\r\t FGW GED on " << dataset << ": "  << progress_bar << std::flush;
	std::size_t num_correct_classifications{0};
	runtime = 0.0;
	for (ged::GEDGraph::GraphID g_id : graph_ids) {
		double distance_to_closest_graph{std::numeric_limits<double>::infinity()};
		ged::GEDGraph::GraphID closest_graph{std::numeric_limits<ged::GEDGraph::GraphID>::max()};
		for (ged::GEDGraph::GraphID h_id : graph_ids) {
			if (g_id == h_id) {
				continue;
			}
			env.run_method(g_id, h_id);
			runtime += env.get_runtime(g_id, h_id);
			if (env.get_upper_bound(g_id, h_id) < distance_to_closest_graph) {
				distance_to_closest_graph = env.get_upper_bound(g_id, h_id);
				closest_graph = h_id;
			}
			progress_bar.increment();
			std::cout << "\r\t FGW GED on " << dataset << ": "  << progress_bar << std::flush;
		}
		if (env.get_graph_class(g_id) == env.get_graph_class(closest_graph)) {
			num_correct_classifications++;
		}
	}
	std::cout << "\n";
	ratio = static_cast<double>(num_correct_classifications) / static_cast<double>(graph_ids.size());
	runtime /= static_cast<double>(num_runs);
}

int main(int argc, char* argv[]) {
	std::vector<std::string> datasets;
	for (int i{1}; i < argc; i++) {
		datasets.push_back(std::string(argv[i]));
	}
	if (datasets.empty()) {
		datasets = {"Letter", "GREC", "pah", "AIDS", "Protein", "Mutagenicity"};
	}
	std::string result_filename("../output/classification_ratios.csv");
	std::ofstream result_file(result_filename.c_str());
	result_file << "dataset,ratio_fgwd,runtime_fgwd,ratio_ged,runtime_ged\n";
	result_file.close();
	double ratio{0};
	double runtime{0};
	for (auto dataset : datasets) {
		result_file.open(result_filename.c_str(),std::ios_base::app);
		result_file << dataset << ",";
		result_file.close();
		try {
			test_fgw_ged_on_dataset(dataset, ratio, runtime);
			result_file.open(result_filename.c_str(),std::ios_base::app);
			result_file << ratio << "," << runtime << ",";
			result_file.close();
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error computing FGW GED on " << dataset << ".\n";
		}
		try {
			test_ged_on_dataset(dataset, ratio, runtime);
			result_file.open(result_filename.c_str(),std::ios_base::app);
			result_file << ratio << "," << runtime << "\n";
			result_file.close();
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error computing standard GED on " << dataset << ".\n";
		}
	}
	return 0;
}

