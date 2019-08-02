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
 * @file median_letter_demo.cpp
 * @brief GEDLIB implementation of median graph computation for Letter (HIGH) graphs.
 * @details Illustrates jow to use GEDLIB for implementing the algorithm suggested in:
 * - N. Boria, S. Bougleux, B. Ga&uuml;z&egrave;re, L. Brun:
 *   &ldquo;Generalized median graph via iterative alternate minimizations&rdquo;
 *   https://doi.org/10.1007/978-3-030-20081-7_10
 */

/* Use GEDLIB as shared library for comparing graphs in GXL file format.
 */
#define GXL_GEDLIB_SHARED

/* Include the main header file.
 *
 * This is the only header you have to include for using GEDLIB.
 */
#include "../../src/env/ged_env.hpp"

/* Threshold used for comparisons.
 */
extern const double epsilon{0.0001};

/* Saves a Letter graph given as ged::ExchangeGraph as a GXL file.
 *
 * Helper function used by main().
 */
void save_letter_graph_as_gxl_file(ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & median, const std::string & gxl_file_name) {
	std::ofstream gxl_file(gxl_file_name.c_str());
	gxl_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	gxl_file << "<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n";
	gxl_file << "<gxl xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
	gxl_file << "<graph id=\"Z_HIGH_median\" edgeids=\"false\" edgemode=\"undirected\">\n";
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		gxl_file << "<node id=\"_" << i << "\">";
		gxl_file << "<attr name=\"x\"><float>" << median.node_labels.at(i).at("x") << "</float></attr>";
		gxl_file << "<attr name=\"y\"><float>" << median.node_labels.at(i).at("y") << "</float></attr>";
		gxl_file << "</node>\n";
	}
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		for (std::size_t j{i + 1}; j < median.num_nodes; j++) {
			if (median.adj_matrix[i][j] == 1) {
				gxl_file << "<edge from=\"_" << i << "\" to=\"_" << j << "\"/>\n";
			}
		}
	}
	gxl_file << "</graph></gxl>\n";
	gxl_file.close();
}

/* Saves a Letter graph given as ged::ExchangeGraph as a TikZ file.
 *
 * Helper function used by main().
 */
void save_letter_graph_as_tikz_file(ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & median, const std::string & tikz_file_name) {
	std::ofstream tikz_file(tikz_file_name.c_str());
	tikz_file << "\\documentclass[tikz]{standalone}\n";
	tikz_file << "\\begin{document}\n";
	tikz_file << "\\begin{tikzpicture}\n";
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		tikz_file << "\\node[anchor=base,fill,circle,inner sep=1pt] (" << i << ") at (" << median.node_labels.at(i).at("x") << "," << median.node_labels.at(i).at("y") << ") {};\n";
	}
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		for (std::size_t j{i + 1}; j < median.num_nodes; j++) {
			if (median.adj_matrix[i][j] == 1) {
				tikz_file << "\\draw (" << i << ") -- (" << j << ");\n";
			}
		}
	}
	tikz_file << "\\end{tikzpicture}\n";
	tikz_file << "\\end{document}\n";
	tikz_file.close();
}

/* Computes the geometric median (used for updating the node labels).
 *
 * Helper function used by update_median_graph(). Implements Weiszfeld's algorithm.
 */
ged::GXLLabel compute_geometric_median(const ged::GXLLabel & old_label, const std::vector<ged::GXLLabel> & node_labels) {

	/* Compute mean label to be used as initial solution and transform labels into ccordinates.
	 */
	std::vector<std::pair<double, double>> node_labels_as_coords;
	double label_x, label_y;
	double sum_x{0};
	double sum_y{0};
	for (const auto & node_label : node_labels) {
		label_x = std::stod(node_label.at("x"));
		label_y = std::stod(node_label.at("y"));
		sum_x += label_x;
		sum_y += label_y;
		node_labels_as_coords.emplace_back(label_x, label_y);
	}
	std::pair<double, double> median(sum_x / static_cast<double>(node_labels.size()), sum_y / static_cast<double>(node_labels.size()));

	/* Run main while-loop of Weiszfeld's Algorithm.
	 */
	double delta{1.0};
	std::size_t num_itrs{0};
	bool all_equal{false};
	while ((delta > epsilon) and (num_itrs++ < 100) and (not all_equal)) {
		std::pair<double, double> numerator(0, 0);
		double denominator{0};
		for (const auto & node_label_as_coord : node_labels_as_coords) {
			double norm{0};
			norm += (node_label_as_coord.first - median.first) * (node_label_as_coord.first - median.first);
			norm += (node_label_as_coord.second - median.second) * (node_label_as_coord.second - median.second);
			norm = std::sqrt(norm);
			if (norm > 0) {
				numerator.first += node_label_as_coord.first / norm;
				numerator.second += node_label_as_coord.second / norm;
				denominator += 1.0 / norm;
			}
		}
		if (denominator == 0) {
			all_equal = true;
		}
		else {
			std::pair<double, double> new_median(numerator.first / denominator, numerator.second / denominator);
			delta = std::abs(median.first - new_median.first) + std::abs(median.second - new_median.second);
			median = new_median;
		}
	}

	/* Transform solution to ged::GXLLabel and return it.
	 *
	 * Only return the new node label if it is significantly different from the old node label. Otherwise, return the old
	 * node label.
	 */
	delta = std::abs(median.first - std::stod(old_label.at("x"))) + std::abs(median.second - std::stod(old_label.at("y")));
	if (delta > epsilon) {
		ged::GXLLabel median_label;
		median_label["x"] = std::to_string(median.first);
		median_label["y"] = std::to_string(median.second);
		return median_label;
	}
	return old_label;
}


/* Update the median graph.
 *
 * Helper function used by main(). For details, consult the paper by Boria et al.
 */
bool update_median_graph(ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & median,
		const std::vector<ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>> & graphs,
		const std::vector<ged::NodeMap> & node_maps, const std::vector<ged::GEDGraph::GraphID> & graph_ids) {

	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> old_median(median);

	/* Update the node labels of the median.
	 *
	 * For details, consult the paper by Boria et al.
	 */
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		std::vector<ged::GXLLabel> labels_images;
		for (auto graph_id : graph_ids) {
			ged::GEDGraph::NodeID k{node_maps.at(graph_id).image(i)};
			if (k != ged::GEDGraph::dummy_node()) {
				labels_images.emplace_back(graphs.at(graph_id).node_labels.at(k));
			}
		}
		if (labels_images.size() > 0) {
			median.node_labels[i] = compute_geometric_median(median.node_labels.at(i), labels_images);
		}
	}

	/* Update the (unlabeled) edges of the median.
	 *
	 * For details, consult the paper by Boria et al.
	 */
	median.edge_labels.clear();
	ged::GXLLabel dummy_edge_label;
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		for (std::size_t j{0}; j < median.num_nodes; j++) {
			std::size_t num_subs{0};
			for (auto graph_id : graph_ids) {
				ged::GEDGraph::NodeID k{node_maps.at(graph_id).image(i)};
				ged::GEDGraph::NodeID l{node_maps.at(graph_id).image(j)};
				if (k != ged::GEDGraph::dummy_node() and l != ged::GEDGraph::dummy_node()) {
					num_subs += graphs.at(graph_id).adj_matrix.at(k).at(l);
				}
			}
			if (2 * num_subs >= graph_ids.size()) {
				median.adj_matrix[i][j] = 1;
				median.adj_matrix[j][i] = 1;
				median.edge_labels[std::make_pair(i, j)] = dummy_edge_label;
				median.edge_labels[std::make_pair(j, i)] = dummy_edge_label;
			}
			else {
				median.adj_matrix[i][j] = 0;
				median.adj_matrix[j][i] = 0;
			}
		}
	}

	return not (median == old_median);

}

/* The main method.
 */
int main(int argc, char* argv[]) {

	/* Initialize the environment.
	 *
	 * Since the graphs are given as GXL files, we use ged::GXLNodeID a.k.a. std::string
	 * as type for the node IDs and ged::GXLLabel a.k.a. std::map<std::string, std::string>
	 * as type for the node and edge labels. It does not matter that the edges are unlabeled,
	 * we can model this by using an empty ged::GXLLabel dummy_label.
	 */
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;

	/* Use edit costs for Letter graphs defined in ged::Letter<GXLLabel, GXLLabel>.
	 *
	 * Like most predefined edit costs, the edit cost functions for the Letter graphs
	 * are parameterized via constants. In the case of Letter, there are three of them:
	 * the cost of inserting or deleting a node (node_ins_del_cost, default = 0.9),
	 * the cost of inserting or deleting an edge (edge_ins_del_cost, default = 1.7),
	 * and the importance of node edit operations w.r.t. edge edit operations (alpha, default = 0.75).
	 * You can put the constants to non-default values (e.g. node_ins_del_cost = 0.7,
	 * edge_ins_del_cost = 1.5, and alpha = 0.25) by using initializer lists like so:
	 *
	 * env.set_edit_costs(ged::Options::EditCosts::Letter, {0.7, 1.5, 0.25});
	 *
	 * You can also define your own edit costs. For this, define your own implementation
	 * CustomEditCosts<ged::UserNodeLabel, ged::UserEdgeLabel> of the abstract class
	 * ged::EditCosts<ged::UserNodeLabel, ged::UserEdgeLabel>. Assume that your custom
	 * edit costs class has a constructor without arguments. Then you can employ your
	 * edit costs like so:
	 *
	 * CustomEditCosts<ged::GXLLabel, ged::GXLLabel> custom_edit_costs;
	 * env.set_edit_costs(&custom_edit_costs);
	 */
	env.set_edit_costs(ged::Options::EditCosts::LETTER);

	/* Paths to an XML file that specifies which graphs should be loaded and to the directory that contains them.
	 *
	 * The graphs must be given as GXL files, and the collection file must respect
	 * the document type {GEDLIB_ROOT}/data/collections/GraphCollection.dtd.
	 */
	std::string letter_class("A");
	if (argc > 1) {
		letter_class = std::string(argv[1]);
	}
	std::string collection_file("../collections/Letter_" + letter_class + ".xml");
	std::string graph_dir("../../data/datasets/Letter/HIGH/");

	/* Load the GXL graphs into the environment.
	 *
	 * Letter (HIGH) contains graphs with labeled nodes and labeled edges, which we specify with
	 * the third and fourth argument. The method returns a vector with the IDs of all graphs
	 * that are contained in the environment. It is sometimes the case that the edit cost functions
	 * do not use all of the attributes given in the GXL files. In this case, you can tell GEDLIB
	 * to ignore the irrelevant attributes by calling env.load_gxl_graphs() with a fifth argument
	 * const std::unordered_set<std::string> & irrelevant_node_attributes and a sixth argument
	 * const std::unordered_set<std::string> & irrelevant_edge_attributes. For instance, if your
	 * edit cost functions for the Letter graphs do not use the "x" attribute of the node labels, you
	 * can tell GEDLIB to ignore it like so:
	 *
	 * std::unordered_set<std::string> irrelevant_node_attributes, irrelevant_edge_attributes;
	 * irrelevant_node_attributes.insert("x");
	 * std::vector<ged::GEDGraph::GraphID> graph_ids = env.load_gxl_graphs(graph_dir, collection_file,
	 *		ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::UNLABELED,
	 *		irrelevant_node_attributes, irrelevant_edge_attributes);
	 */
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::UNLABELED));

	/* Add empty graph to be used later for median.
	 *
	 * It is recommendable to allocate space for all graphs that have to be added later on.
	 * The reason is that whenever a new graph is added after initialization, all graphs contained
	 * in the environment have to be re-initialized.
	 */
	ged::GEDGraph::GraphID median_id{env.add_graph("median_Letter_HIGH_" + letter_class, letter_class)};


	/* Initialize the environment.
	 *
	 * There are four initialization types: EAGER_WITHOUT_SHUFFLED_COPIES, EAGER_WITH_SHUFFLED_COPIES,
	 * LAZY_WITHOUT_SHUFFLED_COPIES, and LAZY_WITH_SHUFFLED_COPIES. Eager initialization means that
	 * all edit costs are pre-computed and stored in a matrix. If lazy initialization is selected,
	 * the edit cost functions are evaluated on the fly whenever they are needed. This can slow down
	 * the methods significantly, so whenever you have enough memory for doing so, you should use
	 * eager initialization. If you choose initialization with shuffled copies, shuffled copies of all
	 * graphs are constructed which are used whenever a graph is compared to itself.
	 */
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

	/* Select the method that should be used for (approximately) computing GED.
	 *
	 * Most methods can be used with various options that can be specified via strings of the form
	 * "[--<option> <arg>] [...]". For example, we here tell GEDLIB that we want to run IPFP in six threads
	 * from 40 different randomly generated initial solutions. Moreover, some methods can optionally be
	 * initialized such that, at runtime, some parts have already been pre-computed.
	 * This can be done by calling env.init_method() after having called env.set_method().
	 */
	std::string ipfp_options("--threads 6 --initial-solutions 5 --initialization-method RANDOM");
	env.set_method(ged::Options::GEDMethod::IPFP, ipfp_options);

	/* Get all graphs that have been loaded into the environment in easily managable format.
	 *
	 * The structure ged::ExchangeGraph<ged::UserNodeID, ged::UserNodeLabel, ged::UserEdgeLabel>
	 * allows the user to easily inspect and modify the graphs contained in the environment.
	 * Loading them into memory is necessary only if you have to inspect the graphs after
	 * initializing the environment.
	 */
	std::vector<ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>> graphs;
	for (auto graph_id : graph_ids) {
		graphs.emplace_back(env.get_graph(graph_id));
	}

	/* Compute the set median.
	 */
	ged::GEDGraph::GraphID set_median_id{0};
	std::vector<double> sums_of_distances({std::numeric_limits<double>::infinity()});
	ged::ProgressBar progress(graph_ids.size());
	std::cout << "\rComputing set median: " << progress << std::flush;
	for (auto g_id : graph_ids) {
		double sum_dists{0};
		for (auto h_id : graph_ids) {
			env.run_method(g_id, h_id);
			sum_dists += env.get_upper_bound(g_id, h_id);
		}
		if (sum_dists < sums_of_distances.at(0)) {
			sums_of_distances[0] = sum_dists;
			set_median_id = g_id;
		}
		progress.increment();
		std::cout << "\rComputing set median: " << progress << std::flush;
	}
	std::cout << "\n";

	/* Get the node maps from the set median to all other graphs.
	 */
	std::vector<ged::NodeMap> node_maps;
	for (auto graph_id : graph_ids) {
		node_maps.emplace_back(env.get_node_map(set_median_id, graph_id));
	}

	/* Get the ExchangeGraph representation of the set median and set the median ID to the ID of a newly added graph.
	 *
	 * Since we have called env.add_graph(), the environment is again uninitialized. So we have to keep in mind that
	 * we must call env.init() again before calling env.run_method().
	 */
	ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median(graphs.at(set_median_id));
	median.id = median_id;


	/* Save the set median as GXL file and as TikZ file.
	 */
	std::string gxl_file_name("../output/set_median_Letter_HIGH_" + letter_class + ".gxl");
	save_letter_graph_as_gxl_file(median, gxl_file_name);
	std::string tikz_file_name("../output/set_median_Letter_HIGH_" + letter_class + ".tex");
	save_letter_graph_as_tikz_file(median, tikz_file_name);


	/* Main while loop.
	 */
	bool median_was_modified{true};
	bool node_maps_were_modified{true};
	std::size_t num_itrs{0};
	std::cout << "Running alternate minimization:\n";
	while (median_was_modified or node_maps_were_modified) {

		/* Update the median graph.
		 */
		median_was_modified = update_median_graph(median, graphs, node_maps, graph_ids);

		/* Load modified median into the environment.
		 *
		 * We have to re-initialize the environment after calling env.clear_graph() or env.load_exchange_graph().
		 */
		env.load_exchange_graph(median, median_id);
		env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

		/* Compute induced costs of old node maps w.r.t. the updated median.
		 */
		for (auto graph_id : graph_ids) {
			env.compute_induced_cost(median_id, graph_id, node_maps.at(graph_id));
		}

		/* Update the node maps.
		 */
		node_maps_were_modified = false;
		for (auto graph_id : graph_ids) {
			env.run_method(median_id, graph_id);
			if (env.get_upper_bound(median_id, graph_id) < node_maps.at(graph_id).induced_cost() - epsilon) {
				node_maps[graph_id] = env.get_node_map(median_id, graph_id);
				node_maps_were_modified = true;
			}
		}

		/* Compute sum of distances of the current median and the current node maps.
		 */
		sums_of_distances.emplace_back(0);
		std::cout << "\tIteration " << ++num_itrs << ": ";
		for (auto graph_id : graph_ids) {
			sums_of_distances[num_itrs] += node_maps.at(graph_id).induced_cost();
		}
		if (sums_of_distances[num_itrs] < sums_of_distances[num_itrs - 1] - epsilon) {
			std::cout << "OLD-SOD>NEW-SOD";
		}
		else if (sums_of_distances[num_itrs] > sums_of_distances[num_itrs - 1] + epsilon) {
			std::cout << "OLD-SOD<NEW-SOD";
		}
		else {
			std::cout << "OLD-SOD=NEW-SOD";
		}
		std::cout << ", MEDIAN-MODIFIED=" << median_was_modified << ", NODE-MAPS-MODIFIED=" << node_maps_were_modified << "\n";
	}

	/* Print information about obtained solution.
	 */
	std::cout << "Number of iterations: " << num_itrs << "\n";
	std::cout << "Sum of distances for set median: " << sums_of_distances[0] << "\n";
	std::cout << "Sum of distances for generalized median: " << sums_of_distances[num_itrs] << "\n";

	/* Save the obtained solution as GXL file and as TikZ file.
	 */
	gxl_file_name = "../output/gen_median_Letter_HIGH_" + letter_class + ".gxl";
	save_letter_graph_as_gxl_file(median, gxl_file_name);
	tikz_file_name = "../output/gen_median_Letter_HIGH_" + letter_class + ".tex";
	save_letter_graph_as_tikz_file(median, tikz_file_name);

}

