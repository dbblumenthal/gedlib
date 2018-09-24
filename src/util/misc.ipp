/*!
 * 	@file  misc.ipp
 *  @brief Definition of miscellaneous utility functions.
 */

#ifndef SRC_UTIL_MISC_IPP_
#define SRC_UTIL_MISC_IPP_

namespace ged {

namespace util {

void
init_adj_matrix(const GEDGraph & graph, const GEDGraph::SizeTNodeMap & ids_to_nodes, DMatrix & adj_matrix) {
	for (std::size_t row{0}; row < adj_matrix.num_rows(); row++) {
		for (std::size_t col{0}; col < adj_matrix.num_cols(); col++) {
			if (graph.is_edge(ids_to_nodes.at(row), ids_to_nodes.at(col))) {
				adj_matrix(row, col) = 1.0;
			}
			else {
				adj_matrix(row, col) = 0.0;
			}
		}
	}
}

void
parse_config_file(const std::string & filename, std::map<std::string, std::string> & options) {
	std::ifstream config_file(filename);
	std::string line;
	std::size_t line_nr{1};
	while(std::getline(config_file, line)) {
		if (line.at(0) == '#') {
			continue;
		}
		std::string error_msg("Line " + std::to_string(line_nr) + " has invalid format.\nExpected format: \"<key>=<value>\"\nLine " + std::to_string(line_nr) + ": " + line);
		std::istringstream line_stream(line);
		std::string key;
		if (not std::getline(line_stream, key, '=')) {
			throw Error(error_msg);
		}
		std::string value;
		if(not std::getline(line_stream, value)) {
			throw Error(error_msg);
		}
		else {
			options[key] = value;
		}
		line_nr++;
	}
	config_file.close();
}

void
save_as_config_file(const std::string & filename, const std::map<std::string, std::string> & options) {
	std::ofstream config_file(filename);
	for (auto key_value = options.begin(); key_value != options.end(); key_value++) {
		config_file << key_value->first << "=" << key_value->second << "\n";
	}
	config_file.close();
}

template<class Solver>
void
construct_node_map_from_solver(const Solver & solver, const GEDGraph::SizeTNodeMap & g_ids_to_nodes,
		const GEDGraph::SizeTNodeMap & h_ids_to_nodes, NodeMap & matching, std::size_t solution_id) {
	matching.clear();
	std::size_t num_nodes_g{g_ids_to_nodes.size()};
	std::size_t num_nodes_h{h_ids_to_nodes.size()};

	// add deletions and substitutions
	for (std::size_t row{0}; row < num_nodes_g; row++) {
		if (solver.get_assigned_col(row, solution_id) >= num_nodes_h) {
			matching.add_assignment(g_ids_to_nodes.at(row), GEDGraph::dummy_node());
		}
		else {
			matching.add_assignment(g_ids_to_nodes.at(row), h_ids_to_nodes.at(solver.get_assigned_col(row, solution_id)));
		}
	}

	// insertions
	for (std::size_t col{0}; col < num_nodes_h; col++) {
		if (solver.get_assigned_row(col, solution_id) >= num_nodes_g) {
			matching.add_assignment(GEDGraph::dummy_node(), h_ids_to_nodes.at(col));
		}
	}
}

void
init_node_to_id_indices(const GEDGraph & graph, std::map<GEDGraph::GraphID, GEDGraph::NodeSizeTMap> & nodes_to_ids) {
	nodes_to_ids[graph.id()] = GEDGraph::NodeSizeTMap();
	init_node_to_id_indices(graph, nodes_to_ids[graph.id()]);
}

void
init_node_to_id_indices(const GEDGraph & graph, GEDGraph::NodeSizeTMap & nodes_to_ids) {
	std::size_t id{0};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		nodes_to_ids[*node] = id++;
	}
}

void
init_id_to_node_indices(const GEDGraph & graph, std::map<GEDGraph::GraphID, GEDGraph::SizeTNodeMap> & ids_to_nodes) {
	ids_to_nodes[graph.id()] = GEDGraph::SizeTNodeMap();
	init_id_to_node_indices(graph, ids_to_nodes[graph.id()]);
}

void
init_id_to_node_indices(const GEDGraph & graph, GEDGraph::SizeTNodeMap & ids_to_nodes) {
	std::size_t id{0};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		ids_to_nodes[id++] = *node;
	}
}

void
counting_sort(std::vector<LabelID>::iterator first, std::vector<LabelID>::iterator last) {
	// Find maximum label value and range.
	LabelID max_label_val{0};
	std::size_t range{0};
	for (auto label = first; label != last; label++) {
		max_label_val = std::max(max_label_val, *label);
		range++;
	}

	// Compute histograms that store the number of labels in input for each label value.
	std::vector<LabelID> hist(max_label_val + 1, 0);
	for (auto label = first; label != last; label++) {
		hist[*label]++;
	}

	// Compute starting position for each label value;
	std::vector<LabelID> pos(max_label_val + 1, 0);
	for (std::size_t label_val{0}; label_val < max_label_val; label_val++) {
		pos[label_val + 1] = pos.at(label_val) + hist.at(label_val);
	}

	// Compute sorted label vector.
	std::vector<LabelID> sorted_labels(range);
	for (auto label = first; label != last; label++) {
		sorted_labels[pos[*label]++] = *label;
	}

	// Write sorted label vector into input.
	for (auto label = sorted_labels.begin(); label != sorted_labels.end(); label++) {
		*first++ = *label;
	}
}

}

}

#endif /* SRC_UTIL_MISC_IPP_ */
