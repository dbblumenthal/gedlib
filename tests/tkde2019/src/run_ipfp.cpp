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

#include "util.hpp"

void run_on_dataset(const std::string & dataset) {
    ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
    util::setup_environment(dataset, false, env);
    env.set_method(ged::Options::GEDMethod::IPFP, "--threads 10 --initial-solutions 40 --ratio-runs-from-initial-solutions 0.25 --num-randpost-loops 3");
    env.init_method();
    std::size_t num_runs{(env.graph_ids().second * env.graph_ids().second) - env.graph_ids().second};
    ged::ProgressBar progress_bar(num_runs);
    std::cout << "\rRunning IPFP on " << dataset << ": " << progress_bar << std::flush;
    ged::GEDGraph::GraphID closest_graph_id{std::numeric_limits<ged::GEDGraph::GraphID>::max()};
    double distance_to_closest_graph{std::numeric_limits<double>::infinity()};
    double avg_runtime{0};
    double avg_ub{0};
    double avg_classification_ratio{0};
    for (ged::GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++) {
        closest_graph_id = std::numeric_limits<ged::GEDGraph::GraphID>::max();
        distance_to_closest_graph = std::numeric_limits<double>::infinity();
        for (ged::GEDGraph::GraphID h_id = env.graph_ids().first; h_id != env.graph_ids().second; h_id++) {
            if (g_id == h_id) {
                continue;
            }
            env.run_method(g_id, h_id);
            avg_ub += env.get_upper_bound(g_id, h_id);
            avg_runtime += env.get_runtime(g_id, h_id);
            if (env.get_upper_bound(g_id, h_id) < distance_to_closest_graph) {
                distance_to_closest_graph = env.get_upper_bound(g_id, h_id);
                closest_graph_id = h_id;
            }
            progress_bar.increment();
            std::cout << "\rRunning IPFP on " << dataset << ": " << progress_bar << std::flush;
        }
        if (env.get_graph_class(g_id) == env.get_graph_class(closest_graph_id)) {
            avg_classification_ratio += 1.0;
        }
    }
    avg_ub /= static_cast<double>(num_runs);
    avg_runtime /= static_cast<double>(num_runs);
    avg_classification_ratio /= static_cast<double>(env.graph_ids().second);
    std::string result_filename("../output/");
    result_filename += dataset + "__IPFP.csv";
    std::ofstream result_file(result_filename.c_str());
    result_file << "avg_ub,avg_runtime,avg_classification_ratio\n";
    result_file << avg_ub << "," << avg_runtime << "," << avg_classification_ratio << "\n";
    result_file.close();
}

int main(int argc, char* argv[]) {
    std::vector<std::string> datasets{"Letter_HIGH", "pah", "alkane", "S-MOL_NL01", "S-MOL_NL04",
                                      "S-MOL_NL07", "S-MOL_NL10", "S-mao_NL03", "S-mao_NL05", "S-mao_NL07",
                                      "S-mao_NL09", "S-acyclic_NL03", "S-acyclic_NL05", "S-acyclic_NL07",
                                      "S-acyclic_NL09", "AIDS"};
    for (const std::string & dataset : datasets) {
        run_on_dataset(dataset);
    }

}
