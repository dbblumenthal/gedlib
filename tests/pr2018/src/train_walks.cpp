/*!
 * @file train_walks.cpp
 * @brief Trains the method ged::Walks for different datasets.
 * @details The binary built from this file was used for the experiments in the following submission:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Designing heuristics for graph edit distance&rdquo;,
 *   Submitted to PR.
 */

#include "util.hpp"

void train_on_dataset(const std::string & dataset) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	util::setup_environment(dataset, true, env);

	// Initialize the method.
	env.set_method(ged::Options::GEDMethod::WALKS, util::init_options(dataset, "walks"));
	env.init_method();
}

int main(int argc, char* argv[]) {
	std::vector<std::string> datasets;
	util::setup_datasets(datasets);
	for (auto dataset : datasets) {
		try {
			train_on_dataset(dataset);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << ".\n";
		}
	}
	return 0;
}
