/*!
 * @file train_ring.cpp
 * @brief Trains the method ged::Ring for different datasets.
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

	// Learn the parameters.
	std::vector<std::string> led_methods{"GAMMA", "LSAPE_GREEDY", "LSAPE_OPTIMAL"};
	for (auto led_method : led_methods) {
		std::cout << "\n=== " << led_method << " ===\n";
		env.set_method(ged::Options::GEDMethod::RING, util::init_options(dataset, std::string("ring_") + led_method) +  " --led-method " + led_method);
		env.init_method();
	}
}

int main(int argc, char* argv[]) {
	std::vector<std::string> datasets;
	util::setup_datasets(datasets);
	for (auto dataset : datasets) {
		if (dataset == "Letter_HIGH" or dataset == "pah" or dataset == "AIDS") { // todo rm
			continue;					// todo rm
		}								// todo rm
		try {
			train_on_dataset(dataset);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << ".\n";
		}
	}
	return 0;
}
