
#include "catch.hpp"

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"

TEST_CASE("testing on MAO graphs") {
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(ged::Options::EditCosts::CHEM_1);
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs("../../../data/datasets/mao/", "../collections/mao_5.xml"));
	ged::GEDGraph::GraphID g {graph_ids[0]};
	ged::GEDGraph::GraphID h {graph_ids[0]};
	//env.init();
	env.init(ged::Options::InitType::EAGER_WITH_SHUFFLED_COPIES);
	double lower_bound{0.0};
	double exact{0.0};


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

	SECTION("Subgraph") {
		std::cout << "\n=== Subgraph ===\n";
		env.set_method(ged::Options::GEDMethod::SUBGRAPH, "--threads 1 --time-limit-subproblem 0.0001");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 1 --time-limit-subproblem 0.0001\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::SUBGRAPH, "--threads 1 --time-limit-subproblem 0.001");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --time-limit-subproblem 0.001\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
		env.set_method(ged::Options::GEDMethod::SUBGRAPH, "--threads 1 --time-limit-subproblem 0.01");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --time-limit-subproblem 0.01\"\n";
		std::cout << "upper bound = " << env.get_upper_bound(g, h) << ", runtime = " << env.get_runtime(g, h) << "\n";
		CHECK(lower_bound <= env.get_upper_bound(g, h));
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

	SECTION("Exact") {
		std::cout << "\n===running Exact ===\n";
		env.set_method(ged::Options::GEDMethod::EXACT, "--threads 5 --time-limit 1.0");
		env.run_method(g, h);
		std::cout << "\noptions = \"--threads 5 --time-limit 1.0\"\n";
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

	/*SECTION("Ring") {
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




