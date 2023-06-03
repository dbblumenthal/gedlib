
#include "../env/ged_graph.hpp"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace ged {
    extern "C" double
    compute_substitution_cost_cuda_(const GEDGraph &g, const GEDGraph &h, GEDGraph::NodeID i, GEDGraph::NodeID k) const;

    extern "C" double compute_deletion_cost_cuda_(const GEDGraph &g, GEDGraph::NodeID i) const;

    extern "C" double compute_insertion_cost_cuda_(const GEDGraph &h, GEDGraph::NodeID k) const;
}

extern "C" double ged::compute_substitution_cost_cuda_(const GEDGraph &g, const GEDGraph &h, GEDGraph::NodeID i,
                                                       GEDGraph::NodeID k) const {
    // Collect node substitution costs.
    double cost{this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

    // Initialize subproblem.
    DMatrix subproblem(g.degree(i) + 1, h.degree(k) + 1);

    // Collect edge deletion costs.
    std::size_t j{0};
    for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++, j++) {
        subproblem(j, h.degree(k)) = this->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label()) / 2.0;
    }

    // Collect edge insertion costs.
    std::size_t l{0};
    for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++, l++) {
        subproblem(g.degree(i), l) = this->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl)) / 2.0;
    }
    j = 0;

    // Collect edge relabelling costs.
    for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++, j++) {
        l = 0;
        for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++, l++) {
            subproblem(j, l) = this->ged_data_.edge_cost(g.get_edge_label(*ij), h.get_edge_label(*kl)) / 2.0;
        }
    }

    // Solve subproblem.
    LSAPESolver subproblem_solver(&subproblem);
    subproblem_solver.set_model(this->lsape_model_);
    subproblem_solver.solve();

    // Update and return overall substitution cost.
    cost += subproblem_solver.minimal_cost();
    return cost;
}