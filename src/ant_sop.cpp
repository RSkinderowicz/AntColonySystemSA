#include "ant_sop.h"
#include <iostream>

using namespace std;


AntSOP::AntSOP(const SOPInstance &problem) :
    Ant(problem.get_dimension()),
    problem_(problem)
{
    const auto n = problem_.get_dimension();

    incoming_count_.resize(n);
    const auto &c = problem_.get_incoming_edges_count();
    ::copy(c.begin(), c.end(), incoming_count_.begin());

    available_.resize(n);
    for (auto u = 0u; u < n; ++u) {
        available_[u] = (incoming_count_[u] == 0);
    }
}


void AntSOP::move_to(uint32_t node) {
    // The first node in case of the SOP is always the source, i.e. node 0
    if (visited_.empty()) {
        node = 0;
    }

    assert(available_[node]);
    assert(incoming_count_[node] == 0);

    Ant::move_to(node);
    available_[node] = false;
    // Now decrease all incoming edges count starting with node v
    for (auto u : problem_.get_outgoing_edges(node)) {
        if (--incoming_count_[u] == 0) {
            available_[u] = true;
        }
    }
}
