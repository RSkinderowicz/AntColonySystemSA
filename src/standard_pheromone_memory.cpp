#include "standard_pheromone_memory.h"
#include <cassert>

using namespace std;

/**
 * Fill pheromone memory with initial_pheromone values.
 */
StandardPheromoneMemory::StandardPheromoneMemory(uint32_t size, double initial_pheromone, bool is_symmetric)
    : PheromoneMemory(size, initial_pheromone, is_symmetric) {

    assert(size > 2);
    assert(initial_pheromone > 1.0e-20);

    trails_.resize(size_);
    for (auto &row : trails_) {
        row.resize(size_, initial_pheromone_);
    }
    for (auto i = 0u; i < size_; ++i) {
        trails_[i][i] = 0.0;
    }
}


/**
 * Fills all pheromone trails with initial_pheromone value.
 */
void StandardPheromoneMemory::reset() {
    PheromoneMemory::reset();
    for (auto &row : trails_) {
        fill(begin(row), end(row), initial_pheromone_);
    }
    for (auto i = 0u; i < size_; ++i) {
        trails_[i][i] = 0.0;
    }
}


void StandardPheromoneMemory::set(uint32_t i, uint32_t j, double value) {
    assert(i < size_);
    assert(j < size_);

    trails_[i][j] = value;
    if (is_symmetric_) {
        trails_[j][i] = value;
    }
}


/*
 * Updates the value according to the given evap_ratio and returns
 * true if the value has changed.
 */
bool StandardPheromoneMemory::update_pheromone(uint32_t u, uint32_t v, double evap_ratio, double pheromone_delta) {
    auto old = get(u, v);
    auto updated = old * (1 - evap_ratio) + evap_ratio * pheromone_delta;
    set(u, v, updated);
    return (updated != old);
}
