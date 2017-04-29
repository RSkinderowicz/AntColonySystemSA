#ifndef PHEROMONE_MEMORY_H
#define PHEROMONE_MEMORY_H

#include <iostream>
#include <memory>
#include <cassert>
#include "rng.h"

/*
 * Base pheromone memory class with almost all methods abstract and no specific
 * storage (i.e. matrix) for the pheormone values specified.
*/
class PheromoneMemory {
public:

    PheromoneMemory(uint32_t size, double initial_pheromone, bool is_symmetric)
        : size_(size),
          initial_pheromone_(initial_pheromone),
          is_symmetric_(is_symmetric) {
    }

    virtual ~PheromoneMemory() {}

    virtual void reset() {}


    virtual double get(uint32_t i, uint32_t j) const = 0;


    virtual bool update_pheromone(uint32_t u, uint32_t v, double evap_ratio, double pheromone_delta) = 0;


    virtual bool local_update(const std::vector<uint32_t> &visited,
                              uint32_t node_index, double evap_ratio) {
        auto u = visited[node_index];
        auto v = visited[ (node_index + 1) % visited.size() ];
        return update_pheromone(u, v, evap_ratio, initial_pheromone_);
    }


    /*
     * Performs global pheromone update.
     *
     * If 'changed' container is specified then a list of edges (node pairs) is
     * returned.
     */
    virtual void global_update(const std::vector<uint32_t> &visited,
                               double route_length,
                               double evap_ratio,
                               std::vector<uint32_t> *changed = nullptr) {
        auto delta = 1.0 / route_length;
        auto &nodes = visited;
        auto n = nodes.size();
        for (auto i = 0u; i < n; ++i) {
            auto u = nodes[i];
            auto v = nodes[ (i+1) % n ];
            update_pheromone(u, v, evap_ratio, delta);
            if (changed != nullptr) {
                changed->push_back(u);
                changed->push_back(v);
            }
        }
    }


    void print_debug() {
        std::cout << "Initial pheromone: " << initial_pheromone_ << std::endl;
    }


    virtual void mix(PheromoneMemory *mem, RNG &rng) {
        throw std::runtime_error("PheromoneMemory::mix - Not implemented yet!");
    }

protected:

    uint32_t size_;
    double initial_pheromone_;
    bool is_symmetric_;
};

#endif /* ifndef PHEROMONE_MEMORY_H */
