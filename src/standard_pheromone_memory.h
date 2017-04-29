#ifndef STANDARD_PHEROMONE_MEMORY_H
#define STANDARD_PHEROMONE_MEMORY_H

#include <vector>
#include <cstdint>
#include <cassert>
#include "rng.h"
#include "pheromone_memory.h"

/**
 * An implementation of the standard ACO pheromone memory using a matrix as a
 * storage.
 */
class StandardPheromoneMemory : public PheromoneMemory {
public:

    StandardPheromoneMemory(uint32_t size, double initial_pheromone, bool is_symmetric);


    void reset() override;


    inline double get(uint32_t i, uint32_t j) const final {
        assert(i < size_);
        assert(j < size_);
        return trails_[i][j];
    }


    void set(uint32_t i, uint32_t j, double value);

    /*
     * Updates the value according to the given evap_ratio and returns
     * true if the value has changed.
     */
    bool update_pheromone(uint32_t u, uint32_t v, double evap_ratio, double pheromone_delta) override;


    const std::vector<double>& get_trails(uint32_t node) const {
        return trails_.at(node);
    }

private:

    std::vector<std::vector<double>> trails_;
    bool is_symmetric_ = true;
};


#endif /* ifndef STANDARD_PHEROMONE_MEMORY_H */
