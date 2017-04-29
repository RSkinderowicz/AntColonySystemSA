#ifndef LRU_PHEROMONE_MEMORY_H
#define LRU_PHEROMONE_MEMORY_H

#include "rng.h"
#include "pheromone_memory.h"


/**
 * Implementation of the selective pheromone memory (SPM) in which the least
 * recently added pheromone trail is replaced with a new one if the limit per
 * node is * reached.
 */
class LRUPheromoneMemory : public PheromoneMemory {
public:
    // Pheromone trails are stored in buckets which are linear sequences of
    // (pheromone, node) pairs
    using Bucket = std::vector<std::pair<uint32_t, double>>;


    LRUPheromoneMemory(uint32_t size, double initial_pheromone, bool is_symmetric,
                       uint32_t bucket_size = 8);


    void reset() override;


    double get(uint32_t i, uint32_t j) const override;


    void set(uint32_t i, uint32_t j, double value);


    bool update_pheromone(uint32_t u, uint32_t v,
                          double evap_ratio, double pheromone_delta) override {
        double updated = get(u, v) * (1 - evap_ratio) + evap_ratio * pheromone_delta;
        set(u, v, updated);
        return true;
    }


    void forget(uint32_t count, RNG &rng);


    double get_initial_pheromone() const { return initial_pheromone_; }


    const Bucket &get_bucket(uint32_t index) const {
        return buckets_.at(index);
    }


    void mix(PheromoneMemory *other, RNG &rng) override ;

private:

    void set_helper(uint32_t u, uint32_t v, double value) noexcept;

    const uint32_t bucket_size_;
    std::vector<Bucket> buckets_;
};

#endif /* ifndef SPM_LRU_PHEROMONE_STORAGE_H */
