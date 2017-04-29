#include "lru_pheromone_memory.h"

#include <iostream>

using namespace std;

LRUPheromoneMemory::LRUPheromoneMemory(uint32_t size,
                                       double initial_pheromone,
                                       bool is_symmetric,
                                       uint32_t bucket_size)
    : PheromoneMemory(size, initial_pheromone, is_symmetric),
      bucket_size_(bucket_size) {
    reset();
}


void LRUPheromoneMemory::reset() {
    buckets_.resize(size_);
    for (auto &el : buckets_) {
        el.clear();
    }
}


double LRUPheromoneMemory::get(uint32_t u, uint32_t v) const {
    assert(u < buckets_.size());

    for (auto &el : buckets_[u]) {
        if (el.first == v) {
            return el.second;
        }
    }
    return initial_pheromone_;
}


void LRUPheromoneMemory::set(uint32_t u, uint32_t v, double value) {
    set_helper(u, v, value);
    if (is_symmetric_) {
        set_helper(v, u, value);
    }
}


/**
 * A helper method to insert a new trail into a bucket for the node v.
 *
 * New trail is inserted at the beginning of the bucket. An updated trail stays
 * at the same place. If the bucket is full, the last recently inserted trail is
 * removed from the bucket.
 */
void LRUPheromoneMemory::set_helper(uint32_t u, uint32_t v, double value) noexcept {
    assert(u < buckets_.size());
    assert(v < buckets_.size());

    auto &bucket = buckets_[u];

    for (auto &el : bucket) {
        if (el.first == v) {
            el.second = value; // Just update
            return ;           // and finish
        }
    }
    if (bucket.size() == bucket_size_) { // Is bucket full?
        bucket.pop_back();
    }
    // Insert new trail at the beginning to prevent it from quick removal
    bucket.insert(bucket.begin(), make_pair(v, value));
}


void LRUPheromoneMemory::forget(uint32_t count, RNG &rng) {
    const auto size = buckets_.size();

    auto sample = rng.random_sample(0u, size-1u, count);

    for (auto i : sample) {
        auto &bucket = buckets_[i];
        const auto m = bucket.size();
        if (m > 0) {
            bucket.pop_back();
        }
    }
}


void LRUPheromoneMemory::mix(PheromoneMemory *mem, RNG &rng) {
    auto &other = dynamic_cast<LRUPheromoneMemory&>(*mem);
    assert(buckets_.size() == other.buckets_.size());
    auto exchange_rate = 0.5;
    for (auto i = 0u; i < buckets_.size(); ++i) {
        auto &bucket_a = buckets_[i];
        auto &bucket_b = other.buckets_[i];
        auto m = std::min(bucket_a.size(), bucket_b.size());
        for (auto j = 0u; j < m; ++j) {
            if (rng.random() < exchange_rate) {
                std::swap(bucket_a[j], bucket_b[j]);
            }
        }
    }
}
