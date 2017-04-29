#include "rng.h"

RNG::RNG() :
    seed_(rd_()),
    generator_(seed_),
    real_unit_distribution_(0.0, 1.0)
{
}


RNG::RNG(uint32_t seed) :
    seed_(seed),
    generator_(seed),
    real_unit_distribution_(0.0, 1.0)
{
}


/*
 * Returns a random sample taken from the set { start, start+1, ...,
 * end_inclusive }.
 */
std::vector<uint32_t>
RNG::random_sample(uint32_t start, uint32_t end_inclusive, uint32_t sample_size) {
    assert(sample_size <= end_inclusive - start);
    std::vector<uint32_t> sample;
    sample.reserve(sample_size);
    for (auto i = 0u; i < sample_size; ++i) {
        sample.push_back(start++);
    }
    // replace elements with gradually decreasing probability
    for ( ; start <= end_inclusive; ++start) {
        auto j = rand_uint(0, start);  // important: inclusive range
        if (j < sample_size) {
            sample[j] = start;
        }
    }
    return sample;
}
