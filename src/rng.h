#ifndef RNG_H
#define RNG_H

#include <random>
#include <algorithm>
#include <cassert>


/*
 * RNG - random number generator - a wrapper class to simplify random number
 * generation.
 */
class RNG {
public:
    RNG();

    RNG(uint32_t seed);

    /* Returns a random number in the range [0.0, 1.0] */
    double random() { return real_unit_distribution_(generator_); }


    uint32_t rand_uint(uint32_t start, uint32_t end_inclusive) {
        auto dist = std::uniform_int_distribution<uint32_t>(start, end_inclusive);
        return dist(generator_);
    }


    uint32_t get_seed() const { return seed_; }


    template<class RandomIt>
    void shuffle(RandomIt first, RandomIt last) {
        std::shuffle(first, last, generator_);
    }


    /*
     * Returns a random sample taken from the set { start, start+1, ...,
     * end_inclusive }.
     */
    std::vector<uint32_t> random_sample(uint32_t start, uint32_t end_inclusive,
                                        uint32_t sample_size);

private:
    std::random_device rd_; // We need this to initialize our generator
    const uint32_t seed_;
    std::mt19937 generator_;
    std::uniform_real_distribution<double> real_unit_distribution_;
};

#endif /* ifndef RNG_H */
