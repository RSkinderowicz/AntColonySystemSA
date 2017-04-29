#ifndef K_OPT_H
#define K_OPT_H

#include "tsp_instance.h"
#include "local_search.h"
#include "rng.h"


/*
 * A functor for performing a modified version of 2-opt local search algorithm
 * with dont-look-bits heuristic to cut down the runtime.
 */
class Opt2 : public LocalSearch {
public:

    Opt2(TSPInstance &problem, RNG &rng);


    /*
     * Applies the underlying heuristic to the ants' solutions.
     * Solutions may or may not be improved by the local search.
     */
    void apply(std::vector<std::unique_ptr<Ant>> &ants,
               uint32_t count,
               const Ant *best_so_far=nullptr) override;


    void log_statistics() override {}


    /**
     * Should be called before a new problem instance is solved using this
     * object
     */
    void reset() override {};

private:

    bool try_improve(Ant &ant);


    /**
     * Tries to improve the solution, returns true if succedded.
     **/
    bool apply_2_opt(std::vector<uint32_t> &route);


    RNG &rng_;
    TSPInstance &problem_;
    std::vector<std::vector<uint32_t>> nn_lists_;
    std::vector<uint8_t> dont_look_;
    std::vector<uint32_t> cust_pos_;
    std::vector<uint32_t> random_order_;
};


/*
 * A functor for performing a modified version of 3-opt local search algorithm
 * with dont-look-bits heuristic to cut down the runtime.
 */
class Opt3 : public LocalSearch {
public:

    Opt3(TSPInstance &problem, RNG &rng);


    /*
     * Applies the underlying heuristic to the ants' solutions.
     * Solutions may or may not be improved by the local search.
     */
    void apply(std::vector<std::unique_ptr<Ant>> &ants,
               uint32_t count,
               const Ant *best_so_far=nullptr) override;

    void log_statistics() override {}

    /**
     * Should be called before a new problem instance is solved using this
     * object
     */
    void reset() override {};

private:

    bool try_improve(Ant &ant);


    /**
     * Tries to improve the solution, returns true if succedded.
     **/
    bool apply_3_opt(std::vector<uint32_t> &route);


    RNG &rng_;
    TSPInstance &problem_;
    std::vector<std::vector<uint32_t>> nn_lists_;
    std::vector<uint8_t> dont_look_;
    std::vector<uint32_t> cust_pos_;
    std::vector<uint32_t> random_order_;
};

#endif
