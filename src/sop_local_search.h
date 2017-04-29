#ifndef SOPLOCALSEARCH_H
#define SOPLOCALSEARCH_H

#include "local_search.h"
#include "sop_instance.h"
#include "ant.h"
#include "rng.h"
#include <queue>
#include "sa.h"


/**
 * An implementation of the SOP-3-exchange heuristic for the SOP.
 */
class SOPLocalSearch : public LocalSearch {
public:
    SOPLocalSearch(const SOPInstance &problem, RNG &rng);

    /*
     * Applies the underlying heuristic to the ants' solutions.
     * Solutions may or may not be improved by the local search.
     */
    virtual void apply(std::vector<std::unique_ptr<Ant>> &ants,
                       uint32_t count,
                       const Ant *best_so_far = nullptr) override;


    virtual void log_statistics() override;

    /**
     * Should be called before a new problem instance is solved using this
     * object
     */
    virtual void reset() override;

protected:

    /**
     * Tries to improve the given route. Returns true if succeeded, false
     * otherwise.
     */
    virtual bool improve_route(std::vector<uint32_t> &route,
                              const std::vector<uint32_t> *best_so_far);


    /** Looks for an improving move, i.e. a swapping of two adjacent route
     *  segments.
     *
     *  If the move had been found it is applied using apply_move method,
     *  otherwise the route remains unchanged.
     *
     *  Returns true is succeded, false otherwise.
     */
    virtual bool find_move(std::vector<uint32_t> &route, uint32_t begin_h);

    /** Looks for an improving swap of a two adjacent segments. The search starts at
     * begin_h index and goes to the right.
     */
    virtual bool find_move_forward(std::vector<uint32_t> &route, uint32_t begin_h);

    /** Looks for an improving swap of a two adjacent segments. The search starts at
     * begin_h index and goes to the left, i.e. the second of the swapped segments
     * has to preceed the segment which starts at begin_h index.
     */
    virtual bool find_move_backward(std::vector<uint32_t> &route, uint32_t begin_h);

    /** Applies the move to the route. Swaps two segments defined by the h, i, j
     *  indices. End-nodes for the newly created edges are added to the don't push
     *  stack.
     */
    virtual void apply_move(std::vector<uint32_t> &route,
                            uint32_t h, uint32_t i, uint32_t j);

    /** Don't push stack allows to focus the search only on the recently changed
     * fragments of a solution.
     *
     * If best_so_far solution is given, the search focuses only on the nodes in
     * route * which are displaced relative to the best_so_far.
     */
    virtual void init_dont_push_stack(std::vector<uint32_t> &route,
                                      const std::vector<uint32_t> *best_so_far);

    virtual void dont_push_stack_push(uint32_t element);

    virtual uint32_t dont_push_stack_pop();

    virtual bool try_improve(Ant &ant, const Ant *best_so_far);

    virtual bool is_move_accepted(int cost_decrease, int max_cost_decrease);


    const SOPInstance &problem_;
    RNG &rng_;
    std::vector<uint32_t> mark_;
    std::vector<uint32_t> dont_push_stack_;
    std::vector<uint32_t> on_stack_;
    uint32_t h_max_;
    uint32_t h_min_;
    std::vector<uint32_t> changes_per_application_;
};



/*
 * A hybrid of the SOP-3-exchange with the Simulated Annealing algorithm to
 * improve the convergence.
 *
 * The exponential cooling schedule is used in this version.
 */
class SOPLocalSearchSA final : public SOPLocalSearch {
public:

    SOPLocalSearchSA(const SOPInstance &problem,
                     RNG &rng,
                     double sa_cooling_ratio,
                     double sa_initial_accept_prob,
                     double sa_min_temperature_multiplier);

    virtual void log_statistics() override;

    /**
     * Should be called before a new problem instance is solved using this
     * object
     */
    virtual void reset() override;

protected:

    bool try_improve(Ant &ant, const Ant *best_so_far) override;

    /**
     * In the SA version we can allow a worse move to be accepted with a probability
     * inversely proportional to the difference between the current and the proposed
     * solutions.
     *
     * cost_decrease - cost change introduced by the considered move
     * max_cost_decrease - cost change introduced by the best of the currently
     *                     considered moves
     */
    bool is_move_accepted(int cost_decrease, int max_cost_decrease) override;

    /** Don't push stack allows to focus the search only on the recently changed
     * fragments of a solution.
     *
     * If best_so_far solution is given, the search focuses only on the nodes in
     * route * which are displaced relative to the best_so_far.
     */
    virtual void init_dont_push_stack(std::vector<uint32_t> &route,
                                      const std::vector<uint32_t> *best_so_far);

    // Simulated Annealing related parameters grouped for convenience and
    // clarity
    struct SAControl {
        uint64_t sampling_moves_ = 100u;
        uint64_t sampling_better_moves_ = 0u;
        double sampling_total_delta_ = 0.0;
        double init_prob_ = 0.01;
        double cooling_factor_ = 0.9999;
        double T_ = std::numeric_limits<double>::min();
    };

    SAControl sa_ctrl_;
    const SAControl initial_sa_;
    uint64_t worse_moves_accepted_ = 0u;
    std::vector<double> init_deltas_;
    std::unique_ptr<SimulatedAnnealing> sa_ = nullptr;
};


#endif // SOPLOCALSEARCH_H
