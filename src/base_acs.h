#ifndef BASE_ACS_H
#define BASE_ACS_H

#include "rng.h"
#include "ant.h"
#include "ant_sop.h"
#include "local_search.h"
#include "stop_condition.h"



/** Abstract base class for ACS-related algorithms.
 */
class BaseACS {
public:

    BaseACS(RNG &rng)
        : rng_(rng) {}


    /** Initializes the computations */
    virtual void run_begin(std::shared_ptr<StopCondition> stop_condition) = 0;

    /** Executes only a single iteration of the alg. */
    virtual void run_next_iteration() = 0;

    /** Executes a complete algorithm until a stop condition is met. */
    virtual void run() = 0;

    RNG &get_rng() { return rng_; }


    void set_local_search(std::shared_ptr<LocalSearch> ls) {
        local_search_ = ls;
    }

    std::shared_ptr<LocalSearch> get_local_search() {
        return local_search_;
    }

protected:


    /** A helper method which returns a pointer to the ant with the smallest
     *  solution value */
    Ant* get_best_ant(std::vector<std::unique_ptr<Ant>> &ants) {
        assert(!ants.empty());
        auto &best_ant = *min_element(ants.begin(), ants.end(),
                          [=] (std::unique_ptr<Ant> &a, std::unique_ptr<Ant> &b) {
                                  return a->get_value() < b->get_value();
                          });
        return best_ant.get();
    }


    RNG &rng_;
    std::shared_ptr<LocalSearch> local_search_ = nullptr;
};



Ant *create_ant(const ProblemModel<TSPInstance> &problem) {
    return new Ant(problem.get_dimension());
}


Ant *create_ant(const ProblemModel<SOPInstance> &problem) {
    return new AntSOP(problem.get_problem_instance());
}

#endif /* BASE_ACS_H */
