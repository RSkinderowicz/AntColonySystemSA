#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include <vector>
#include <memory>
#include "ant.h"
#include "stop_condition.h"


/**
 * Base, abstract class for all the local search heuristics.
 */
class LocalSearch {
public:

    virtual ~LocalSearch() {}


    /*
     * Applies the underlying heuristic to the ants' solutions.
     * Solutions may or may not be improved by the local search.
     *
     * The LS is applied to the first count elements.
     */
    virtual void apply(std::vector<std::unique_ptr<Ant>> &ants,
                       uint32_t count,
                       const Ant *best_so_far = nullptr) = 0;


    /**
     * As the LS may take a long time we should have an ability to stop it if
     * the stop condition is a time-limit one.
     */
    virtual void set_stop_condition(std::shared_ptr<StopCondition> stop_condition) {
        stop_condition_ = stop_condition;
    }


    virtual void log_statistics() = 0;


    /**
     * Should be called before a new problem instance is solved using this
     * object
     */
    virtual void reset() = 0;

protected:


    virtual bool stop_forced() {
        return stop_condition_ != nullptr && stop_condition_->is_reached();
    }

    std::shared_ptr<StopCondition> stop_condition_ = nullptr;
};


#endif
