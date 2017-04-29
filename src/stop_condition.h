#ifndef STOP_CONDITION_H
#define STOP_CONDITION_H

#include <cstdint>
#include <chrono>
#include <cassert>

/*
 * An auxiliary class to implement various algorithm stop conditions, i.e. based
 * on the number of iterations made, time limit, etc.
 */
class StopCondition {
public:
    /**
     * This should be called before the first use of the other methods.
     */
    virtual void init() = 0;

    /**
     * This should be called after each iteration to allow for the condition to
     * update its internal counter.
     *
     * Paramter n can be used to update to iteration counter based on the
     * number of solutions constructed etc.
     */
    virtual void update(uint64_t n = 0) = 0;

    /**
     * This method returns true when the stopping criterion has been reached,
     * i.e. the maximum number of iteration has been performed or the time limit
     * has been exceeded
     */
    virtual bool is_reached() = 0;

    /**
     * Returns current iteration number
     */
    virtual uint64_t get_iteration() const = 0;


    /**
     * Returns a number from [0, 1] indicating the progress of the computations.
     */
    virtual double get_progress() const = 0;

    /**
     * Sets a linked condition. The link condition gets updated every time this
     * condition is updated.
     */
    virtual void set_linked_condition(StopCondition *condition) = 0;
};



/**
 * This class allows to stop the algorithm after a certain, fixed number of
 * iterations were performed. By iteration we understand a completion of
 * steps_per_iter steps, i.e. each iteration consists of 'steps_per_iter' steps.
 * For example, one step could be equal to construction of a single solution.
 */
class FixedIterationsStopCondition : public StopCondition {
public:

    FixedIterationsStopCondition(uint64_t max_iterations,
                                 uint64_t steps_per_iter) :
            steps_count_(0),
            max_iterations_(max_iterations),
            steps_per_iter_(steps_per_iter) {
        assert(steps_per_iter > 0);
    }

    virtual void init() override {
        steps_count_ = 0;
    }

    /*
     * Updates the internal counters. Should be called after some steps have
     * been completed.
     */
    virtual void update(uint64_t steps_completed = 1) override {
        if (get_iteration() < max_iterations_) {
            steps_count_ += steps_completed;
            if (linked_condition_) {
                linked_condition_->update(steps_completed);
            }
        }
    }

    virtual bool is_reached() override {
        return get_iteration() >= max_iterations_;
    }

    virtual uint64_t get_iteration() const override {
        return steps_count_ / steps_per_iter_;
    }


    /**
     * Returns a number from [0, 1] indicating the progress of the computations.
     */
    virtual double get_progress() const override {
        if (linked_condition_ != nullptr) {
            return linked_condition_->get_progress();
        }
        return get_iteration() / double(max_iterations_);
    }


    void set_linked_condition(StopCondition *condition) override {
        linked_condition_ = condition;
    }

private:

    uint64_t steps_count_;
    uint64_t max_iterations_;
    uint64_t steps_per_iter_;
    // Linked condition gets updated every time this condition is updated
    StopCondition *linked_condition_ = nullptr;
};


/*
 * Allows to stop alg. after a fixed time limit was reached.
 */
class TimeoutStopCondition : public FixedIterationsStopCondition {
public:

    TimeoutStopCondition(double max_seconds, uint64_t steps_per_iter);

    virtual void init();

    virtual bool is_reached();

    /**
     * Returns a number from [0, 1] indicating the progress of the computations.
     */
    virtual double get_progress() const override;

private:

    double max_seconds_;
    std::chrono::steady_clock::time_point start_time_;
};


/*
 * Works similarily to FixedIterationsStopCondition but the is_reached method
 * always returns false to allow for the algorithm to run indefinitely.
 */
class NeverStopCondition final : public FixedIterationsStopCondition {
public:

    NeverStopCondition(uint64_t steps_per_iter);

    bool is_reached() override { return false; }

    virtual double get_progress() const override {
        return 0.0;
    }
};

#endif
