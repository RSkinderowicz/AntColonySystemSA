#include <algorithm>
#include <limits>
#include "stop_condition.h"


TimeoutStopCondition::TimeoutStopCondition(double max_seconds, uint64_t steps_per_iter) :
    FixedIterationsStopCondition(std::numeric_limits<uint64_t>::max(), steps_per_iter),
    max_seconds_(std::max(0.0, max_seconds)) {
}


void TimeoutStopCondition::init() {
    FixedIterationsStopCondition::init();
    start_time_ = std::chrono::steady_clock::now();
}


bool TimeoutStopCondition::is_reached() {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(now - start_time_);
    return (elapsed.count() / 1.0e6) >= max_seconds_;
}


double TimeoutStopCondition::get_progress() const {
    const auto now = std::chrono::steady_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(now - start_time_);
    return (elapsed.count() / 1.0e6) / max_seconds_;
}


NeverStopCondition::NeverStopCondition(uint64_t steps_per_iter) :
    FixedIterationsStopCondition(std::numeric_limits<uint64_t>::max(), steps_per_iter)
{}
