/**

The source code is licensed under the [MIT
License](http://opensource.org/licenses/MIT):

Copyright 2017 Rafał Skinderowicz (rafal.skinderowicz@us.edu.pl)

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/


#ifndef SA_H
#define SA_H

#include "rng.h"
#include <iostream>
#include <vector>


/** Simple impl. of the Simulated Annealing algorithm with the exponential
 *  cooling schedule.
 */
class SimulatedAnnealing {
public:
    struct SAControl {
        double T_ = std::numeric_limits<double>::min();
        uint64_t improving_moves_count_ = 0u;
        double initial_temperature_ = std::numeric_limits<double>::min();
        double temp_init_delta_sum_ = 0;
        double min_temperature_ = 0;
        uint64_t worse_accepted_ = 0;
        uint64_t worse_rejected_ = 0;
    };


    SimulatedAnnealing(RNG &rng,
             double sa_cooling_ratio,
             double sa_initial_accept_prob,
             double sa_min_temperature_multiplier) :
        rng_(rng),
        cooling_ratio_(sa_cooling_ratio),
        initial_accept_prob_(sa_initial_accept_prob),
        min_temperature_multiplier_(sa_min_temperature_multiplier)
    {
    }


    void reset() {
        ctrl_ = initial_ctrl_;
    }


    bool is_active() const { return ctrl_.T_ != std::numeric_limits<double>::min(); }

    /**
     * Initializes the temperature based on values of sample solutions.
     * The sample should be large (> 100). According to the central limit
     * theorem the mean of a large sample of independent random variables is
     * approximately normally distributed. The mean and the stdev are used to
     * calc. the max. value of the delta with 99.7% confidence which is
     * subsequently used to calc. the initial temp. so that the probability of
     * accepting the max delta equals initial_accept_prob_.
     */
    void init_temperature(std::vector<double> sol_values) {
        assert(!sol_values.empty());
        // Draw a relatively large random sample approximating the distribution of
        // differences between values of randomly generated solutions which
        // values are passed in sol_values.
        std::vector<double> deltas;
        const auto sample_size = 10000u;
        const auto N = sol_values.size();
        while (deltas.size() < sample_size) {
            auto fst = sol_values[rng_.rand_uint(0, N-1)];
            auto snd = sol_values[rng_.rand_uint(0, N-1)];
            deltas.push_back(fabs(fst - snd));
        }
        init_temperature_from_deltas(deltas);
    }


    void init_temperature_from_deltas(const std::vector<double> deltas) {
        auto &v = deltas;
        const auto sum = std::accumulate(std::begin(v), std::end(v), 0.0);
        const auto mean = sum / v.size();
        double accum = 0.0;
        std::for_each (std::begin(v), std::end(v), [&](const double d) {
            accum += (d - mean) * (d - mean);
        });
        const auto stdev = sqrt(accum / (v.size()-1));

        //std::cout << "mean uphill move: " << mean << "\tstdev: " << stdev << std::endl;
        ctrl_.initial_temperature_ = -(mean + 3 * stdev) / log(initial_accept_prob_);
        ctrl_.T_ = ctrl_.initial_temperature_;
        ctrl_.min_temperature_ = min_temperature_multiplier_ * ctrl_.initial_temperature_;
    }


    void reset_temperature_to_initial() {
        ctrl_.T_ = ctrl_.initial_temperature_;
    }


    void update_temperature() {
        if (ctrl_.T_ > 0) {
            // Do not allow the temperature to fall too much
            ctrl_.T_ = std::max(ctrl_.min_temperature_,
                              ctrl_.T_ * cooling_ratio_);
        }
    }


    /** Returns true if a solution with the candidate value should be accepted
     *  given that a current solution has the current value.
     *
     *  If the candidate < current true is returned,
     *  otherwise the candidate solution may be accepted with some prob.
     *  depending on the temperature.
     */
    bool try_accept_solution(double candidate, double current) {
        return try_accept_solution_by_delta(candidate - current);
    }


    /** Returns true if the solution for which the difference between its value
     *  and the current solution is delta should be accepted.
     *  A better solution (delta < 0) is always accepted,
     *  a worse one according to the Metropolis criterion,
     *  i.e. r < exp(-delta/T)
     */
    bool try_accept_solution_by_delta(double delta) {
        bool accept = false;
        auto &T = ctrl_.T_;

        if (delta < 0) { // Always accept a better solution
            accept = true;
            ++ctrl_.improving_moves_count_;
        } else if (T > 0) { // Annealing phase, worse (or not better) solution
            const auto prob = calc_worse_solution_accept_prob(delta);
            if (rng_.random() < prob) {
                accept = true;
                ctrl_.worse_accepted_++;
            } else {
                ctrl_.worse_rejected_++;
            }
        }
        return accept;
    }


    double calc_worse_solution_accept_prob(double delta) {
        assert( delta >= 0 );
        return (ctrl_.T_ > 0) ?  exp(-delta / ctrl_.T_) : 0.0;
    }


    SAControl get_control_parameters() const {
        return ctrl_;
    }

private:

    RNG &rng_;
    double cooling_ratio_ = 0.999;
    double initial_accept_prob_ = 0.9;
    double min_temperature_multiplier_ = 1e-4;
    const SAControl initial_ctrl_;
    SAControl ctrl_;
};

#endif
