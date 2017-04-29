#include <iostream>
#include "sop_local_search.h"
#include "log.h"

using namespace std;


/**
 * Swaps two adjacent subsequences at indices
 *  [start .. mid-1] [mid .. end-1] of the
 * vector v resulting in a new subsequence with order of elements
 *  [mid .. end-1] [start .. mid-1]
 * The end_1 and end_2 are inclusive indices.
 */
void swap_adj_subsequences(std::vector<uint32_t> &v, uint32_t start, uint32_t mid, uint32_t end) {
    assert(end <= v.size());
    assert(start < mid && mid < end);

    std::rotate( std::begin( v ) + start,
                 std::begin( v ) + mid,
                 std::begin( v ) + end );
}


SOPLocalSearch::SOPLocalSearch(const SOPInstance &problem, RNG &rng)
    : problem_(problem),
      rng_(rng)
{}


void SOPLocalSearch::init_dont_push_stack(vector<uint32_t> &route, const vector<uint32_t> *best_so_far) {
    const auto n = problem_.get_dimension();
    dont_push_stack_.clear();
    on_stack_.resize(n);

    if (best_so_far == nullptr || best_so_far->empty()) {
        dont_push_stack_.resize(n);
        for (auto i = 0u; i < n; ++i) {
            dont_push_stack_[i] = i;
            on_stack_[i] = 1;
        }
    } else {
        // Focus search only on new edges
        const auto &best = *best_so_far;
        vector<uint32_t> indices(n);
        auto i = 0u;
        for (auto node : best) {
            indices.at(node) = i++;
        }
        for (auto j = 1u; j < n-1; ++j) {
            const auto prev = route.at(j - 1);
            const auto node = route.at(j);
            const auto next = route.at(j + 1);

            const auto idx = indices.at( node );
            const auto best_prev = ( idx > 0 ) ? best.at( idx - 1 ) : 0;
            const auto best_next = ( idx+1 < n ) ? best.at( idx + 1 ) : 0;

            // Is node connected to different neighbours?
            if (prev != best_prev || next != best_next) {
                dont_push_stack_push(prev);
            }
        }
    }
    // This randomly changes the order of the elements in the don't push stack
    // It seems to help with some instances.
    rng_.shuffle(begin(dont_push_stack_), end(dont_push_stack_));
}


void SOPLocalSearch::dont_push_stack_push(uint32_t element) {
    assert(element < problem_.get_dimension());
    if (on_stack_[element] == 0) {
        on_stack_[element] = 1;
        dont_push_stack_.push_back(element);
    }
}


uint32_t SOPLocalSearch::dont_push_stack_pop() {
    assert( !dont_push_stack_.empty() );

    auto last = dont_push_stack_.back();
    dont_push_stack_.pop_back();
    on_stack_[last] = 0;

    return last;
}


bool SOPLocalSearch::improve_route(vector<uint32_t> &route, const vector<uint32_t> *best_so_far) {
    mark_.resize(route.size());
    init_dont_push_stack(route, best_so_far);
    h_min_ = 1;
    h_max_ = route.size()-2;

    bool success = false;
    auto moves_made = 0u;
    while (!dont_push_stack_.empty() && !stop_forced()) {
        const auto h = dont_push_stack_pop();
        const auto h_idx = find(route.begin(), route.end(), h) - route.begin();
        const auto move_found = find_move(route, h_idx);
        success |= move_found;
        moves_made += (uint64_t)move_found;
    }
    changes_per_application_.push_back(moves_made);
    return success;
}


bool SOPLocalSearch::find_move(vector<uint32_t> &route, uint32_t begin_h) {
    return find_move_forward(route, begin_h) || find_move_backward(route, begin_h);
}


bool SOPLocalSearch::find_move_forward(vector<uint32_t> &route, uint32_t begin_h) {
    const auto n = route.size();
    assert(begin_h < n);

    const auto Sentinel = numeric_limits<uint32_t>::max();

    if (begin_h > n - 2 || begin_h >= h_max_) {
        return false;
    }

    fill(mark_.begin(), mark_.end(), Sentinel);
    const auto h = begin_h;
    const auto _h = route[h];
    const auto _h_1 = route[h + 1];
    assert(_h_1 < n);
    const auto dist_h_h1 = problem_.get_distance(_h, _h_1);

    for (auto i = h + 1; i < n-1; ++i) {
        const auto _i_1 = route[i + 1];
        if (mark_[_i_1] == h) {
            break ;
        }
        const auto _i = route[i];
        assert(_i_1 < n);
        const auto dist_i_i1 = problem_.get_distance(_i, _i_1);
        const auto dist_h_i1 = problem_.get_distance(_h, _i_1);

        // Labelling procedure
        for (auto &node : problem_.get_outgoing_edges(_i)) {
            mark_[node] = h;
        }

        auto best_j = Sentinel;
        const auto h_i_move = dist_h_h1 + dist_i_i1 - dist_h_i1;
        int max_cost_decrease = 0;

        for (auto j = i + 1, _j = route[j]; mark_[_j] != h; ++j) {
            const auto _j_1 = route[j + 1];
            assert(_j_1 < n);
            const auto cost_decrease = h_i_move
                                     + problem_.get_distance(_j, _j_1)
                                     - problem_.get_distance(_j, _h_1)
                                     - problem_.get_distance(_i, _j_1);

            if (is_move_accepted(cost_decrease, max_cost_decrease)) {
                max_cost_decrease = cost_decrease;
                best_j = j;
            }
            _j = _j_1;
        }
        if (best_j != Sentinel ) {
            apply_move(route, h+1, i+1, best_j+1);
            return true;
        }
    }

    return false;
}


bool SOPLocalSearch::find_move_backward(vector<uint32_t> &route, uint32_t begin_h) {
    assert(begin_h < route.size());

    if (begin_h < 1 || begin_h <= h_min_) {
        return false;
    }

    const auto Sentinel = numeric_limits<uint32_t>::max();
    fill(mark_.begin(), mark_.end(), Sentinel);

    const auto h = begin_h;
    const auto _h = route[h];
    const auto _h_1 = route[h - 1];
    const auto dist_h1_h = problem_.get_distance(_h_1, _h);

    for (auto i = h - 1; i > 0; --i) {
        const auto _i_1 = route[i - 1];
        if (mark_[_i_1] == h) {
            break ;
        }

        const auto _i = route[i];
        assert(_i < route.size());
        const auto dist_i1_i = problem_.get_distance(_i_1, _i);
        const auto dist_i1_h = problem_.get_distance(_i_1, _h);

        // Labelling procedure
        for (auto &node : problem_.get_incoming_edges(_i)) {
            mark_[node] = h;
        }

        auto best_j = Sentinel;
        const auto h_i_move = dist_h1_h + dist_i1_i - dist_i1_h;
        int max_cost_decrease = 0;

        for (auto j = i - 1, _j = route[j]; mark_[_j] != h; --j) {
            const auto _j_1 = route[j - 1];

            assert(_j < route.size());
            const auto cost_decrease = h_i_move
                                     + problem_.get_distance(_j_1, _j)
                                     - problem_.get_distance(_h_1, _j)
                                     - problem_.get_distance(_j_1, _i);

            if (is_move_accepted(cost_decrease, max_cost_decrease)) {
                max_cost_decrease = cost_decrease;
                best_j = j;
            }
            _j = _j_1;
        }
        if (best_j != Sentinel) {
            //apply_move(route, best_j, i, h);

            dont_push_stack_push(route[h]);
            dont_push_stack_push(route[h-1]);
            dont_push_stack_push(route[i]);
            dont_push_stack_push(route[i-1]);
            dont_push_stack_push(route[best_j]);
            dont_push_stack_push(route[best_j-1]);

            swap_adj_subsequences(route, best_j, i, h);

            return true;
        }
    }
    return false;
}


/*
 * Applies the move to the route. Swaps two segments defined by the h, i, j
 * indices.  End-nodes for the newly created edges are added to the don't push
 * stack.
 */
void SOPLocalSearch::apply_move(std::vector<uint32_t> &route,
                    uint32_t h, uint32_t i, uint32_t j) {

    dont_push_stack_push(route[j-1]);
    dont_push_stack_push(route[j]);
    dont_push_stack_push(route[i-1]);
    dont_push_stack_push(route[i]);
    dont_push_stack_push(route[h-1]);
    dont_push_stack_push(route[h]);

    swap_adj_subsequences(route, h, i, j);
}


/*
 * Applies the underlying heuristic to the ants' solutions.
 * Solutions may or may not be improved by the local search.
 */
void SOPLocalSearch::apply(std::vector<std::unique_ptr<Ant>> &ants,
                           uint32_t count,
                           const Ant *best_so_far) {
    for (auto i = 0u; i < count; ++i) {
        try_improve(*ants[i], best_so_far);
    }
}



bool SOPLocalSearch::try_improve(Ant &ant, const Ant *best_so_far) {
    auto route = ant.get_visited();
    const auto success = (best_so_far != nullptr)
                         ? improve_route(route, &best_so_far->get_visited())
                         : improve_route(route, nullptr);
    if (success) {
        ant.set_visited(route);
    }
    return success;
}


bool SOPLocalSearch::is_move_accepted(int cost_decrease, int max_cost_decrease) {
    return cost_decrease > max_cost_decrease;
}


void SOPLocalSearch::log_statistics() {
    const auto sum = std::accumulate(std::begin(changes_per_application_),
                                     std::end(changes_per_application_), 0u);

    log_add("SOPLocalSearch_total_changes", sum);
    const auto mean = (sum == 0) ? 0.0
                                 : sum / (double)changes_per_application_.size();
    log_add("SOPLocalSearch_mean_changes", mean);
}


/**
 * Should be called before a new problem instance is solved using this
 * object
 */
void SOPLocalSearch::reset() {
    changes_per_application_.clear();
}


//==============================================================================
// SOPLocalSearchSA
//==============================================================================


SOPLocalSearchSA::SOPLocalSearchSA(const SOPInstance &problem, RNG &rng,
             double sa_cooling_ratio,
             double sa_initial_accept_prob,
             double sa_min_temperature_multiplier)
    : SOPLocalSearch(problem, rng)
{
    sa_.reset( new SimulatedAnnealing(
            rng,
            sa_cooling_ratio,
            sa_initial_accept_prob,
            sa_min_temperature_multiplier));
}


void SOPLocalSearchSA::init_dont_push_stack(vector<uint32_t> &route,
                                            const vector<uint32_t> *best_so_far) {
    const auto n = problem_.get_dimension();
    dont_push_stack_.clear();
    on_stack_.resize(n);

    if (best_so_far == nullptr || best_so_far->empty()) {
        dont_push_stack_.resize(n);
        for (auto i = 0u; i < n; ++i) {
            dont_push_stack_[i] = i;
            on_stack_[i] = 1;
        }
    } else {
        // Focus search only on new edges
        const auto &best = *best_so_far;
        vector<uint32_t> indices(n);
        auto i = 0u;
        for (auto node : best) {
            indices.at(node) = i++;
        }
        for (auto j = 1u; j < n-1; ++j) {
            const auto prev = route.at(j - 1);
            const auto node = route.at(j);
            const auto next = route.at(j + 1);

            const auto idx = indices.at( node );
            const auto best_prev = ( idx > 0 ) ? best.at( idx - 1 ) : 0;
            const auto best_next = ( idx+1 < n ) ? best.at( idx + 1 ) : 0;

            // Is node connected to different neighbours?
            if (prev != best_prev || next != best_next) {
                dont_push_stack_push(prev);
            }
        }
        if (dont_push_stack_.empty()) {
            for (auto i = 0u; i < n; ++i) {
                if (rng_.random() < 0.1) {
                    dont_push_stack_push(i);
                }
            }
        }
    }
    // This randomly changes the order of the elements in the don't push stack
    // It seems to help with some instances.
    rng_.shuffle(begin(dont_push_stack_), end(dont_push_stack_));
}



/**
 * In the SA version we can allow a worse move to be accepted with a probability
 * inversely proportional to the difference between the current and the proposed
 * solutions.
 *
 * cost_decrease - cost change introduced by the considered move
 * max_cost_decrease - cost change introduced by the best of the currently
 *                     considered moves
 */
bool SOPLocalSearchSA::is_move_accepted(int cost_decrease, int max_cost_decrease) {
    bool accept = false;
    const auto delta = cost_decrease - max_cost_decrease;

    if (delta < 0) {
        if (init_deltas_.size() < 10000) {
            init_deltas_.push_back(-delta);
        } else if (!sa_->is_active()) {
            sa_->init_temperature_from_deltas(init_deltas_);
        }
    }

    if (delta > 0) {
        accept = true; // Always accept an uphill (better) move
        ++sa_ctrl_.sampling_better_moves_;
    } else if ( delta == 0 ) {
        // If delta = 0 the move does not change solution's value, but the
        // change in the configuration could be useful On the other hand, we do
        // not want accept all such changes to not slow the algorithm.
        accept = (rng_.random() < 0.1);
    } else if (sa_->is_active()) {
        accept = sa_->try_accept_solution_by_delta( -delta );
        sa_->update_temperature();
        worse_moves_accepted_ += (uint64_t)accept;
    }
    return accept;
}


bool SOPLocalSearchSA::try_improve(Ant &ant, const Ant *best_so_far) {
    sa_ctrl_ = initial_sa_; // reset SA parameters
    sa_->reset();
    //init_deltas_.clear();
    return SOPLocalSearch::try_improve(ant, best_so_far);
}


void SOPLocalSearchSA::log_statistics() {
    SOPLocalSearch::log_statistics();
    log_add("SOPLocalSearchSA_worse_moves_accepted", worse_moves_accepted_);
}



/**
 * Should be called before a new problem instance is solved using this
 * object
 */
void SOPLocalSearchSA::reset() {
    SOPLocalSearch::reset();
    worse_moves_accepted_ = 0;
}


/**
 * This is a single-method version of the SOP-3-exchange heuristic
 */
uint32_t sop3_exchange(std::vector<uint32_t> &tour, const SOPInstance &instance) {
    const auto n = tour.size();

    uint32_t bh, bi, bj;
    vector<uint32_t> mark(tour.size(), 0);
    vector<uint32_t> stack;
    vector<bool> used(n);
    vector<uint32_t> tmp;
    uint32_t count_h = 1;

    //init stack
    stack.clear();
    stack.resize(n-1);
    iota(begin(stack), end(stack), 0u);

    fill_n(begin(used), n-1, true);

    auto improvements = 0u;
    bool gainExists = true;

    while (gainExists) {
        gainExists = false;
        //randomizeStack(stack);

        while (!gainExists && !stack.empty() ) {
            const auto h = stack.back();
            stack.pop_back();
            assert(h < n);

            ++count_h;
            auto bestGain = 0;

            used.at(h) = false;
            for (auto i = h + 1; i+1 < n; i++) {
                auto bestj = n;
                auto cannotAdd = false;

                for (auto node : instance.get_outgoing_edges(tour[i])) {
                    mark[node] = count_h;
                }
                for (auto j = i + 1u; j+1 < n && !cannotAdd; j++) {
                    auto res = 0;

                    //feasible
                    if (mark[tour[j]] == count_h) {
                        cannotAdd = true;
                    } else {
                        //calculate gain
                        res+= instance.get_distance(tour.at(h), tour.at(h+1));
                        res+= instance.get_distance(tour.at(i), tour.at(i+1));
                        res+= instance.get_distance(tour.at(j), tour.at(j+1));

                        res-= instance.get_distance(tour.at(h), tour.at(i+1));
                        res-= instance.get_distance(tour.at(j), tour.at(h+1));
                        res-= instance.get_distance(tour.at(i), tour.at(j+1));
                        if (res > bestGain) {
                            bestj = j;
                            bestGain = res;
                            bh = h;
                            bi = i;
                            bj = j;
                        }
                    }
                }
                if (bestj != n) break;
            }

            //Apply the exchange
            if (bestGain > 0) {
                ++improvements;

                swap_adj_subsequences(tour, bh+1, bi+1, bj+1);

                gainExists = true;
                assert(bi+1 < n);
                assert(bh+1 < n);
                assert(bj+1 < n);
                if (!used.at(bj+1)) { stack.push_back(bj+1); used.at(bj+1) = true; }
                if (!used.at(bj)) { stack.push_back(bj); used.at(bj) = true; }
                if (!used.at(bi+1)) { stack.push_back(bi+1); used.at(bi+1) = true; }
                if (!used.at(bi)) { stack.push_back(bi); used.at(bi) = true; }
                if (!used.at(bh+1)) { stack.push_back(bh+1); used.at(bh+1) = true; }
                if (!used.at(bh)) { stack.push_back(bh); used.at(bh) = true; }

                continue;
            }

            ++count_h;

            //Backward search
            bestGain = 0;
            if (h + 1 == n || h < 2) {
                continue;
            }
            for (auto i = h - 1; i >= 1; --i) {
                auto bestj = n;
                auto cannotAdd = false;

                for (auto node : instance.get_incoming_edges(tour.at(i+1))) {
                    mark.at(node) = count_h;
                }
                //Special case for backward direction
                //if (mark.at(tour.at(i)) == count_h) {
                    //cannotAdd = true;
                //}
                for (auto j = i - 1, o = i; o >= 1 && !cannotAdd; --j, --o) {
                    auto res = 0;

                    //feasible
                    if (mark.at(tour.at(j+1))==count_h) {
                        cannotAdd = true;
                    } else {
                        //calculate gain
                        assert(h+1 < n);
                        // [...j] [j+1...i] [i+1..h] [h+1 ...]
                        // [...j] [i+1..h] [j+1...i] [h+1 ...]
                        res+= instance.get_distance(tour.at(h), tour.at(h+1));
                        res+= instance.get_distance(tour.at(i), tour.at(i+1));
                        res+= instance.get_distance(tour.at(j), tour.at(j+1));

                        res-= instance.get_distance(tour.at(i), tour.at(h+1));
                        res-= instance.get_distance(tour.at(j), tour.at(i+1));
                        res-= instance.get_distance(tour.at(h), tour.at(j+1));

                        if (res > bestGain) {
                            bestj = j;
                            bestGain = res;
                            bh = h;
                            bi = i;
                            bj = j;
                        }
                    }
                }
                if (bestj != n) break;
            }

            //Apply the exchange
            if (bestGain != 0) {
                ++improvements;

                swap_adj_subsequences(tour, bj+1, bi+1, bh+1);

                gainExists = true;

                assert(bh+1 < n);
                assert(bi+1 < n);
                assert(bj+1 < n);
                if (!used.at(bh+1)) { stack.push_back(bh+1); used.at(bh+1) = true; }
                if (!used.at(bh)) { stack.push_back(bh); used.at(bh) = true; }
                if (!used.at(bi+1)) { stack.push_back(bi+1); used.at(bi+1) = true; }
                if (!used.at(bi)) { stack.push_back(bi); used.at(bi) = true; }
                if (!used.at(bj+1)) { stack.push_back(bj+1); used.at(bj+1) = true; }
                if (!used.at(bj)) { stack.push_back(bj); used.at(bj) = true; }

            }
        }
    }
    return improvements;
}
