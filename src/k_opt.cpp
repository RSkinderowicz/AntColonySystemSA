#include <algorithm>
#include <cassert>
#include <limits>
#include <random>
#include <iostream>
#include "k_opt.h"

using namespace std;


static void reverse_subroute_inner(std::vector<uint32_t> &route,
                                   std::vector<uint32_t> &pos,
                                   uint32_t beg, uint32_t end) {
    auto i = pos[beg];
    auto j = pos[end];

    assert(i <= j);

    while (i < j) {
        auto c1 = route[i];
        auto c2 = route[j];
        route[i] = c2;
        route[j] = c1;
        pos[c1] = j;
        pos[c2] = i;
        ++i;
        --j;
    }
}


static void reverse_subroute_outer(std::vector<uint32_t> &route,
                                   std::vector<uint32_t> &pos,
                                   uint32_t beg, uint32_t end) {
    const auto N = route.size();
    auto i = pos[beg];
    auto j = pos[end];
    //assert(i >= j);

    auto mid = (N - i) + j + 1;
    mid /= 2;

    for (auto k = 0u; k < mid; ++k) {
        auto c1 = route[i];
        auto c2 = route[j];
        route[i] = c2;
        route[j] = c1;
        pos[c1] = j;
        pos[c2] = i;
        i = (i + 1) % N;
        j = (j > 0 ? j-1 : N - 1);
    }
}

/*
a b c d e f g
....|.....|..
a g c d e f b


a b c d e f g
......|.....|
c b a d e f g
*/


static void reverse_subroute(std::vector<uint32_t> &route,
                             std::vector<uint32_t> &pos,
                             uint32_t beg, uint32_t end) {
    auto i = pos[beg];
    auto j = pos[end];
    if (i <= j) {
        reverse_subroute_inner(route, pos, beg, end);
    } else {
        reverse_subroute_outer(route, pos, beg, end);
    }
}


static void reverse_subroute_flip(std::vector<uint32_t> &route,
                                  std::vector<uint32_t> &pos,
                                  uint32_t beg, uint32_t end) {
    const auto N = route.size();
    auto i = pos[beg];
    auto j = pos[end];
    if (i <= j) {
        if ((j - i) <= N/2 + 1) {
            reverse_subroute_inner(route, pos, beg, end);
        } else {
            auto end_next = route[(j + 1) % N];
            auto beg_prev = route[ i > 0 ? i-1 : N-1 ];
            reverse_subroute(route, pos, end_next, beg_prev);
        }
    } else {
        if ((i - j) <= N/2 + 1) {
            reverse_subroute_outer(route, pos, beg, end);
        } else {
            auto end_next = route[(j + 1) % N];
            auto beg_prev = route[ i > 0 ? i-1 : N-1 ];
            reverse_subroute(route, pos, end_next, beg_prev);
        }
    }
}


Opt2::Opt2(TSPInstance &problem, RNG &rng) :
    rng_(rng),
    problem_(problem) {
    const auto dim = problem.get_dimension();
    const auto nn_size = 32u;
    nn_lists_.reserve(dim);
    for (auto i = 0u; i < dim; ++i) {
        nn_lists_.push_back(problem.get_nearest_neighbors(i, nn_size));
    }
}


void Opt2::apply(std::vector<std::unique_ptr<Ant>> &ants,
                 uint32_t count,
                 const Ant * /*best_so_far*/) {
    for (auto i = 0u; i < count; ++i) {
        try_improve(*ants[i]);
    }
}


bool Opt2::try_improve(Ant &ant) {
    vector<uint32_t> route = ant.get_visited();
    if (apply_2_opt(route)) {
        ant.set_visited(route);
        return true;
    }
    return false;
}


bool Opt2::apply_2_opt(std::vector<uint32_t> &route) {
    auto improvements_count = 0u;
    bool improvement = true;
    const auto N = route.size();

    dont_look_.resize(N+1u);
    fill(begin(dont_look_), end(dont_look_), 0);

    cust_pos_.resize(N+1u);
    for (auto i = 0u; i < route.size(); ++i) {
        cust_pos_[ route[i] ] = i;
    }

    random_order_.resize(N);
    iota(random_order_.begin(), random_order_.end(), 0);
    rng_.shuffle(random_order_.begin(), random_order_.end());

    auto get_distance = [&] (uint32_t a, uint32_t b) { return problem_.get_distance(a, b); };

    while (improvement) {
        improvement = false;

        for (auto cust_i = 0u; cust_i < N; ++cust_i) {
            auto c1 = random_order_[cust_i];
            if (dont_look_[c1] == 1) {
                continue;
            }
            auto pos_c1 = cust_pos_[c1];
            auto c1_next = route[(pos_c1+1) % N];
            double radius = get_distance(c1, c1_next);
            bool improve_node = false;
            uint32_t h1, h2, h3, h4;

            // look through c1's nearest neighbours list
            auto const &nn_list = nn_lists_.at(c1);
            const auto NNCount = nn_list.size();
            for (auto j = 0u; (improve_node == false) && (j < NNCount); ++j) {
                auto c2 = nn_list[j];

                double c1_c2_dist = get_distance(c1, c2);

                if (c1_c2_dist < radius) { // improvement possible
                    auto c2_next = route[ (cust_pos_[c2]+1) % N ];
                    double gain = (c1_c2_dist + get_distance(c1_next, c2_next))
                                - (radius + get_distance(c2, c2_next));
                    if (gain < 0) {
                        h1 = c1;
                        h2 = c1_next;
                        h3 = c2;
                        h4 = c2_next;
                        improve_node = true;
                    }
                } else {
                    // each next c2 will be farther then current c2 so there is no
                    // sense in checking c1_c2_dist < radius condition
                    break;
                }
            }
            if (improve_node == false) {
                auto c1_prev = (cust_pos_[c1] > 0) ? route[ cust_pos_[c1]-1 ] : route[N-1];
                radius = get_distance(c1_prev, c1);

                for (auto j = 0u; (improve_node == false) && (j < NNCount); ++j) {
                    auto c2 = nn_list[j];
                    double c1_c2_dist = get_distance(c1, c2);

                    if (c1_c2_dist < radius) { // improvement possible
                        auto c2_prev = (cust_pos_[c2] > 0) ? route[ cust_pos_[c2]-1 ] : route[N-1];

                        if (c2_prev == c1 || c1_prev == c2) {
                            continue;
                        }
                        double gain = (c1_c2_dist + get_distance(c1_prev, c2_prev))
                                      - (radius + get_distance(c2_prev, c2));
                        if (gain < 0) {
                            h1 = c1_prev;
                            h2 = c1;
                            h3 = c2_prev;
                            h4 = c2;
                            improve_node = true;
                        }
                    } else {
                        // each next c2 will be farther then current c2 so there is no
                        // sense in checking c1_c2_dist < radius condition
                        break;
                    }
                }
            }

            if (improve_node) {
                improvement = true;
                ++improvements_count;

                if (cust_pos_[h3] < cust_pos_[h1]) {
                    swap(h1, h3);
                    swap(h2, h4);
                }
                assert( cust_pos_[h2] < cust_pos_[h3] );

                if (cust_pos_[h3] - cust_pos_[h2] < N/2 + 1) {
                    // reverse inner part from cust_pos[h2] to cust_pos[h3]
                    reverse_subroute_inner(route, cust_pos_, h2, h3);
                } else {
                    // reverse outer part from cust_pos[h4] to cust_pos[h1]
                    if ( cust_pos_[h4] > cust_pos_[h1] ) {
                        reverse_subroute_outer(route, cust_pos_, h4, h1);
                    } else {
                        reverse_subroute_inner(route, cust_pos_, h4, h1);
                    }
                }
                dont_look_[h1] = dont_look_[h2] = dont_look_[h3] = dont_look_[h4] = 0;

            } else { // no improvement
                dont_look_[c1] = 1;
            }
        }
    }
    return improvements_count > 0;
}


//==============================================================================
// Opt3 class
//==============================================================================

Opt3::Opt3(TSPInstance &problem, RNG &rng) :
    rng_(rng),
    problem_(problem) {
    const auto dim = problem.get_dimension();
    const auto nn_size = 40u; // TODO this should be a parameter?
    nn_lists_.reserve(dim);
    for (auto i = 0u; i < dim; ++i) {
        nn_lists_.push_back(problem.get_nearest_neighbors(i, nn_size));
    }
}


void Opt3::apply(std::vector<std::unique_ptr<Ant>> &ants,
                 uint32_t count,
                 const Ant * /*best_so_far*/) {
    for (auto i = 0u; i < count; ++i) {
        try_improve(*ants[i]);
    }
}


bool Opt3::try_improve(Ant &ant) {
    vector<uint32_t> route = ant.get_visited();
    if (apply_3_opt(route)) {
        ant.set_visited(route);
        return true;
    }
    return false;
}


bool Opt3::apply_3_opt(std::vector<uint32_t> &route) {
    //int improvements_count = 0;
    bool improvement = true;
    const auto N = route.size();

    dont_look_.resize(N+1);
    fill(begin(dont_look_), end(dont_look_), 0);

    cust_pos_.resize(N+1);
    for (auto i = 0u; i < route.size(); ++i) {
        cust_pos_[ route[i] ] = i;
    }
    random_order_.reserve(N+1);
    iota(random_order_.begin(), random_order_.end(), 0);
    rng_.shuffle(random_order_.begin(), random_order_.end());

    auto get_distance = [&] (uint32_t a, uint32_t b) { return problem_.get_distance(a, b); };

    /*
       We restric our attention to 6-tuples of cities (a, b, c, d, e, f)

       Each 3-opt move results in removal of the edges (a,b) (c,d) (e,f)

       There are two possibilities for new edges to connect the route together again

       A) edges (a,d) (c,e) (f,b) are created
       B) edges (a,d) (e,b) (c,f) are created

       Ad. A)

       This move needs reversal of subroutes (f..a) and (d..e)

       Ad. B)

       This move needs reversal of three subroutes: (f..a) (b..c) (d..e)

       There is also 2-opt move possible

       edges (a,b) (c,d) are replaced with (a,c) (b,d)

       This move needs reversal of subroute (b..c) or (d..a)

       We need to check only tuples in which d(a,b) > d(a,d) and 

       d(a,b) + d(c,d) > d(a,d) + d(c,e)

       or 

       d(a,b) + d(c,d) > d(a,d) + d(e,b)
       */
    auto total_improvements = 0u;
    while (improvement) {
        improvement = false;

        for (auto cust_i = 0u; cust_i < N; ++cust_i) {
            //a = random_order_[cust_i];
            auto a = cust_i;
            if (dont_look_[a] == 1) {
                continue;
            }
            auto pos_a = cust_pos_[a];
            auto b = route[(pos_a+1) % N];
            auto radius = get_distance(a, b);
            auto dist_a_b = radius;
            bool improve_node = false;
            uint32_t c, d, e, f;

            // look for c in a's nearest neighbours list
            auto const &nn_list = nn_lists_.at(a);
            const auto NNCount = nn_list.size();

            // search for 2-opt move
            for (auto j = 0u; (improve_node == false) && (j < NNCount); ++j) {
                c = nn_list[j];
                auto dist_a_c = get_distance(a, c);
                if (dist_a_c < dist_a_b) { // 2-Opt possible
                    auto pos_c = cust_pos_[c];
                    d = route[(pos_c + 1) % N];
                    auto gain = (dist_a_b + get_distance(c,d)) -
                        (dist_a_c + get_distance(b, d));
                    if (gain > 0) {
                        reverse_subroute_flip(route, cust_pos_, b, c);
                        dont_look_[a] = dont_look_[b] = dont_look_[c] = dont_look_[d] = 0;

                        improve_node = true;
                    }
                } else {
                    // the rest of neighbours are further than c
                    break ;
                }
            }
            auto prev_a = route[ pos_a > 0 ? pos_a-1 : N-1 ];
            auto dist_prev_a_a = get_distance(prev_a, a);
            for (auto j = 0u; (improve_node == false) && (j < NNCount); ++j) {
                c = nn_list[j];
                auto dist_a_c = get_distance(a, c);

                if (dist_a_c < dist_prev_a_a) { // 2-Opt possible
                    auto pos_c = cust_pos_[c];
                    auto prev_c = route[ pos_c > 0 ? pos_c - 1 : N-1 ]; // d is now a predecessor of c

                    if (prev_c == a || prev_a == c)
                        continue ;

                    auto gain = (dist_prev_a_a + get_distance(prev_c, c)) -
                        (dist_a_c + get_distance(prev_a, prev_c));

                    if (gain > 0) {
                        reverse_subroute_flip(route, cust_pos_, c, prev_a);
                        dont_look_[prev_a] = dont_look_[a] = dont_look_[prev_c] = dont_look_[c] = 0;
                        improve_node = true;
                    }
                } else {
                    // the rest of neighbours are further than c
                    break ;
                }
            }

            // search for 3-opt move
            for (auto j = 0u; (improve_node == false) && (j < NNCount); ++j) {
                c = nn_list[j];
                auto pos_c = cust_pos_[c];
                d = route[(pos_c + 1) % N];

                if (d == a) {
                    continue ;
                }

                auto dist_a_d = get_distance(a, d);

                if (dist_a_d < dist_a_b) { // improvement possible -> now find e

                    auto const &nn_list_c = nn_lists_.at(c);
                    // look for e in c's neighbours list
                    for (auto k = 0u; (improve_node == false) && (k < NNCount); ++k) {
                        e = nn_list_c[k];
                        auto pos_e = cust_pos_[e];
                        f = route[ (cust_pos_[e] + 1) % N ];
                        // node e has to lay between nodes c and a, i.e. a..c..e
                        if ( (f == a) ||
                                !( (pos_a < pos_c && (pos_c < pos_e || pos_e < pos_a)) ||
                                    (pos_a > pos_c && (pos_c < pos_e && pos_e < pos_a)) ) )
                            continue;

                        // now check two possibilities
                        auto dist_c_d = get_distance(c, d);
                        auto dist_e_f = get_distance(e, f);

                        // A) edges (a,d) (c,e) (f,b) are created
                        auto gain = dist_a_b + dist_c_d + dist_e_f -
                            (dist_a_d + get_distance(c, e) + get_distance(f, b));
                        if (gain > 0) {
                            // This move needs reversal of subroutes (f..a) and (d..e)
                            reverse_subroute(route, cust_pos_, d, e);
                            reverse_subroute(route, cust_pos_, f, a);
                            dont_look_[a] = dont_look_[b] = dont_look_[c] = 0;
                            dont_look_[d] = dont_look_[e] = dont_look_[f] = 0;
                            improve_node = true;
                            ++total_improvements;
                        }
                        // B) edges (a,d) (e,b) (c,f) are created
                        if (!improve_node) {
                            gain = dist_a_b + dist_c_d + dist_e_f -
                                (dist_a_d + get_distance(e, b) + get_distance(c, f));
                        }
                        if (!improve_node && gain > 0) {
                            // This move needs reversal of three subroutes: (f..a) (b..c) (d..e)
                            reverse_subroute(route, cust_pos_, f, a);
                            reverse_subroute(route, cust_pos_, b, c);
                            reverse_subroute(route, cust_pos_, d, e);
                            dont_look_[a] = dont_look_[b] = dont_look_[c] = 0;
                            dont_look_[d] = dont_look_[e] = dont_look_[f] = 0;
                            improve_node = true;
                            ++total_improvements;
                        }
                        /**
                         * Stop searching if edge (c,e) gets too long
                         */
                        if (!improve_node && (dist_a_b + dist_c_d) < get_distance(c, e))
                            break ;
                    }
                } else {
                    // each next c will be farther then current so there is no
                    // sense in further search
                    break;
                }
            }
            if (improve_node) {
                improvement = true;
            } else {
                auto prev_a = pos_a > 0 ? route[pos_a-1] : route[(pos_a+1) % N];
                if (get_distance(prev_a, a) < dist_a_b) {
                    dont_look_[a] = 1;
                }
            }
        }
    }
    return total_improvements != 0;
}
