#ifndef ACS_NODE_SELECTION_IMPL_H
#define ACS_NODE_SELECTION_IMPL_H

#include "tsp_instance.h"
#include "ant.h"
#include "rng.h"
#include "lru_pheromone_memory.h"
#include "standard_pheromone_memory.h"
#include "aco_distance_heuristic.h"

#include <iostream>


/*
 *  ACS knowledge model consists of both pheromone memory and additional
 *  heuristic information about the problem.
*/
template<typename Problem>
class ACSNodeSelectionImpl {
public:
    using Problem_type = Problem;


    ACSNodeSelectionImpl(const Problem_type &problem,
                         uint32_t cand_list_size) :
        problem_(problem)
    {
        nn_lists_ = problem_.get_nn_lists(cand_list_size);
    }


    void reset() {
    }


    uint32_t get_first_unvisited_nn(const Ant *ant) {
        auto pos = ant->get_position();
        assert(pos < nn_lists_.size());
        for (auto node : nn_lists_[pos]) {
            if (ant->is_available(node)) {
                return node;
            }
        }
        return problem_.get_dimension();
    }


    void get_candidate_nodes(const Ant *ant, std::vector<uint32_t> & out_list) {
        out_list.clear();
        auto pos = ant->get_position();
        for (auto node : nn_lists_[pos]) {
            if (ant->is_available(node)) {
                out_list.push_back(node);
            }
        }
    }


    /**
     * Returns a next node for the ant selected based on the pseudo-random
     * proportional selection.
     */
    uint32_t select_next_node(StandardPheromoneMemory &ph_mem, Ant *ant, double q0, RNG &rng) {
        auto &cand = next_node_candidates_;
        get_candidate_nodes(ant, cand);

        const auto start_node = ant->get_position();
        uint32_t next_node = start_node;

        if (!cand.empty()) {
            auto q = rng.random();

            if (q < q0) { // Greedy selection of the node with a maximum product
                          // of heuristic knowledge and pheromone
                double max_prod = 0.0;
                double pher_in_max_prod = 0.0;
                for (auto node : cand) {
                    // We exploit the fact that the heuristic values for nn are non-increasing
                    // and read heuristic value only if necessary
                    const auto pher = ph_mem.get(start_node, node);
                    if (pher > pher_in_max_prod) {
                        const auto prod = pher * problem_.get_heuristic(start_node, node);
                        if (prod > max_prod) {
                            max_prod = prod;
                            next_node = node;
                            pher_in_max_prod = pher;
                        }
                    }
                }
            } else { // Select randomly with the probability proportional to the
                     // product of heuristic knowledge and pheromone
                double total = 0.0;
                for (auto node : cand) {
                    total += ph_mem.get(start_node, node) * problem_.get_heuristic(start_node, node);
                }
                const auto r = rng.random() * total;
                auto acc = 0.0;
                next_node = cand.back();
                for (auto node : cand) {
                    acc += ph_mem.get(start_node, node) * problem_.get_heuristic(start_node, node);
                    if (acc >= r) {
                        next_node = node;
                        break ;
                    }
                }
            }
        } else { // All nearest neighbours of 'curr_node' are already in the route
            next_node = select_node_greedy(ph_mem, ant);
        }
        assert(ant->is_available(next_node));
        return next_node;
    }


    /**
     * Returns a next node for the ant selected based on the pseudo-random
     * proportional selection without the mixed ACS approach with q_0
     */
    uint32_t select_next_node_p(StandardPheromoneMemory &ph_mem, Ant *ant, RNG &rng) {
        auto &cand = next_node_candidates_;
        get_candidate_nodes(ant, cand);

        const auto start_node = ant->get_position();
        uint32_t next_node = start_node;

        if (!cand.empty()) {
            // Select randomly with the probability proportional to the
            // product of heuristic knowledge and pheromone
            double total = 0.0;
            for (auto node : cand) {
                total += ph_mem.get(start_node, node) * problem_.get_heuristic(start_node, node);
            }
            const auto r = rng.random() * total;
            auto acc = 0.0;
            next_node = cand.back();
            for (auto node : cand) {
                acc += ph_mem.get(start_node, node) * problem_.get_heuristic(start_node, node);
                if (acc >= r) {
                    next_node = node;
                    break ;
                }
            }
        } else { // All nearest neighbours of 'curr_node' are already in the route
            //next_node = select_node_greedy(ph_mem, ant);
            const auto &unvisited = ant->get_unvisited();
            double total = 0.0;
            for (auto node : unvisited) {
                if (ant->is_available(node)) {
                    total += ph_mem.get(start_node, node) * problem_.get_heuristic(start_node, node);
                }
            }
            const auto r = rng.random() * total;
            auto acc = 0.0;
            next_node = unvisited.back();
            for (auto node : unvisited) {
                if (ant->is_available(node)) {
                    acc += ph_mem.get(start_node, node) * problem_.get_heuristic(start_node, node);
                    if (acc >= r) {
                        next_node = node;
                        break ;
                    }
                }
            }
        }
        assert(ant->is_available(next_node));
        return next_node;
    }


    /**
     * Return an unvisited node to which an edge with a maximum product of
     * pheromone trail and heuristic knowledge leads.
     */
    uint32_t select_node_greedy(StandardPheromoneMemory &ph_mem, Ant *ant) {
        const auto sentinel = problem_.get_dimension();
        const auto start_node = ant->get_position();
        const auto &unvisited = ant->get_unvisited();
        double max_prod = 0.0;
        auto next_node = sentinel;
        for (auto node : unvisited) {
            const auto prod = ph_mem.get(start_node, node) * problem_.get_heuristic(start_node, node);
            if (prod > max_prod && ant->is_available(node)) {
                max_prod = prod;
                next_node = node;
            }
        }
        return next_node;
    }


    uint32_t select_next_node(LRUPheromoneMemory &ph_mem, Ant *ant, double q0, RNG &rng) {
        const auto curr_node = ant->get_position();
        assert(curr_node < problem_.get_dimension());
        uint32_t next_node = curr_node;
        const auto first_unvisited_nn = get_first_unvisited_nn(ant);

        if (first_unvisited_nn != problem_.get_dimension()) {
            auto q = rng.random();

            if (q < q0) {
                next_node = first_unvisited_nn;
                const auto first_nn_heuristic = problem_.get_heuristic(curr_node, first_unvisited_nn);
                const auto first_nn_trail = ph_mem.get(curr_node, first_unvisited_nn);
                auto max_prod = first_nn_trail * first_nn_heuristic;

                for (auto &el : ph_mem.get_bucket(curr_node)) {
                    const auto node = el.first;
                    if (el.second > first_nn_trail && ant->is_available(node)) {
                        const auto prod = el.second * problem_.get_heuristic(curr_node, node);
                        if (prod > max_prod) {
                            max_prod = prod;
                            next_node = node;
                            //first_nn_trail = std::max(first_nn_trail, prod / first_nn_heuristic);
                        }
                    }
                }
            } else {
                auto &cand = next_node_candidates_;
                get_candidate_nodes(ant, cand);
                double total = 0.0;
                for (auto node : cand) {
                    total += ph_mem.get(curr_node, node) * problem_.get_heuristic(curr_node, node);
                }
                const auto r = rng.random() * total;
                auto acc = 0.0;
                next_node = cand.back();
                for (auto node : cand) {
                    acc += ph_mem.get(curr_node, node) * problem_.get_heuristic(curr_node, node);
                    if (acc >= r) {
                        next_node = node;
                        break ;
                    }
                }
            }
        } else { // All nearest neighbours of 'curr_node' are already in the route
            next_node = select_node_greedy(ph_mem, ant);
        }
        assert(ant->is_available(next_node));
        return next_node;
    }


    uint32_t select_next_node_p(LRUPheromoneMemory &ph_mem, Ant *ant, RNG &rng) {
        std::cout << "select_next_node_p not implemented yet" << std::endl;
        return 0;
    }


    /**
     * Returns an unvisited node to which an edge with a maximum product of
     * pheromone trail and heuristic knowledge leads.
     *
     * This is a version for the selective pheromone memory (SPM).
     */
    template<typename SPM>
    uint32_t select_node_greedy(SPM &ph_mem, Ant *ant) {
        const auto curr_node = ant->get_position();
        auto max_product = 0.0;
        const auto sentinel = problem_.get_dimension();
        auto next_node = sentinel;
        // First check only the edges for which there is a pheromone trail
        // in trails array
        for (auto &el : ph_mem.get_bucket(curr_node)) {
            auto node = el.first;
            if (ant->is_available(node)) {
                auto prod = el.second * problem_.get_heuristic(curr_node, node);
                if (prod > max_product) {
                    max_product = prod;
                    next_node = node;
                }
            }
        }
        auto &nn = nn_lists_.at(curr_node);
        auto threshold = ph_mem.get_initial_pheromone() * problem_.get_heuristic(curr_node, nn.back());
        // We can't get more than threshold value for any of the rest of the
        // edges because all the trails have initial_pheromone pheromone
        // levels and the heuristic value is no larger than the heuristic
        // value of the last nearest neighbor
        if (next_node == curr_node || threshold > max_product) {
            // The rest of edges do have a default pheromone value, i.e.
            // initial_pheromone_. Find one with the heuristic value.
            auto max_weight = 0.0;
            auto &unvisited = ant->get_unvisited();
            auto max_weight_node = unvisited.front();
            for (auto node : unvisited) {
                auto w = problem_.get_heuristic(curr_node, node);
                if (w > max_weight && ant->is_available(node)) {
                    max_weight = w;
                    max_weight_node = node;
                }
            }
            if (max_weight * ph_mem.get_initial_pheromone() > max_product
                    && ant->is_available(max_weight_node)) {
                next_node = max_weight_node;
            }
        }
        return next_node;
    }

private:

    const Problem_type &problem_;
    std::vector<std::vector<uint32_t>> nn_lists_;
    std::vector<uint32_t> next_node_candidates_;
};


#endif /* ifndef ACS_NODE_SELECTION_IMPL_H */
