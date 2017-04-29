#ifndef PROBLEM_MODEL_H
#define PROBLEM_MODEL_H


#include "aco_distance_heuristic.h"
#include <vector>


/**
 * ProblemModel glues together a problem instance with a problem specific
 * heuristic.
 */
template<typename ProblemType>
class ProblemModel {

    using Heuristic = ACODistanceHeuristic<ProblemType>;

public:

    template<typename Parameters>
    ProblemModel(ProblemType &instance, Parameters &params)
        : instance_(instance),
          heuristic_( new Heuristic(instance, params.beta_) )
    {
    }


    uint32_t get_dimension() const { return instance_.get_dimension(); }


    double get_distance(uint32_t u, uint32_t v) const { return instance_.get_distance(u, v); }


    double get_heuristic(uint32_t u, uint32_t v) const {
        return heuristic_->get(u, v);
    }


    double eval_solution(const std::vector<uint32_t> &route) const {
        return instance_.eval_solution(route);
    }


    bool is_solution_valid(const std::vector<uint32_t> &sol) const  {
        return instance_.is_solution_valid(sol);
    }


    /*
    For the given city returns the list of its nn_count closest neighbours.
    */
    std::vector<uint32_t> get_nearest_neighbors(uint32_t city,
                                                uint32_t nn_count) const {
        return instance_.get_nearest_neighbors(city, nn_count);
    }


    const ProblemType &get_problem_instance() const {
        return instance_;
    }


    /**
     * Returns a vector of vectors with nearest neighbours for subsequent nodes.
     */
    std::vector<std::vector<uint32_t>> get_nn_lists(uint32_t cand_list_size) const {
        std::vector<std::vector<uint32_t>> nn_lists;
        // For each city create a list of its 'n' closest neighbours
        const auto n = get_dimension();
        nn_lists.reserve(n);
        for (auto i = 0u; i < n; ++i) {
            nn_lists.push_back( get_nearest_neighbors(i, cand_list_size) );
        }
        return nn_lists;
    }

private:

    const ProblemType &instance_;
    std::unique_ptr<Heuristic> heuristic_;
};


#endif
