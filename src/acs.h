#ifndef ACS_H
#define ACS_H

#include "acs_parameters.h"
#include "tsp_instance.h"
#include "sop_instance.h"
#include "progress_visitor.h"
#include <iostream>
#include "problem_model.h"
#include "base_acs.h"


/*
 * Implementation of (standard) the Ant Colony System algorithm for the TSP
 * problem. In this version the pheromone memory is stored separately what
 * allows for a different implementations of the memory.
 */
template<typename NodeSelectionImpl, typename PheromoneMemory>
class ACS : public BaseACS {
public:
    using Problem_type = typename NodeSelectionImpl::Problem_type;
    using Parameters = ACSParameters;

    ACS(Parameters param,
            std::shared_ptr<NodeSelectionImpl> node_sel_impl,
            std::shared_ptr<PheromoneMemory> ph_mem,
            Problem_type &problem,
            RNG &rng)
        : BaseACS(rng),
          param_(param),
          node_sel_impl_(node_sel_impl),
          pheromone_memory_(ph_mem),
          problem_(problem),
          global_best_ant_(problem.get_dimension())
    {
    }


    void run_begin(std::shared_ptr<StopCondition> stop_condition) override {
        stop_condition_ = stop_condition;
        stop_condition_->init();

        node_sel_impl_->reset();
        pheromone_memory_->reset();

        global_best_ant_.reset();
        global_best_found_iter_ = 0;

        if (progress_visitor_) {
            progress_visitor_->started(*this);
        }
    }


    void run_next_iteration() override {
        build_ants_solutions();
        select_global_best_ant();
        pheromone_memory_->global_update(global_best_ant_.get_visited(),
                                         global_best_ant_.get_value(),
                                         param_.rho_);

        if (progress_visitor_) {
            progress_visitor_->iteration_done(*this);
        }
        stop_condition_->update(param_.ants_count_);
    }


    void run() override {
        while (!stop_condition_->is_reached()) {
            run_next_iteration();
        }
        if (progress_visitor_) {
            progress_visitor_->finished(*this);
        }
    }


    const Ant &get_global_best_ant() const { return global_best_ant_; }

    /*
     * Returns the number of iteration during which the current global best
     * solution was found.
     */
    uint64_t get_global_best_found_iter() const { return global_best_found_iter_; }


    uint64_t get_current_iter() const {
        if (!stop_condition_) {
            throw std::runtime_error("ACS::get_current_iter: No stop condition given");
        }
        return stop_condition_->get_iteration();
    }


    void set_progress_visitor(std::shared_ptr<ProgressVisitor<ACS<NodeSelectionImpl, PheromoneMemory>>> visitor) {
        progress_visitor_ = visitor;
    }


    uint32_t get_ants_count() const { return param_.ants_count_; }


private:

    void build_ants_solutions() {
        const auto dimension = problem_.get_dimension();

        ants_.clear();
        for (auto i = 0u; i < param_.ants_count_; ++i) {
            ants_.push_back( std::unique_ptr<Ant>{ create_ant(problem_) } );
        }
        // Place each ant on a random start node
        for (auto &ant : ants_) {
            ant->move_to(rng_.rand_uint(0u, dimension-1));
        }

        for (auto i = 0u; i < dimension-1; ++i) {
            for (auto &ant : ants_) {
                perform_ant_move(ant.get());
                local_update(ant.get(), i, param_.phi_);
            }
        }
        // Local evaporation for last (closing) edge of each tour
        for (auto &ant : ants_) {
            assert(ant->has_complete_solution());
            local_update(ant.get(), dimension-1, param_.phi_);
        }
        apply_local_search();
        eval_ants_solutions();
    }


    bool local_update(const Ant *ant, uint32_t node_index, double evap_ratio) {
        return pheromone_memory_->local_update(ant->get_visited(),
                                               node_index, evap_ratio);
    }


    void eval_ants_solutions() {
        for (auto &ant : ants_) {
            ant->set_value( problem_.eval_solution( ant->get_visited() ));
        }
    }


    void select_global_best_ant() {
        // Now select best ant
        using namespace std;
        auto best_ant = get_best_ant(ants_);
        if (global_best_ant_.get_value() > best_ant->get_value()) {
            if (!problem_.is_solution_valid(best_ant->get_visited())) {
                std::cerr << "Invalid solution!" << std::endl;
                exit(EXIT_FAILURE);
            }
            global_best_ant_ = *best_ant;
            global_best_found_iter_ = stop_condition_->get_iteration();
            if (param_.verbose_) {
                std::cout << "Global best ant: " << global_best_ant_.get_value() << std::endl;
            }
        }
    }


    /*
     *  Move ant to a next node selected using pseudorandom proportional rule.
     *  Current implementation uses also candidates lists as proposed by M. Dorigo
     */
    void perform_ant_move(Ant *ant) {
        auto next_node = node_sel_impl_->select_next_node(*pheromone_memory_,
                                                            ant, param_.q0_, rng_);
        assert(next_node < problem_.get_dimension());
        ant->move_to(next_node);
    }


    void apply_local_search() {
        if (local_search_) {
            local_search_->apply(ants_, ants_.size(), &global_best_ant_);
        }
    }


    Parameters param_;
    std::shared_ptr<NodeSelectionImpl> node_sel_impl_;
    std::shared_ptr<PheromoneMemory> pheromone_memory_;
    Problem_type &problem_;
    Ant global_best_ant_;
    uint64_t global_best_found_iter_;
    std::vector<std::unique_ptr<Ant>> ants_; // Pointers to Ants of the current iteration
    std::shared_ptr<ProgressVisitor<ACS<NodeSelectionImpl, PheromoneMemory>>> progress_visitor_ = nullptr;
    std::shared_ptr<StopCondition> stop_condition_ = nullptr;
};

#endif /* ifndef ACS_H */
