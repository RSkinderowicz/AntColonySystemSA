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

#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <string>
#include <fstream>


#include "cmd_line_args.h"
#include "tsp_instance.h"
#include "standard_pheromone_memory.h"
#include "lru_pheromone_memory.h"
#include "acs_node_selection_impl.h"
#include "acs.h"
#include "tsplib.h"
#include "k_opt.h"
#include "sop_instance.h"
#include "sop_local_search.h"
#include "log.h"
#include "acs_sa.h"
#include "fs_utils.h"
#include "problem_model.h"
#include "eacs.h"
#include "eacs_sa.h"

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>


void handle_SIGSEV(int sig) {
    void *array[10];

    // get void*'s for all entries on the stack
    auto size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

using namespace std;


uint32_t CCounter = 0;

/** A helper class for tracking the progress of the ACS algorithm, and gathering
 *  some basic statistics like runtime and best solution found.
 */
template<typename ACS_t>
class RecordBasicRunProgressACS : public ProgressVisitor<ACS_t> {
public:

    RecordBasicRunProgressACS() :
        best_solution_value_(std::numeric_limits<double>::max()),
        iter_when_best_found_(0) {
        log_open_dict();
    }


    virtual void started(ACS_t &/*alg*/) override {
        start_ = std::chrono::steady_clock::now();
    }


    virtual void iteration_done(ACS_t &alg) override {
        if (alg.get_current_iter() % best_value_record_period_ == 0) {
            best_values_track_.push_back(alg.get_global_best_ant().get_value());
        }
    }


    virtual void finished(ACS_t &alg) override {
        end_ = std::chrono::steady_clock::now();
        best_solution_value_ = alg.get_global_best_ant().get_value();
        best_solution_ = alg.get_global_best_ant().get_visited();
        iter_when_best_found_ = alg.get_global_best_found_iter();
        auto elapsed = chrono::duration_cast<chrono::microseconds>(end_ - start_);
        auto calc_time = elapsed.count() / 1.0e6;
        cout << "Calc. time: " << calc_time << " sec" << endl;
        log_add("calc_time", calc_time);
        log_add("best_solution_value", best_solution_value_);
        log_add("iter_when_best_found", iter_when_best_found_);
        log_add("best_solution", best_solution_);
        log_add("iter_made", alg.get_current_iter());
        log_add("best_solutions_track", best_values_track_);

        if (alg.get_local_search() != nullptr) {
            alg.get_local_search()->log_statistics();
        }
        log_close_dict();
    }

protected:

    double best_solution_value_;
    std::vector<uint32_t> best_solution_;
    uint32_t iter_when_best_found_;
    std::chrono::time_point<std::chrono::steady_clock> start_;
    std::chrono::time_point<std::chrono::steady_clock> end_;
    uint64_t best_value_record_period_ = 10;
    std::vector<double> best_values_track_;
};


/** This is for reporting the progress & stats of ACS combined with the
 *  Simulated Annealing */
template<typename ACS_t>
class RecordBasicRunProgressACS_SA : public RecordBasicRunProgressACS<ACS_t> {
    using ParentClass = RecordBasicRunProgressACS<ACS_t>;
public:
    RecordBasicRunProgressACS_SA() {
    }

    virtual void started(ACS_t &alg) override {
        ParentClass::started(alg);
    }


    virtual void iteration_done(ACS_t &alg) override {
        ParentClass::iteration_done(alg);
    }


    virtual void finished(ACS_t &alg) override {
        const auto sa = alg.get_simulated_annealing_control();
        const auto param = alg.get_parameters();
        // Simulated annealing related parameters & stats.
        log_open_dict("simulated_annealing");
        log_add("cooling_ratio", param.sa_cooling_ratio_);
        log_add("initial_accept_prob", param.sa_initial_accept_prob_);
        log_add("initial_temperature", sa.initial_temperature_);
        log_add("min_temperature", sa.min_temperature_);
        log_add("min_temperature_multiplier_", param.sa_min_temperature_multiplier_);
        log_add("temp_init_moves",
                param.sa_temp_init_moves_);
        log_add("worse_accepted", sa.worse_accepted_);
        log_add("worse_rejected", sa.worse_rejected_);
        log_close_dict();

        ParentClass::finished(alg);
    }
};


std::shared_ptr<TSPInstance> read_tsp_instance(std::string path) {
    auto tsp_instance = TSPInstance::load_from_tsplib_file(path);
    auto it = TSPLIB_optima.find(tsp_instance->get_name());
    if (it != TSPLIB_optima.end()) {
        tsp_instance->set_optimum_value(it->second);
    }
    return tsp_instance;
}


std::shared_ptr<SOPInstance> read_sop_instance(std::string path) {
    auto instance = SOPInstance::load_from_file(path);
    return instance;
}


/** Default progress visitor */
template<typename Algorithm>
void create_progress_visitor(Algorithm &alg) {
    alg.set_progress_visitor( make_shared<RecordBasicRunProgressACS<Algorithm>>() );
}


/** Progress visitior for the ACS with the Simulated Annealing */
template<typename A, typename B>
void create_progress_visitor(ACS_SA<A,B> &alg) {
    alg.set_progress_visitor( make_shared<RecordBasicRunProgressACS_SA<ACS_SA<A,B>>>() );
}


/*
 * A structure used to carry necessary info about which version of the algorithm
 * to run etc.
 */
template<typename Problem,
         typename PheromoneMemory,
         template<typename, typename> class ACOAlgorithm>
class Experiment {
public:

    using Problem_type = Problem;
    using PheromoneMemory_type = PheromoneMemory;
    using NodeSelectionImpl_type = ACSNodeSelectionImpl<ProblemModel<Problem_type>>;
    using Algorithm_type = ACOAlgorithm<NodeSelectionImpl_type, PheromoneMemory_type>;

    using Parameters_t = typename Algorithm_type::Parameters;


    Experiment(std::shared_ptr<Problem> problem,
               RNG &rng) :
        problem_(problem),
        rng_(rng) {
        initial_pheromone_ = calc_initial_pheromone(*problem);
    }


    template<typename Alg>
    void execute(std::shared_ptr<Alg> alg,
                        std::shared_ptr<StopCondition> stop_condition,
                        uint32_t trials,
                        std::shared_ptr<LocalSearch> local_search = nullptr) {
        vector<double> results;

        auto start = chrono::steady_clock::now();

        alg->set_local_search(local_search);
        if (local_search != nullptr) {
            local_search->set_stop_condition(stop_condition);
        }

        log_open_list("runs");
        for (auto i = 0u; i < trials; ++i) {
            cout << "\n\nTrial: " << (i+1) << endl;
            create_progress_visitor(*alg);
            alg->run_begin(stop_condition);
            alg->run();
            results.push_back( alg->get_global_best_ant().get_value() );

            if (local_search != nullptr) {
                local_search->reset();
            }
        }
        log_close_list();

        auto end = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::microseconds>(end - start);
        auto elapsed_sec = elapsed.count() / 1.0e6;
        log_add("total_time_sec", elapsed_sec);

        if (!results.empty()) {
            auto mean = accumulate(results.cbegin(), results.cend(), 0.0) / results.size();
            auto min = *min_element(results.cbegin(), results.cend());
            auto max = *max_element(results.cbegin(), results.cend());
            cout << "Mean.: " << mean << "\tmin: " << min << "\tmax: " << max << endl;

            auto rel_err = problem_->calc_relative_error(min);
            cout << "Min. relative error [%]: " << rel_err * 100 << endl;
            cout << "Mean relative error [%]: " << problem_->calc_relative_error(mean) * 100 << endl;

            log_add("solution_value_min", min);
            log_add("solution_value_max", max);
            log_add("solution_value_mean", mean);
            log_add("solution_value_min_rel_err", rel_err);
        }
    }


    void run(std::shared_ptr<StopCondition> stop_condition,
             uint32_t trials, std::shared_ptr<LocalSearch> ls = nullptr) {

        Parameters_t params;
        set_parameters(params);
        execute(create_acs(params), stop_condition, trials, ls);
    }

    void run_sa_version(std::shared_ptr<StopCondition> stop_condition,
                              uint32_t trials,
                              std::shared_ptr<LocalSearch> ls = nullptr) {
        Parameters_t params;
        set_parameters(params);
        //set_sa_parameters(params);

        execute(create_acs(params), stop_condition, trials, ls);
    }

    std::shared_ptr<PheromoneMemory_type> create_ph_mem() const {
        std::shared_ptr<PheromoneMemory_type> ptr;
        create_ph_mem(ptr);
        return ptr;
    }


    template<typename PheromoneMemory_type>
    void create_ph_mem(std::shared_ptr<PheromoneMemory_type> &ptr) const {
        ptr.reset( new PheromoneMemory_type(problem_->get_dimension(),
                                            initial_pheromone_,
                                            problem_->is_symmetric()) );
    }


    void create_ph_mem(std::shared_ptr<LRUPheromoneMemory> &ptr) const {
        uint32_t bucket_size = (uint32_t)get_cmd_line_arg("--lru_bucket");
        ptr.reset( new LRUPheromoneMemory(problem_->get_dimension(), initial_pheromone_,
                                          problem_->is_symmetric(), bucket_size) );
    }


    template<typename T>
    void set_acs_parameters(T &params) {
        params.ants_count_ = (uint32_t)get_cmd_line_arg("--ants");
        auto q0 = (double)get_cmd_line_arg("--q0");
        if (q0 > 1) {
            const auto size = problem_->get_dimension();
            q0 = (size - q0) / size;
        }
        params.q0_ = q0;

        auto cand_list_size = (int)get_cmd_line_arg("--cand_list");
        if (cand_list_size == -1) {
            cand_list_size = (int)(problem_->get_dimension() - 1);
        }
        if (cand_list_size <= 0
                || cand_list_size > (int)problem_->get_dimension()-1) {
            throw runtime_error("Invalid candidate list size: " + to_string(cand_list_size));
        }
        params.cand_list_size_ = (uint32_t)cand_list_size;

        auto beta = (double)get_cmd_line_arg("--beta");
        if (beta < 0.0 || beta > 10.0) {
            throw runtime_error("Invalid beta parameter value (valid range [0, 10.0]): " + to_string(beta));
        }
        params.beta_ = beta;

        // Global evaporation rate
        auto rho = (double)get_cmd_line_arg("--rho");
        if (!(rho > 0.0 && rho < 1.0)) {
            throw runtime_error("Invalid rho (global evaporation ratio) parameter value (valid range (0, 1)): " + to_string(rho));
        }
        params.rho_ = rho;
        // Local evaporation rate
        auto phi = (double)get_cmd_line_arg("--phi");
        if (!(phi > 0.0 && phi < 1.0)) {
            throw runtime_error("Invalid phi (local evaporation ratio) parameter value (valid range (0, 1)): " + to_string(phi));
        }
        params.phi_ = phi;

        params.verbose_ = ((int)get_cmd_line_arg("--verbosity") >= 2);

        log_add("ants", params.ants_count_);
        log_add("q0", params.q0_);
        log_add("cand_list_size", params.cand_list_size_);
        log_add("beta", params.beta_);
        log_add("rho", params.rho_);
        log_add("phi", params.phi_);

        if (model_ == nullptr) {
            model_.reset( new ProblemModel<Problem>(*problem_, params) );
        }
    }


    void set_parameters(ACSParameters &params) {
        set_acs_parameters(params);
    }


    /** Sets SA related params */
    void set_parameters(ACSSAParameters &params) {
        set_acs_parameters(params); // Set common parameters
        // Now the SA related
        params.sa_cooling_ratio_ = (double)get_cmd_line_arg("--sa-cooling-ratio");
        params.sa_initial_accept_prob_ = (double)get_cmd_line_arg("--sa-init-accept-prob");
        params.sa_temp_init_moves_ = (uint32_t)get_cmd_line_arg("--sa-init-moves");
        params.sa_min_temperature_multiplier_ =
            (double)get_cmd_line_arg("--sa-min-temp-multiplier");
    }


    std::shared_ptr<Algorithm_type> create_acs(Parameters_t params) {
        auto ph_mem = create_ph_mem();
        auto km = make_shared<NodeSelectionImpl_type>(*model_, params.cand_list_size_);
        return make_shared<Algorithm_type>(params, km, ph_mem,
                                           *model_,  rng_);
    }


    static double calc_initial_pheromone(Problem &problem) {
        cout << "calc_initial_pheromone: " << endl;
        auto start = std::chrono::steady_clock::now();
        auto greedy_val = problem.eval_solution( problem.build_greedy_solution() );
        auto end = std::chrono::steady_clock::now();
        cout << "end  "<< endl;
        auto elapsed = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "Greedy sol. calc. time: " << elapsed.count() / 1.0e6 << " sec" << endl;
        return 1.0 / (problem.get_dimension() * greedy_val);
    }


    std::shared_ptr<Problem> problem_;
    std::unique_ptr<ProblemModel<Problem>> model_;
    RNG &rng_;
    double initial_pheromone_;
};



std::shared_ptr<LocalSearch> create_ls(TSPInstance &problem, string ls_name, RNG &rng) {
    std::shared_ptr<LocalSearch> ls = nullptr;
    if (ls_name == "2opt") {
        ls.reset( new Opt2(problem, rng) );
    } else if (ls_name == "3opt") {
        ls.reset( new Opt3(problem, rng) );
    } else if (ls_name != "none") {
        throw runtime_error("Unknown local search algorithm.");
    }
    return ls;
}


std::shared_ptr<LocalSearch> create_ls(SOPInstance &problem, string ls_name, RNG &rng) {
    std::shared_ptr<LocalSearch> ls = nullptr;
    if (ls_name == "sop3exchange") {
        ls.reset( new SOPLocalSearch(problem, rng) );
    } else if (ls_name == "sop3exchange-sa") {

        auto cooling_ratio = (double)get_cmd_line_arg("--sop-ls-sa-cooling-ratio");
        auto initial_accept_prob  = (double)get_cmd_line_arg("--sop-ls-sa-init-accept-prob");
        auto min_temperature_multiplier  = (double)get_cmd_line_arg("--sop-ls-sa-min-temp-multiplier");

        ls.reset( new SOPLocalSearchSA(problem, rng,
                                       cooling_ratio,
                                       initial_accept_prob,
                                       min_temperature_multiplier) );
    } else if (ls_name != "none") {
        throw runtime_error("Unknown local search algorithm. Available [sop3exchange, sop3exchange-sa]");
    }
    return ls;
}


template<typename Experiment>
void start_experiment(Experiment &exp, std::string algorithm_name) {
    const auto trials = (uint32_t)get_cmd_line_arg("--trials");
    const auto ls_name = (string)get_cmd_line_arg("--ls");
    auto ls = create_ls(*exp.problem_, ls_name, exp.rng_);

    const auto ants_count = (uint32_t)get_cmd_line_arg("--ants");
    auto iterations = (uint32_t)get_cmd_line_arg("--iter");
    const auto iter_mult = (uint32_t)get_cmd_line_arg("--iter_mult");
    const auto timeout = (double)get_cmd_line_arg("--timeout");
    shared_ptr<StopCondition> stop_condition = nullptr;

    if (timeout > 0) {
        stop_condition.reset(new TimeoutStopCondition(timeout, ants_count));
    } else if (iter_mult > 0) {
        iterations = exp.problem_->get_dimension() * iter_mult / ants_count;
        stop_condition.reset(new FixedIterationsStopCondition(iterations, ants_count));
    } else if (iterations > 0) {
        stop_condition.reset(new FixedIterationsStopCondition(iterations, ants_count));
    } else {
        throw runtime_error("No stop condition given");
    }


    if (algorithm_name == "acs"
            || algorithm_name == "eacs") {
        exp.run(stop_condition, trials, ls);
    } else if (algorithm_name == "acs_sa"
            || algorithm_name == "eacs_sa") {
        exp.run_sa_version(stop_condition, trials, ls);
    } else {
        throw runtime_error("start_experiment(): Unknown algorithm");
    }
    log_add("trials", trials);
    log_add("iter_per_trial", iterations);
    log_add("local_search", ls_name);
    log_add("timeout", timeout);
}


template<typename Problem, typename PhMem>
void run_experiment_helper(std::shared_ptr<Problem> instance,
                           RNG &rng,
                           std::string alg_to_run) {
    if (alg_to_run == "acs") {
        Experiment<Problem, PhMem, ACS> exp { instance, rng };
        start_experiment(exp, alg_to_run);
    } else if (alg_to_run == "acs_sa") {
        Experiment<Problem, PhMem, ACS_SA> exp { instance, rng };
        start_experiment(exp, alg_to_run);
    } else if (alg_to_run == "eacs") {
        Experiment<Problem, PhMem, EACS> exp { instance, rng };
        start_experiment(exp, alg_to_run);
    } else if (alg_to_run == "eacs_sa") {
        Experiment<Problem, PhMem, EACS_SA> exp { instance, rng };
        start_experiment(exp, alg_to_run);
    } else {
        throw runtime_error("Unknown algorithm");
    }
    log_problem_info(instance);
}


template<typename Problem>
void run_experiment(std::shared_ptr<Problem> instance,
                    RNG &rng,
                    std::string alg_to_run,
                    std::string ph_mem) {
    if (ph_mem == "std") {
        run_experiment_helper<Problem, StandardPheromoneMemory>(instance, rng, alg_to_run);
    } else if (ph_mem == "lru") {
        run_experiment_helper<Problem, LRUPheromoneMemory>(instance, rng, alg_to_run);
    } else {
        throw runtime_error("Unknown pheromone memory");
    }
}



unique_ptr<RNG> init_random_number_generator() {
    unique_ptr<RNG> prng;
    const auto seed = (uint32_t)get_cmd_line_arg("--seed");
    if (seed) {
        prng.reset(new RNG(seed));
    } else {
        prng.reset(new RNG());
    }
    log_add("rng_seed", prng->get_seed());
    return prng;
}


std::string get_default_results_dir() {
    std::time_t t = std::time(NULL);
    auto now_tm = std::gmtime(&t);

    string path;
    path = "results/" + to_string(now_tm->tm_year + 1900)
         + "/" + to_string(now_tm->tm_mon+1)
         + "/" + to_string(now_tm->tm_mday);
    return path;
}



std::string create_results_dir() {
    auto dir = (string)get_cmd_line_arg("--res_dir");

    if (dir == "auto") { // Create the directory
        dir = get_default_results_dir();
    }
    auto res_dir_path = dir;

    if (is_path_relative(res_dir_path)) {
        res_dir_path = get_cwd() + '/' + dir;
    }
    if (!dir_exists(res_dir_path)) {
        cout << "Creating results directory: " << res_dir_path << endl;
        create_dir_path(res_dir_path);
    }
    return res_dir_path;
}


// pla7397_2015-05-04_20_29_15-2639.js
std::string get_results_filename(const std::string &prefix, RNG &rng) {
    std::time_t now = std::time(NULL);
    auto t = std::gmtime(&now);

    string name = prefix;
    if (name.size() > 0) {
        name += "_";
    }
    name += to_string(t->tm_year + 1900);
    name += "-" + to_string(t->tm_mon+1);
    name += "-" + to_string(t->tm_mday);
    name += "_" + to_string(t->tm_hour);
    name += "_" + to_string(t->tm_min);
    name += "_" + to_string(t->tm_sec);
    name += "_" + to_string(rng.rand_uint(0, 1000));
    name += ".js"; // JSON
    return name;
}


std::string get_results_filepath(const std::string &prefix, RNG &rng) {
    std::string path;
    do {
        auto fname = get_results_filename(prefix, rng);
        path = create_results_dir();
        path += "/";
        path += fname;
    } while (file_exists(path));
    return path;
}



enum class ProblemType { TSP, SOP };


template<typename Problem>
void log_problem_info(std::shared_ptr<Problem> instance) {
    log_add("problem_name", instance->get_name());
    log_add("problem_size", instance->get_dimension());
    log_add("problem_optimum", instance->get_optimum_value());
}



std::string to_string(ProblemType t) {
    switch(t) {
    case ProblemType::TSP: return "tsp";
    case ProblemType::SOP: return "sop";
    default: return "unknown";
    }
}


ProblemType get_problem_type_from_file_path(std::string path) {
    auto ends_with = [](std::string str, std::string suffix) {
        auto pos = str.rfind(suffix);
        if (pos != std::string::npos) {
            return pos + suffix.size() == str.size();
        }
        return false;
    };

    ProblemType problem_type { ProblemType::TSP };
    if (ends_with(path, ".tsp")) {
        problem_type = ProblemType::TSP;
    } else if (ends_with(path, ".sop")) {
        problem_type = ProblemType::SOP;
    }
    return problem_type;
}


static const char USAGE[] =
R"(ACS framework for solving TSP and SOP problems.

Usage:
    acs [options]

Options:
    --alg=Name                  Name of the algorithm to run [default: acs]
    --ants=Int                  Number of ants [default: 20]
    --beta=Real                 beta ACO parameter [default: 3]
    --cand_list=Int             Size of the candidates list (if -1 then n-1 is used) [default: 32]
    -h --help                   Show this screen.
    --iter=Int                  Number of alg. iterations [default: 1000]
    --iter_mult=Int             Number of alg. iteraions = problem_size * iter_mult / ants [default: 0]
    --label=Name                Optional label for the experiment, written to the results file [default: none]
    --lru_bucket=Int            Size of the bucket in the LRU pheromone memory [default: 8]
    --ls=Name                   Name of the local search alg. to use [default: none]
    --phi=Real                  Local pheromone evaporation ratio [default: 0.01]
    --phmem=Name                Pheromone memory type to use [default: std]
    --q0=Real                   q0 parameter (if > 1 then q := (n - q0)/n) [default: 20]
    --res_dir=Path              Path where the results should be written [default: auto]
    --rho=Real                  Global pheromone evaporation ratio [default: 0.1]
    --seed=Int                  Initial seed for the PRNG (if 0 then default) [default: 0]
    --test=Path                 Path to problem instance.
    --timeout=Real              Algorithm timeout - 0 means infinite [default: 0]
    --trials=Int                Number of trials [default: 1]
    --verbosity=Int             Limits the number of messages printed to console, from 0 to 2 [default: 1]
    --sa-cooling-ratio=Real     Sim. Annealing temp. cooling ratio [default: 0.999]
    --sa-init-accept-prob=Real  Sim. Annealing initial accept prob. [default: 0.9]
    --sa-init-moves=Int         Sim. Annealing temp. calc. non-improving moves [default: 1000]
    --sa-min-temp-multiplier=Real         Sim. Annealing minimum temp. multiplier [default: 0.0]
    --sop-ls-sa-cooling-ratio=Real        Cooling ratio for the SOP-3-exchange-SA [default: 0.99]
    --sop-ls-sa-init-accept-prob=Real     Cooling init. accept prob. for the SOP-3-exchange-SA [default: 0.1]
    --sop-ls-sa-min-temp-multiplier=Real  Min. temp. multiplier for the SOP-3-exchange-SA [default: 0.0]
)";


int main(int argc, char *argv[]) {
    signal(SIGSEGV, handle_SIGSEV);   // install handler to print stack trace

    init_cmd_line_args(USAGE, argc, argv);

    log_init();
    log_open_dict("cmd_line_args");
    for (auto &arg : get_cmd_line_args()) {
        log_add(arg.first, arg.second);
    }
    log_close_dict();

    unique_ptr<RNG> prng = init_random_number_generator();
    auto &rng = *prng;

    const auto alg_to_run = (string)get_cmd_line_arg("--alg");
    const auto phmem = (string)get_cmd_line_arg("--phmem");

    const auto path = (string)get_cmd_line_arg("--test");
    ProblemType problem_type = get_problem_type_from_file_path(path);

    string instance_name;
    if (problem_type == ProblemType::TSP) {
        auto instance = read_tsp_instance(path);
        instance_name = instance->get_name();
        run_experiment(instance, rng, alg_to_run, phmem);
    } else if (problem_type == ProblemType::SOP) {
        auto instance = read_sop_instance(path);
        instance_name = instance->get_name();
        run_experiment(instance, rng, alg_to_run, phmem);
    }

    log_add("problem_type", to_string(problem_type));
    log_add("problem_path", path);
    log_add("algorithm", alg_to_run);
    log_add("ph_mem", phmem);
    log_add("label", (string)get_cmd_line_arg("--label"));

    if ((int)get_cmd_line_arg("--verbosity") >= 2) {
        cout << "Log:" << log_get_contents() << endl;
    }
    create_results_dir();
    const auto prefix = to_string(problem_type)
                + "_" + instance_name
                + "_" + alg_to_run;
    const auto res_file_path = get_results_filepath(prefix, rng);
    cout << "Saving results to a file: " << res_file_path << endl;
    ofstream out(res_file_path);
    if (out) {
        out << log_get_contents();
    }
    return 0;
}
