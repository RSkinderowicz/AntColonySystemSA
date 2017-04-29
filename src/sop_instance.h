#ifndef SOP_INSTANCE_H
#define SOP_INSTANCE_H

#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>
#include <iostream>
#include <map>

/**
 * A class for the Sequential Ordering Problem definition.
 */
class SOPInstance {
public:
    static const uint32_t LARGE_INSTANCE_THRESHOLD = 6000;


    SOPInstance(const std::vector<std::vector<int>> &dist_matrix);


    static std::shared_ptr<SOPInstance> load_from_file(std::string path);


    uint32_t get_dimension() const { return dimension_; }

    /*
     * SOP instances are by definition assymetric
     */
    bool is_symmetric() const { return false; }

    /*
     * Returns relative error value or positive infinity if the optimum value is
     * not known.
     */
    double calc_relative_error(double sol_value);


    int get_distance(uint32_t u, uint32_t v) const {
        assert(u < dimension_);
        assert(v < dimension_);
        return dist_matrix_[u][v];
    }

    /*
    For the given city returns the list of its nn_count closest neighbours.
    */
    std::vector<uint32_t> get_nearest_neighbors(uint32_t city,
                                                uint32_t nn_count) const {
        assert(city < dimension_);
        assert(nn_count < dimension_);

        std::vector<uint32_t> nodes;
        for (auto i = 0u; i < dimension_; ++i) {
            if (i != city && get_distance(city, i) >= 0) {
                nodes.push_back(i);
            }
        }
        if (nn_count > nodes.size()) {
            nn_count = nodes.size();
        }
        auto self = this;
        std::partial_sort(std::begin(nodes), std::begin(nodes) + nn_count, std::end(nodes),
             [&self, &city] (uint32_t a, uint32_t b) {
                 return self->get_distance(city, a) < self->get_distance(city, b);
             });

        nodes.resize(nn_count);
        return nodes;
    }


    /*
     * Name of the problem instance as read from the problem description file.
     */
    std::string get_name() const { return name_; }


    void set_optimum_value(double value) {
        optimum_value_ = value;
    }


    double eval_solution(const std::vector<uint32_t> &route) const;


    const std::vector<uint32_t>& get_incoming_edges(uint32_t node) const { return incoming_.at(node); }


    const std::vector<uint32_t>& get_incoming_edges_count() const { return incoming_count_; }


    const std::vector<uint32_t>& get_outgoing_edges(uint32_t node) const { return outgoing_.at(node); }


    std::vector<uint32_t> build_greedy_solution() const;


    bool is_solution_valid(const std::vector<uint32_t> &sol) const;


    bool is_large_instance() const { return dimension_ > LARGE_INSTANCE_THRESHOLD; }


    double get_optimum_value() const { return optimum_value_; }

private:

    static std::shared_ptr<SOPInstance> read_using_tsplib_format(std::istream &in);

    static std::shared_ptr<SOPInstance> read_using_soplib_format( std::istream &in );


    uint32_t dimension_;
    std::vector<std::vector<int>> dist_matrix_;
    std::string name_; // optional instance name
    double optimum_value_ = -1;

    // outgoing_[v] = a list of nodes u such that edge (v, u) means that
    // node v has to be before node u in the solution
    std::vector<std::vector<uint32_t>> outgoing_;
    // incoming_[v] = a list of nodes u such that edge (u, v) means that
    // node u has to be before node v in the solution
    std::vector<std::vector<uint32_t>> incoming_;
    // incoming_count_[v] = a number of incoming edges for the node v
    std::vector<uint32_t> incoming_count_;
};



#endif /* ifndef SOP_INSTANCE_H */

