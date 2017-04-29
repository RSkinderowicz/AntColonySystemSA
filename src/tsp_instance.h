#ifndef TSP_INSTANCE_H
#define TSP_INSTANCE_H

#include <vector>
#include <memory>
#include <cassert>
#include <algorithm>

/**
 * A class managing a definition of the TSP instance.
 */
class TSPInstance {
public:
    enum class EdgeWeightType { EUC_2D, GEO };
    using Point2D = std::pair<double, double>;

    static const uint32_t LARGE_INSTANCE_THRESHOLD = 6000;


    TSPInstance(std::vector<Point2D> coords, EdgeWeightType weight_type,
                std::string name);


    static std::shared_ptr<TSPInstance> load_from_tsplib_file(std::string path);


    uint32_t get_dimension() const { return dimension_; }


    double get_distance(uint32_t u, uint32_t v) const {
        if (use_dist_matrix_) {
            return dist_matrix_[u][v];
        }
        return calc_distance(coords_[u], coords_[v]);
    }


    double eval_solution(const std::vector<uint32_t> &route) const {
        auto dist = 0.0;
        auto u = route.back();
        for (auto v : route) {
            dist += get_distance(u, v);
            u = v;
        }
        return dist;
    }

    /*
    For the given city returns the list of its nn_count closest neighbours.
    */
    std::vector<uint32_t> get_nearest_neighbors(uint32_t city,
                                                uint32_t nn_count) const {
        assert(city < dimension_);
        assert(nn_count < dimension_);

        std::vector<uint32_t> nodes;
        nodes.reserve(dimension_ - 1);
        for (auto i = 0u; i < dimension_; ++i) {
            if (i != city) { nodes.push_back(i); }
        }
        auto self = this;
        std::partial_sort(std::begin(nodes), std::begin(nodes) + nn_count, std::end(nodes),
             [&] (uint32_t a, uint32_t b) { return self->get_distance(city, a) < self->get_distance(city, b); });
        nodes.resize(nn_count);
        return nodes;
    }


    /*
    Builds and returns a complete TSP solution build starting from node 0 and
    going though the closest unvisited node, and so on.
    */
    std::vector<uint32_t> build_greedy_solution() {
        std::vector<uint32_t> nodes(dimension_-1);
        std::iota(std::begin(nodes), std::end(nodes), 1);
        std::vector<uint32_t> route;
        route.push_back(0);
        auto curr = *route.begin();
        auto self = this;
        while (!nodes.empty()) {
            auto closest = std::min_element(std::begin(nodes), std::end(nodes),
                    [&](uint32_t a, uint32_t b) {
                        return self->get_distance(curr, a) < self->get_distance(curr, b);
                    });
            route.push_back(*closest);
            curr = *closest;
            nodes.erase(closest);
        }
        return route;
    }


    /*
     * Name of the problem instance as read from the problem description file.
     */
    std::string get_name() const { return name_; }


    /*
     * We assume that the TSP is symmetric, i.e. distance between any two cities
     * (a, b) is equal to distance between (b, a).
     *
     * The assymetric version, ATSP, does not have this property.
     */
    bool is_symmetric() const { return true; }


    /*
     * Returns relative error value or positive infinity if the optimum value is
     * not known.
     */
    double calc_relative_error(double sol_value) const;


    void set_optimum_value(double value) {
        optimum_value_ = value;
    }


    bool is_large_instance() const { return dimension_ >= LARGE_INSTANCE_THRESHOLD; }


    bool is_solution_valid(const std::vector<uint32_t>& ) const  { return true; }


    double get_optimum_value() const { return optimum_value_; }


    const std::vector<std::vector<double>>& get_distance_matrix() const {
        return dist_matrix_;
    }

private:

    void init_distance_matrix();

    double calc_distance(Point2D, Point2D) const;


    std::vector<Point2D> coords_;
    uint32_t dimension_ = 0;
    EdgeWeightType edge_weight_type_;
    bool use_dist_matrix_ = false;
    std::vector<std::vector<double>> dist_matrix_;
    std::string name_; // optional instance name
    double optimum_value_ = -1;
};


#endif /* ifndef TSP_INSTANCE_H */
