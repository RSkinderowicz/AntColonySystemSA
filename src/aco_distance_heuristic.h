#ifndef ACO_DISTANCE_HEURISTIC
#define ACO_DISTANCE_HEURISTIC

#include <vector>
#include <cmath>

/*
 * ACO algorithms solving TSP and related problems construct solutions to the
 * problem based on pheromone memory and knowledge about the problem instance in
 * the form of a distance matrix between pairs of cities (nodes). The distances
 * are used to calculate a heuristic value playing a role of attractiveness of
 * an edge (connection) between the pair of cities.
 *
 * ACODistanceHeuristic class simplifies calculation of the heuristic values for
 * the given TSP-like problem instance.
 * Heuristic value for the edge (u, v) is equal to [1/distance(u,v)]^beta.
 */
template<typename Problem>
class ACODistanceHeuristic {
public:

    ACODistanceHeuristic(Problem &problem,
                         double beta) :
        problem_(problem),
        beta_(beta)
    {
        if (!problem_.is_large_instance()) {
            init_heuristic_matrix();
        }
    }


    inline double get(uint32_t u, uint32_t v) const {
        return use_heuristic_matrix_ ? heuristic_matrix_[u][v]
                                     : calc_heuristic(u, v);
    }


    double calc_heuristic(uint32_t i, uint32_t j) const {
        auto dist = problem_.get_distance(i, j);
        if (i != j && dist > 0)
            return 1.0 / fast_pow(dist, beta_);
        return 1.0;
    }

private:

    void init_heuristic_matrix() {
        auto n = problem_.get_dimension();

        heuristic_matrix_.resize(n);
        for (auto i = 0u; i < n; ++i) {
            auto &row = heuristic_matrix_[i];
            for (auto j = 0u; j < n; ++j) {
                row.push_back(calc_heuristic(i, j));
            }
        }
        use_heuristic_matrix_ = true;
    }


    /**
     * fast_pow A bit faster implementation of power when the exponent is a positive integer.
     *
     * @return base to the power of exp
     */
    static double fast_pow(double base, double exp) {
        if ((exp - (int)exp) < 1.0e-10) {
            int iexp = (int)exp;
            switch(iexp) {
                case 1: return base;
                case 2: return base * base;
                case 3: return base * base * base;
                case 4:
                    {
                        const double b2 = base * base;
                        return b2 * b2;
                    }
                case 5:
                    {
                        const double b2 = base * base;
                        return b2 * b2 * base;
                    }
                default:
                    {
                        double result = 1;
                        while (iexp) {
                            if (iexp & 1) {
                                result *= base;
                            }
                            iexp >>= 1;
                            base *= base;
                        }
                        return result;
                    }
            }
        } else {
            return pow(base, exp);
        }
    }


    const Problem &problem_;
    double beta_;
    bool use_heuristic_matrix_ = false;
    std::vector<std::vector<double>> heuristic_matrix_;
};

#endif /* ifndef ACO_DISTANCE_HEURISTIC */
