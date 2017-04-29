#ifndef ACS_PARAMS
#define ACS_PARAMS

#include <cstdint>

/*
 * A collection of the standard ACS parameters to simplify passing to functions
 * & methods.
*/
struct ACSParameters {
    uint32_t ants_count_ = 10;
    double beta_ = 3.0;
    uint32_t cand_list_size_ = 25;
    double rho_ = 0.1; // global evaporation rate
    double phi_ = 0.01; // local evaporation rate
    double q0_ = 0.9;

    // Non-optimization related parameters
    bool verbose_ = true;
};


/*
 * A collection of the standard ACS parameters to simplify passing to functions
 * & methods.
*/
struct ACSSAParameters {
    uint32_t ants_count_;
    double beta_ = 3.0;
    uint32_t cand_list_size_ = 25;
    double rho_ = 0.1; // global evaporation rate
    double phi_ = 0.01; // local evaporation rate
    double q0_ = 0.9;
    // Non-optimization related parameters
    bool verbose_ = true;

    // SA-related parameters
    double sa_cooling_ratio_ = 0.999;
    double sa_initial_accept_prob_ = 0.3;
    double sa_min_temperature_multiplier_ = 1e-4;
    uint32_t sa_temp_init_moves_ = 1000u;
};


#endif /* ifndef ACS_PARAMS */
