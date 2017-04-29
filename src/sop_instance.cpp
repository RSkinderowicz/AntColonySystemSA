#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>

#include "sop_instance.h"
#include "ant_sop.h"

using namespace std;


int get_optimum(const std::string &instance_name);


vector<string> split(const string &str, const string &sep) {
    vector<string> res;

    size_t pos;
    size_t prev_pos = 0;
    while ( (pos = str.find(sep, prev_pos)) != string::npos ) {
        res.push_back(str.substr(prev_pos, pos-prev_pos));
        prev_pos = pos + sep.size();
    }
    res.push_back(str.substr(prev_pos, str.size()-prev_pos));
    return res;
}


string strip(const string &str) {
    size_t beg = 0;
    size_t end = str.size()-1;
    while (beg < end) {

        if (isspace(str[beg])) {
            ++beg;
        } else if (isspace(str[end])) {
            --end;
        } else {
            break;
        }
    }
    return str.substr(beg, end-beg+1);
}


bool read_pair(const string &str, const string &first, string &second) {
    vector<string> tokens = split(str, ":");

    if (strip(tokens[0]) == first && tokens.size() == 2) {
        second = strip(tokens[1]);
        return true;
    }
    return false;
}


/**
 * Returns lower case version of the string.
 */
std::string to_lower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}


SOPInstance::SOPInstance(const std::vector<std::vector<int>> &dist_matrix) :
    dimension_(dist_matrix.size()),
    dist_matrix_(dist_matrix) {

    incoming_count_.resize(dimension_, 0);
    outgoing_.resize(dimension_);
    incoming_.resize(dimension_);

    for (size_t u = 0, n = dimension_; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            if (get_distance(u, v) < 0) { // v should be before u
                incoming_count_[u] += 1;

                outgoing_[v].push_back(u);
                incoming_[u].push_back(v);
            }
        }
    }
}


std::shared_ptr<SOPInstance> SOPInstance::load_from_file(std::string path) {
    std::shared_ptr<SOPInstance> result { nullptr };
    bool is_ok = true;

    cout << "Loading test data from file: " << path << endl;

    ifstream in(path.c_str(), ios::in);

    if (in.is_open()) {
        /*
        Header has the following format:

        NAME: ESC78.sop
        TYPE: SOP
        COMMENT: Received by Norbert Ascheuer / Laureano Escudero
        DIMENSION: 80
        EDGE_WEIGHT_TYPE: EXPLICIT
        EDGE_WEIGHT_FORMAT: FULL_MATRIX 
        EDGE_WEIGHT_SECTION
        80
        ...
        */
        string line;
        getline(in, line);
        string name;
        is_ok = read_pair(line, "NAME", name);
        if ( is_ok ) {
            in.seekg( 0 );
            result = read_using_tsplib_format( in );
        } else {
            in.seekg( 0 );
            result = read_using_soplib_format( in );
            const auto name_with_ext = path.substr(path.find_last_of("/")+1);
            const auto bare_name = name_with_ext.substr(0, name_with_ext.find_last_of('.'));
            result->name_ = bare_name;
        }
        if (result) {
            auto optimum = get_optimum(result->get_name());
            if (optimum != -1) {
                result->set_optimum_value(optimum);
            }
        }

    } else {
        cout << "File does not exist!" << endl;
        is_ok = false;
    }
    in.close();

    return result;
}


std::shared_ptr<SOPInstance> SOPInstance::read_using_tsplib_format( istream &in ) {
    std::shared_ptr<SOPInstance> instance;
    bool is_ok = true;
    /*
    Header has the following format:

    NAME: ESC78.sop
    TYPE: SOP
    COMMENT: Received by Norbert Ascheuer / Laureano Escudero
    DIMENSION: 80
    EDGE_WEIGHT_TYPE: EXPLICIT
    EDGE_WEIGHT_FORMAT: FULL_MATRIX 
    EDGE_WEIGHT_SECTION
    80
    ...
    */
    string line;
    getline(in, line);
    string name;
    is_ok = read_pair(line, "NAME", name);
    cout << "Name: " << name << endl;

    getline(in, line); // skip type
    getline(in, line); // skip comment

    getline(in, line); // dimension
    string dim_str;
    is_ok &= read_pair(line, "DIMENSION", dim_str);

    uint32_t dimension = (uint32_t)stoi(dim_str);
    cout << "Dimension: " << dimension << endl;

    getline(in, line); // now edge weight type
    string ewt;
    is_ok &= read_pair(line, "EDGE_WEIGHT_TYPE", ewt);

    if (to_lower(ewt) != "explicit") {
        is_ok = false;
        cout << "Unknown edge weight type" << endl;
    }

    getline(in, line); // now edge weight type
    string ewf;
    is_ok &= read_pair(line, "EDGE_WEIGHT_FORMAT", ewf);

    if (to_lower(ewf) != "full_matrix") {
        throw runtime_error("SOPInstance::read_using_tsplib_format: Unknown edge weight format");
    }
    getline(in, line); // skip "EDGE_WEIGHT_SECTION" line
    getline(in, line); // skip dimension repeat

    vector<vector<int>> dist_matrix(dimension);

    // now read distance matrix
    for (auto &dist_row : dist_matrix) {
        getline(in, line);
        istringstream ins(line);
        int v;
        for (auto j = 0u; j < dimension; ++j) {
            ins >> v;
            dist_row.push_back(v);
            if (ins.fail()) {
                throw runtime_error("SOPInstance::read_using_tsplib_format: error while reading dist. matrix");
            }
        }
    }
    cout << "Distance matrix read" << is_ok << endl;
    instance = make_shared<SOPInstance>(dist_matrix);
    const auto bare_name = name.substr(0, name.find_last_of('.'));
    instance->name_ = bare_name;
    return instance;
}


std::shared_ptr<SOPInstance> SOPInstance::read_using_soplib_format( std::istream &in ) {
    std::shared_ptr<SOPInstance> instance;
    bool is_ok = true;
    auto non_empty_lines = 0u;
    in.seekg(0);
    string line;

    while (in.good()) {
        getline(in, line);
        if (strip(line).size() > 0) {
            ++non_empty_lines;
        }
    }
    in.clear();
    in.seekg(0);

    auto dimension = non_empty_lines;
    cout << "Dimension: " << dimension << endl;

    vector<vector<int>> dist_matrix(dimension);

    // now read distance matrix
    for (auto &dist_row : dist_matrix) {
        getline(in, line);
        istringstream ins(line);
        int v;
        for (auto j = 0u; j < dimension; ++j) {
            ins >> v;
            dist_row.push_back(v);
            if (ins.fail()) {
                throw runtime_error("SOPInstance::read_using_soplib_format: error while reading dist. matrix");
            }
        }
    }
    cout << "Distance matrix read" << is_ok << endl;
    return make_shared<SOPInstance>(dist_matrix);
}


double SOPInstance::calc_relative_error(double sol_value) {
    if (optimum_value_ > 0.0) {
        return (sol_value - optimum_value_) / optimum_value_;
    } else {
        return std::numeric_limits<double>::infinity();
    }
}


double SOPInstance::eval_solution(const std::vector<uint32_t> &route) const {
    auto dist = 0.0;
    auto it = route.begin();
    auto u = *it;
    // Route is not closed in the case of SOP
    for (++it; it != route.end(); ++it) {
        auto v = *it;
        dist += get_distance(u, v);
        u = v;
    }
    return dist;
}


/**
 * A method for building an arbitrary initial solution. The only requirement
 * is that the solution has to be valid.
 */
std::vector<uint32_t> SOPInstance::build_greedy_solution() const {
    AntSOP ant(*this);

    ant.move_to(0); // source

    std::vector<uint32_t> nodes(dimension_-1);
    std::iota(std::begin(nodes), std::end(nodes), 1);

    while (!nodes.empty()) {
        auto closest = end(nodes);
        auto pos = ant.get_position();
        for (auto it = begin(nodes); it != end(nodes); ++it) {
            if (ant.is_available(*it) &&
                    (closest == end(nodes) ||
                     get_distance(pos, *it) < get_distance(pos, *closest))) {
                closest = it;
            }
        }
        ant.move_to(*closest);
        nodes.erase(closest);
    }
    auto sol = ant.get_visited();
    assert(is_solution_valid(sol));
    return sol;
}


/**
 * @return true if the given solution is valid, false otherwise.
 */
bool SOPInstance::is_solution_valid(const std::vector<uint32_t> &sol) const {
    if (sol.size() != (size_t)dimension_) {
        return false;
    }
    vector<bool> seen(dimension_, 0);
    const auto n = sol.size();
    for (auto i = 0u; i < n; ++i) {
        auto node = sol[i];

        if (node >= dimension_) {
            return false;
        }
        if (seen[node]) {
            return false;
        }
        seen[node] = 1;
    }
    // Now check if all precedence constraints are satisfied
    for (auto i = 0u; i < n; ++i) {
        auto node_a = sol[i];

        for (auto j = i + 1; j < n; ++j) {
            auto node_b = sol[j];

            if (get_distance(node_a, node_b) < 0) {
                // node_b should be before node_a in the solution
                return false;
            }
        }
    }
    return true;
}


/*
 * Best know solutions as of 13.01.2016
 * Based on:
 * - Coupling ant colony systems with strong local searches LM Gambardella, R
 *   Montemanni, D Weyland - European Journal of Operational Research, 2012
 * - Load-dependent and precedence-based models for pickup and delivery problems
 *   L Gouveia, M Ruthmair - Computers & Operations Research, 2015
 */
std::map<std::string, int> SOP_optima = {
    {"esc07", 2125},
    {"esc12", 1675},
    {"esc25", 1681},
    {"esc47", 1288},
    {"esc63", 62},
    {"esc78", 18320},
    {"br17.10", 55},
    {"br17.12", 55},
    {"ft53.1", 7531},
    {"ft53.2", 8026},
    {"ft53.3", 10262},
    {"ft53.4", 14425},
    {"ft70.1", 39313},
    {"ft70.2", 40419},
    {"ft70.3", 42535},
    {"ft70.4", 53530},
    {"kro124p.1", 39420},
    {"kro124p.2", 41336},
    {"kro124p.3", 49499},
    {"kro124p.4", 76103},
    {"p43.1", 28140},
    {"p43.2", 28480},
    {"p43.3", 28835},
    {"p43.4", 83005},
    {"prob.42", 243},
    {"prob42", 243},
    {"prob.100", 1163},
    {"rbg048a", 351},
    {"rbg050c", 467},
    {"rbg109a", 1038},
    {"rbg150a", 1750},
    {"rbg174a", 2033},
    {"rbg253a", 2950},
    {"rbg323a", 3140},
    {"rbg341a", 2568},
    {"rbg358a", 2545},
    {"rbg378a", 2816},
    {"ry48p.1", 15805},
    {"ry48p.2", 16666},
    {"ry48p.3", 19894},
    {"ry48p.4", 31446},
    {"r.200.100.1", 61},
    {"r.200.100.15", 1792},
    {"r.200.100.30", 4216},
    {"r.200.100.60", 71749},
    {"r.200.1000.1", 1404},
    {"r.200.1000.15", 20481},
    {"r.200.1000.30", 41196},
    {"r.200.1000.60", 71556},
    {"r.300.100.1", 26},
    {"r.300.100.15", 3161},
    {"r.300.100.30", 6120},
    {"r.300.100.60", 9726},
    {"r.300.1000.1", 1294},
    {"r.300.1000.15", 29183},
    {"r.300.1000.30", 54147},
    {"r.300.1000.60", 109471},
    {"r.400.100.1", 13},
    {"r.400.100.15", 3906},
    {"r.400.100.30", 8165},
    {"r.400.100.60", 15228},
    {"r.400.1000.1", 1343},
    {"r.400.1000.15", 43268},
    {"r.400.1000.30", 85128},
    {"r.400.1000.60", 140816},
    {"r.500.100.1", 4},
    {"r.500.100.15", 5361},
    {"r.500.100.30", 9665},
    {"r.500.100.60", 18240},
    {"r.500.1000.1", 1316},
    {"r.500.1000.15", 50725},
    {"r.500.1000.30", 98987},
    {"r.500.1000.60", 178212},
    {"r.600.100.1", 1},
    {"r.600.100.15", 5684},
    {"r.600.100.30", 12465},
    {"r.600.100.60", 23293},
    {"r.600.1000.1", 1337},
    {"r.600.1000.15", 57237},
    {"r.600.1000.30", 126798},
    {"r.600.1000.60", 214608},
    {"r.700.100.1", 1},
    {"r.700.100.15", 7311},
    {"r.700.100.30", 14510},
    {"r.700.100.60", 24102},
    {"r.700.1000.1", 1231},
    {"r.700.1000.15", 66837},
    {"r.700.1000.30", 134474},
    {"r.700.1000.60", 245589}
};


int get_optimum(const std::string &instance_name) {
    const auto key = to_lower(instance_name);
    if (SOP_optima.find(key) != SOP_optima.end()) {
        return SOP_optima[key];
    }
    return -1;
}
