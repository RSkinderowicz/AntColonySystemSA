#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <map>

#include "tsp_instance.h"

using namespace std;


static inline string &ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
    return s;
}

static inline string &rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
    return s;
}

static inline string &trim(string &s) {
    return ltrim(rtrim(s));
}


static double dtrunc(double x) { return int(x); }

/*
 *Compute geometric distance between two nodes rounded to next
 *integer for TSPLIB instances. Based on the ACOTSP by Thomas Stuetzle
*/
static long int geo_distance(double x1, double y1, double x2, double y2) {
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    long int dd;

    deg = dtrunc(x1);
    min = x1 - deg;
    lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc(x2);
    min = x2 - deg;
    latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc(y1);
    min = y1 - deg;
    longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc(y2);
    min = y2 - deg;
    longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos(longi - longj);
    q2 = cos(lati - latj);
    q3 = cos(lati + latj);
    dd = (int) (6378.388 * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;
}


TSPInstance::TSPInstance(vector<Point2D> coords, EdgeWeightType weight_type, string name)
    : coords_(coords),
      dimension_(coords.size()),
      edge_weight_type_(weight_type),
      name_(name){

    if (!is_large_instance()) {
        init_distance_matrix();
    }
}


double TSPInstance::calc_distance(Point2D p1, Point2D p2) const {
    switch(edge_weight_type_) {
        case EdgeWeightType::EUC_2D: {
            const auto dx = p2.first - p1.first;
            const auto dy = p2.second - p1.second;
            return (int) ( sqrt(dx * dx + dy * dy) + 0.5 );
        }
        case EdgeWeightType::GEO:
            return geo_distance(p1.first, p1.second, p2.first, p2.second);
        default:
            throw std::runtime_error("TSPInstance::calc_distance: Unknown edge weight type");
    }
    return std::numeric_limits<double>::infinity();
}


void TSPInstance::init_distance_matrix() {
    dist_matrix_.resize(dimension_);
    for (auto &row : dist_matrix_) {
        row.resize(dimension_);
    }
    for (auto i = 0u; i < dimension_; ++i) {
        auto &row = dist_matrix_[i];
        for (auto j = 0u; j < dimension_; ++j) {
            if (i < j) {
                row[j] = calc_distance(coords_.at(i), coords_.at(j));
                dist_matrix_[j][i] = row[j];
            }
        }
    }
    use_dist_matrix_ = true;
}


shared_ptr<TSPInstance> TSPInstance::load_from_tsplib_file(string path) {
    ifstream file(path);

    map<string, string> desc;
    vector<Point2D> coords;
    shared_ptr<TSPInstance> result = nullptr;

    if (file.is_open()) {
        string line;
        bool parse_points = false;
        while (getline(file, line)) {
            //transform(line.begin(), line.end(), line.begin(), ::tolower);

            if (line.find("EOF") != string::npos) {
                parse_points = false;
            } else if (parse_points) {
                istringstream in(line);
                int _;
                double x, y;
                in >> _ >> x >> y;
                coords.push_back( make_pair(x, y) );
            } else if (line.find("NODE_COORD_SECTION") != string::npos) {
                parse_points = true;
            } else if (line.find(":") != string::npos) {
                istringstream in(line);
                string key, value, token;
                in >> key;
                while (in >> token) {
                    if (token != ":") {
                        value += token;
                        value += " ";
                    }
                }
                trim(key);
                trim(value);
                if (key.back() == ':') {
                    key = key.substr(0, key.size() - 1);
                }
                desc[key] = value;
                cout << "key: [" << key << "] value: [" << value << "]" << endl;
            }
        }
        file.close();
    }
    auto name = desc["NAME"];
    auto weight_type_name = desc["EDGE_WEIGHT_TYPE"];
    if (weight_type_name == "EUC_2D") {
        result.reset(new TSPInstance(coords, EdgeWeightType::EUC_2D, name));
    } else if (weight_type_name == "GEO") {
        result.reset(new TSPInstance(coords, EdgeWeightType::GEO, name));
    } else {
        throw runtime_error("Unknown edge weight type: " +
                                 weight_type_name);
    }
    return result;
}


double TSPInstance::calc_relative_error(double sol_value) const {
    if (optimum_value_ > 0.0) {
        return (sol_value - optimum_value_) / optimum_value_;
    } else {
        return std::numeric_limits<double>::infinity();
    }
}
