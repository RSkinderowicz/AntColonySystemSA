
#include "log.h"
#include "json.hpp"

#include <vector>
#include <memory>
#include <cassert>
#include <iostream>

using namespace std;

// for convenience
using json = nlohmann::json;


vector< pair<string, shared_ptr<json>> > json_stack;


json &get_tail() {
    assert( !json_stack.empty() );
    return *(json_stack.back().second);
}


std::string log_get_contents() {
    if (json_stack.empty()) {
        return string("");
    }
    return get_tail().dump(2);
}


/*
 * Inits logging mechanisms - should be called only once.
 */
void log_init() {
    assert( json_stack.empty() );
    shared_ptr<json> dict(new json(json::object()));
    json_stack.push_back(make_pair("root", dict));
}


void log_open_dict(const std::string &key) {
    shared_ptr<json> dict(new json(json::object()));
    json_stack.push_back(make_pair(key, dict));
}


/**
 * Crates a new dictionary and appends it to the active (open) list.
 */
void log_open_dict() {
    log_open_dict("");
}


void log_close_dict() {
    assert( !json_stack.empty() );
    auto last = json_stack.back();
    json_stack.pop_back();
    if (!json_stack.empty()) {
        auto &tail = get_tail();

        if (last.first.empty()) {
            tail.push_back(*last.second);
        } else {
            tail[last.first] = *last.second;
        }
    }
}


void log_open_list(const std::string &key) {
    shared_ptr<json> dict(new json(json::array()));
    json_stack.push_back(make_pair(key, dict));
}


void log_open_list() {
    log_open_list("");
}


void log_close_list() {
    log_close_dict();
}


void log_add(const std::string &key, int value) {
    get_tail()[key] = value;
}


void log_add(const std::string &key, uint32_t value) {
    get_tail()[key] = value;
}


void log_add(const std::string &key, uint64_t value) {
    get_tail()[key] = value;
}


void log_add(const std::string &key, double value) {
    get_tail()[key] = value;
}


void log_add(const std::string &key, const std::string & value) {
    get_tail()[key] = value;
}


void log_add(const std::string &key, const std::vector<double> &vec) {
    get_tail()[key] = vec;
}


void log_add(const std::string &key, const std::vector<int> &vec) {
    get_tail()[key] = vec;
}


void log_add(const std::string &key, const std::vector<uint32_t> &vec) {
    get_tail()[key] = vec;
}


void log_add(int value) {
    get_tail().push_back(value);
}


void log_add(uint32_t value) {
    get_tail().push_back(value);
}


void log_add(double value) {
    get_tail().push_back(value);
}


void log_add(const std::string & value) {
    get_tail().push_back(value);
}


void log_add(const std::vector<double> &vec) {
    get_tail().push_back(vec);
}


void log_add(const std::vector<int> &vec) {
    get_tail().push_back(vec);
}


void log_add(const std::vector<uint32_t> &vec) {
    get_tail().push_back(vec);
}

