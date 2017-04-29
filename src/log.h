#ifndef LOG_H
#define LOG_H 1

/*
 * Utility functions to log algorithms' progress data.
 */

#include <string>
#include <vector>


/**
 * Returns the whole log's contents as a single string.
 * The returned string is the log's data serialized to JSON format.
 */
std::string log_get_contents();

/*
 * Inits logging mechanism - should be called only once.
 */
void log_init();


/**
 * Adds a dictionary to the previously open dictionary.
 * The opened dict becomes the active one, i.e. subsequent add_* calls will
 * refer implicitly to this object.
 */
void log_open_dict(const std::string &key);


/**
 * Crates a new dictionary and appends it to the active (open) list.
 */
void log_open_dict();


/*
 * Closes previously open dictionary. Its parent becomes the current, active
 * one.
 */
void log_close_dict();


void log_open_list(const std::string &key);


void log_open_list();


void log_close_list();


void log_add(const std::string &key, int value);

void log_add(const std::string &key, uint32_t value);

void log_add(const std::string &key, uint64_t value);

void log_add(const std::string &key, double value);

void log_add(const std::string &key, const std::string & value);

void log_add(const std::string &key, const std::vector<double> &vec);

void log_add(const std::string &key, const std::vector<int> &vec);

void log_add(const std::string &key, const std::vector<uint32_t> &vec);

void log_add(int value);

void log_add(uint32_t value);

void log_add(double value);

void log_add(const std::string & value);

void log_add(const std::vector<double> &vec);

void log_add(const std::vector<int> &vec);

void log_add(const std::vector<uint32_t> &vec);

#endif /* ifndef LOG_H 1 */
