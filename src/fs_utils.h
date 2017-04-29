#ifndef FS_UTILS_H
#define FS_UTILS_H

/**
 * A few filesystem related auxiliary functions.
 */

#include <string>


bool is_path_relative(std::string path);


bool file_exists(std::string path);


bool dir_exists(std::string path);


bool crate_dir(std::string path);


/**
 * Simliar to the "mkdir -p this/is/a/path" linux command - creates a nested
 * directories if they do not exist.
 *
 * Throws std::runtime_error if failed.
 */
void create_dir_path(std::string path);


/**
 * Returns current working directory or an empty string if an error has occured.
 */
std::string get_cwd();

#endif
