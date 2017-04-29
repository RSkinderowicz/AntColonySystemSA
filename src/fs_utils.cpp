#include "fs_utils.h"

#include <unistd.h>
#include <sys/stat.h>
#include <algorithm>
#include <stdexcept>



bool is_path_relative(std::string path) {
    return path.size() > 0 && path[0] != '/';
}



bool file_exists(std::string path) {
    return access(path.c_str(), F_OK) != -1;
}


bool dir_exists(std::string path) {
    struct stat st = {0};
    return stat(path.c_str(), &st) != -1;
}


bool crate_dir(std::string path) {
    return mkdir(path.c_str(), 0777) == 0;
}


void create_dir_path(std::string path) {
    auto it = path.begin();
    it = find(it + 1, path.end(), '/');
    while (it != path.end()) {
        auto subpath = std::string(path.begin(), it);
        if (!dir_exists(subpath)) {
            if (!crate_dir(subpath)) {
                throw std::runtime_error(std::string("Cannot create dir: ") + subpath);
            }
        }
        it = find(it+1, path.end(), '/');
    }
    if (!dir_exists(path)) {
        if (!crate_dir(path)) {
            throw std::runtime_error(std::string("Cannot create dir: ") + path);
        }
    }
}



std::string get_cwd() {
    const auto size = 1024u;
    char buf[size];
    auto *res = getcwd(buf, size);
    if (res != nullptr) {
        return std::string(res);
    }
    return std::string();
}
