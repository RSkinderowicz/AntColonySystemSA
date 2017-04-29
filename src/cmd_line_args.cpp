#include <map>
#include <iostream>
#include <sstream>

#include "cmd_line_args.h"
#include "docopt.h"


/*
 * A global map with the command line arguments.
 * It is populated by init_cmd_line_args function.
 */
std::map<std::string, docopt::value> global_cmd_line_args;


/**
 * Parses the command line arguments passed to main function.
 * The arguments are inserted into 'global_cmd_line_args' map.
 */
void init_cmd_line_args(const std::string &usage, int argc, char *argv[]) {
    const char * help_argv[] = { argv[0], "-h" };
    if (argc == 1) {
        argv = const_cast<char**>(help_argv);
        argc = 2;
    }
    global_cmd_line_args = docopt::docopt(usage.c_str(),
                                          { argv + 1, argv + argc },
                                          true,               // show help if requested
                                          "ACO framework 0.1");  // version string
}


void print_cmd_line_args() {
    std::cout << "Command line arguments:" << std::endl;
    for(auto const& arg : global_cmd_line_args) {
        std::cout << '\t' << arg.first << ": " << arg.second << std::endl;
    }
}


std::map<std::string, std::string> get_cmd_line_args() {
    std::map<std::string, std::string> res;
    for(auto const& arg : global_cmd_line_args) {
        std::ostringstream out;
        out << arg.second;
        const auto s = out.str();
        res[arg.first] = s.substr(1, s.size()-2);
    }
    return res;
}


StringValue get_cmd_line_arg(const std::string &name) {
    auto it = global_cmd_line_args.find(name);
    if (it != global_cmd_line_args.end()) {
        return StringValue(it->second.asString());
    }
    throw std::runtime_error("No command line argument named: " + name);
    return StringValue("0");
}
