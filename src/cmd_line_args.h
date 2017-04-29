#ifndef CMD_LINE_ARGS_H
#define CMD_LINE_ARGS_H

#include <string>
#include <map>

/*
 * Holds a string value and allows its conversion to a specific type via
 * overloaded type conversion operator methods.
 */
class StringValue {
public:

    StringValue(const std::string &value) :
        value_(value) {
    }

    operator int() { return std::stoi(value_); }

    operator unsigned int() { return (unsigned int)std::stoul(value_); }

    operator long() { return std::stol(value_); }

    operator long long() { return std::stoll(value_); }

    operator float() { return std::stof(value_); }

    operator double() { return std::stod(value_); }

    operator long double() { return std::stold(value_); }

    operator std::string() { return value_; }

private:
    const std::string &value_;
};


void init_cmd_line_args(const std::string &usage, int argc, char *argv[]);


void print_cmd_line_args();


std::map<std::string, std::string> get_cmd_line_args();


StringValue get_cmd_line_arg(const std::string &name);


#endif /* ifndef CMD_LINE_ARGS_H */
