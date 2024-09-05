
#ifndef MINISPICE_CONFIG_H
#define MINISPICE_CONFIG_H

#include <string>

namespace parser_cmd {
    extern std::string inFileName;
    extern std::string outFileName;
    extern bool help;
    extern bool instrError;

    void parse_args(int argc, char** argv);
}

#endif
