#include "config.hpp"
#include <cstring>

namespace parser_cmd {
    std::string inFileName = "";
    std::string outFileName = "";
    bool help = false;
    bool instrError = false;

    void parse_args(int argc, char** argv) {
        for (int i = 1; i < argc;i++) {
            if (argv[i][0] == '-') {
                if (strcmp(argv[i], "-o") == 0) {
                    if (strcmp(argv[i++], "-") == 0)
                        outFileName = ""; 
                    else
                        outFileName = argv[i];
                } else if (strcmp(argv[i], "-f") == 0) {
                    if (strcmp(argv[i++], "-") == 0)
                        inFileName = "";
                    else
                        inFileName = argv[i];
                } else if (strcmp(argv[i], "-h") == 0) {
                    help = true;
                    break;
                } else {
                    instrError = true;
                }
            } else {
                instrError = true;
            }
        }
    }
}