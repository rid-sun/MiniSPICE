
#ifndef MINISPICE_H
#define MINISPICE_H

#include <memory>
#include "device.hpp"
#include "MiniSPICEMatrix.hpp"
#include "mpi.h"

/****************    宏定义   **********************/

#define MINISPICEERROR NA
#define MINISPICEOK OK
#define NONCONV NA

#define MPI_ENV_0_PID_BEGIN(x) if (x == 0) {
#define MPI_ENV_0_PID_END }

/****************  外部类型声明  ***** **************/

/* parser.cpp */
int parser(std::unique_ptr<Netlist>&, std::unique_ptr<Analysis>&, std::string, std::string);

/* setup.cpp */
int setup(std::unique_ptr<SMPMatrix>&, std::unique_ptr<Netlist>&, std::unique_ptr<Analysis>&, std::string);

/* unsetup.cpp */
int unsetup(std::unique_ptr<SMPMatrix>&, std::unique_ptr<Netlist>&, std::unique_ptr<Analysis>&);

/* load.cpp */
int load(std::unique_ptr<SMPMatrix>&, std::unique_ptr<Netlist>&, std::unique_ptr<Analysis>&);

/* solve.cpp 条件编译的多态^_^ */
int solve(std::unique_ptr<SMPMatrix>&, std::unique_ptr<Netlist>&, std::unique_ptr<Analysis>&, int, int, int);

/* dojob.cpp  */
int dojob(std::unique_ptr<SMPMatrix>&, std::unique_ptr<Netlist>&, std::unique_ptr<Analysis>&, std::string, int);

/* plot.cpp */
int plot(std::unique_ptr<Analysis>& analysis);


#endif
