#!/bin/bash

# ignore when building with dokcer
# source /opt/intel/oneapi/setvars.sh

# build minispice
rm -rf build
mkdir build && cd build
cmake -DUSE_SUPERLU=ON -DCMAKE_C_COMPILER="$(which mpiicx)" -DCMAKE_CXX_COMPILER="$(which mpiicpx)" ..
make
