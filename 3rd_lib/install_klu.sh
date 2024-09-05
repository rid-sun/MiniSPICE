#!/bin/bash

# source /opt/intel/oneapi/setvars.sh
rm -rf build
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=..\
      -DCMAKE_C_COMPILER=icx \
      -DCMAKE_CXX_COMPILER=icpx \
      -DSUITESPARSE_ENABLE_PROJECTS="amd;btf;colamd;klu;camd;ccolamd;cholmod" \
      -DBUILD_STATIC_LIBS=ON \
      -DBLA_VENDOR=Intel10_64lp \
      ..

make install -j8
