#!/bin/bash

# export OMP_THREAD_NUM=1, before usage
# source /opt/intel/oneapi/setvars.sh
rm -fr build; mkdir build; cd build;
export PARMETIS_ROOT="/root/local"
cmake .. \
    -DTPL_PARMETIS_INCLUDE_DIRS="${PARMETIS_ROOT}/include" \
    -DTPL_PARMETIS_LIBRARIES="${PARMETIS_ROOT}/lib/libparmetis.a;${PARMETIS_ROOT}/lib/libmetis.a;${PARMETIS_ROOT}/lib/libGKlib.a" \
    -DTPL_ENABLE_INTERNAL_BLASLIB=OFF \
    -DTPL_ENABLE_COMBBLASLIB=OFF\
    -DCMAKE_C_FLAGS="-std=c99 -O3 -g -DPRNTlevel=0 -DDEBUGlevel=0" \
    -DCMAKE_C_COMPILER=mpiicx \
    -DCMAKE_CXX_COMPILER=mpiicpx \
    -DCMAKE_CXX_FLAGS="-std=c++11" \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DBLA_VENDOR=Intel10_64lp \
    -DCMAKE_LINKER=mpiicpx \
    -Denable_openmp=ON \
    -DTPL_ENABLE_LAPACKLIB=ON \
    -DCMAKE_INSTALL_PREFIX=..

make install -j8
