#!/bin/bash

# ignore when building with dokcer
# source /opt/intel/oneapi/setvars.sh

# install needed 3rd_libs
## superlu_dist
git clone https://github.com/xiaoyeli/superlu_dist.git
cd superlu_dist
bash ../install_superlu.sh
cd ../
## KLU
git clone https://github.com/DrTimothyAldenDavis/SuiteSparse.git
cd SuiteSparse
bash ../install_klu.sh
cd ../
# ## Eigen [no-uasge here]
# wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
# tar -zxvf eigen-3.4.0.tar.gz
# cd ../
