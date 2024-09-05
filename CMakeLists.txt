cmake_minimum_required(VERSION 3.22.0)
project(MiniSPICE VERSION 1.0)

# options
option(USE_EIGEN "USE_EIGEN" OFF)
option(USE_KLU "USE_KLU" OFF)
option(USE_SUPERLU "USE_SUPERLU" OFF)
option(USE_MUMPS "USE_MUMPS" OFF)
option(USE_MINISPICE "USE_MINISPICE" OFF)

# macro num
set(ENABLED_COUNT 0)
if(USE_EIGEN)
    math(EXPR ENABLED_COUNT "${ENABLED_COUNT} + 1")
endif()
if(USE_KLU)
    math(EXPR ENABLED_COUNT "${ENABLED_COUNT} + 1")
endif()
if(USE_SUPERLU)
    math(EXPR ENABLED_COUNT "${ENABLED_COUNT} + 1")
endif()
if(USE_MUMPS)
    math(EXPR ENABLED_COUNT "${ENABLED_COUNT} + 1")
endif()
if(USE_MINISPICE)
    math(EXPR ENABLED_COUNT "${ENABLED_COUNT} + 1")
endif()

#
if(ENABLED_COUNT GREATER 1)
    message(FATAL_ERROR "Cannot enable more than one of USE_EIGEN, USE_KLU, USE_MUMPS, USE_SUPERLU and USE_MINISPICE simultaneously.")
endif()

# 链接工具库-解法器
if(ENABLED_COUNT EQUAL 0)
    add_definitions(-DMINISPICEKLU)
    message("The solver is not defined; defaulting to KLU.")
    include_directories(${CMAKE_SOURCE_DIR}/3rd_lib/SuiteSparse/include/suitesparse)
elseif(USE_EIGEN)
    add_definitions(-DMINISPICEEIGEN)
    include_directories(${CMAKE_SOURCE_DIR}/3rd_lib/eigen-3.4.0)
    message("The solver is EIGEN")
elseif(USE_KLU)
    add_definitions(-DMINISPICEKLU)
    message("The solver is KLU")
    include_directories(${CMAKE_SOURCE_DIR}/3rd_lib/SuiteSparse/include/suitesparse)
elseif(USE_SUPERLU)
    add_definitions(-DMINISPICESUPERLU)
    message("The solver is SUPERLU")
    include_directories(${CMAKE_SOURCE_DIR}/3rd_lib/superlu_dist/include)
elseif(USE_MUMPS)
    add_definitions(-DMINISPICEMUMPS)
    message("The solver is MUMPS")
elseif(USE_MINISPICE)
    add_definitions(-DMINSPICEMINISPICELU)
    message("The solver is MINISPICELU")
endif()

# 设置编译器
find_package(MPI REQUIRED)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    link_libraries(${MPI_LIBRARIES})
else()
    message(FATAL_ERROR "MPI not found. Please source oneAPI setvars.sh.")
endif()

# set(CMAKE_C_COMPILER "mpiicc")
# set(CMAKE_CXX_COMPILER "mpiicx")
# set(CMAKE_C_COMPILER "mpicc")
# set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")

# 添加项目头文件
include_directories(${CMAKE_SOURCE_DIR}/include)

# 创建out目录并复制输出
set(out ${CMAKE_SOURCE_DIR}/bin)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${out})

# 注册子级别目录
add_subdirectory(src)
