# 引入 FindPython3 模块
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

message(STATUS "Python3 interpreter: ${Python3_EXECUTABLE}")
message(STATUS "Python3 include dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 library: ${Python3_LIBRARIES}")

include_directories(${Python3_INCLUDE_DIRS})

# 配置脚本路径
set(PLOT_SCRIPT_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
add_compile_definitions(PLOT_SCRIPT_PATH="${PLOT_SCRIPT_PATH}")

file(GLOB srcs ./*.cpp)
# add_library(plot STATIC ${srcs} ${Python3_LIBRARIES})
add_library(plot STATIC ${srcs})
target_link_libraries(plot PRIVATE ${Python3_LIBRARIES})
