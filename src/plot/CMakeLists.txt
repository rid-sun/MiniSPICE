# 引入 FindPython3 模块
find_package(Python3 REQUIRED)

message(STATUS "Python3 interpreter: ${Python3_EXECUTABLE}")
message(STATUS "Python3 include dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 library: ${Python3_LIBRARIES}")

include_directories(${Python3_INCLUDE_DIRS})

file(GLOB srcs ./*.cpp)
add_library(plot STATIC ${srcs} ${Python3_LIBRARIES})
