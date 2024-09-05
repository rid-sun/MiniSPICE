#
file(GLOB srcs ./*.cpp)
add_library(utils STATIC ${srcs})
