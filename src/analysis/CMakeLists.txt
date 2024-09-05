#
file(GLOB srcs ./*.cpp)
add_library(analysis STATIC ${srcs})
