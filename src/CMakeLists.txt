# 注册子级别目录
add_subdirectory(solver)

# 注册子级别目录
add_subdirectory(plot)

# 注册子级别目录
add_subdirectory(utils)

# Analysis
add_subdirectory(analysis)

# 寻找所有的src文件夹下的cpp文件
file(GLOB srcs ./*.cpp)
add_executable(minispice ${srcs})

target_link_libraries(minispice PRIVATE analysis plot utils solver)

# 指定parser输出路径
set_target_properties(minispice PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
