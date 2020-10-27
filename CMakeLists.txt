cmake_minimum_required(VERSION 3.8)

project(cg_method)
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES eisenstat.c)
add_executable(simple_example ${SOURCE_FILES})
