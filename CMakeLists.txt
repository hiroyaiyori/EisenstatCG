cmake_minimum_required(VERSION 3.8)

project(cg_method)
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard

# Add main.cpp file of project root directory as source file
set(TAKAHASHI_SOURCE takahashi.c)
add_executable(takahashi ${TAKAHASHI_SOURCE})

set(EISENSTAT_SOURCE eisenstat.c)
add_executable(eisenstat ${EISENSTAT_SOURCE})
