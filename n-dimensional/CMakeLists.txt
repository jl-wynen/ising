cmake_minimum_required(VERSION 3.9.0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

project(ising CXX)

add_subdirectory(src)

add_subdirectory(test)
add_test(NAME IsingTest COMMAND ising-test)
enable_testing()
