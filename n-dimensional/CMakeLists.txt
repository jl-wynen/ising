cmake_minimum_required(VERSION 3.9.0)

project(ising CXX)

add_subdirectory(src)

add_subdirectory(test)
add_test(NAME IsingTest COMMAND ising-test)
enable_testing()
