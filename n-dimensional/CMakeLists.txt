cmake_minimum_required(VERSION 3.9.0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(GCC_CLANG_WARNINGS -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wpedantic
 -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
 -Wconversion -Wsign-conversion -Wnull-dereference
 -Wdouble-promotion -Wformat=2)

set(GCC_EXTRA_WARNINGS -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches
 -Wlogical-op -Wuseless-cast)

project(ising CXX)

add_subdirectory(src)

add_subdirectory(test)
add_test(NAME IsingTest COMMAND ising-test)
enable_testing()
