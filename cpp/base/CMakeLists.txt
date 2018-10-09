cmake_minimum_required(VERSION 3.9.0)

project(ising CXX)

set(SOURCE
  ising.cpp)

add_executable(ising ${SOURCE})
set_target_properties(ising PROPERTIES CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON)

target_compile_options(ising PUBLIC -Wall -Wextra)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(ising PUBLIC -Wno-c++98-compat)
endif()

target_link_libraries(ising stdc++fs)
