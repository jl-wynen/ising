cmake_minimum_required(VERSION 3.9.0)

project(ising-dim1 CXX)

set(SOURCE
  ising1d.cpp)

add_executable(ising1d ${SOURCE})
set_target_properties(ising1d PROPERTIES CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)

target_compile_options(ising1d PUBLIC -Wall -Wextra)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(ising1d PUBLIC -Wno-c++98-compat)
endif()

target_link_libraries(ising1d stdc++fs)
