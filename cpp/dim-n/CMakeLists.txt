cmake_minimum_required(VERSION 3.9.0)

project(ising CXX)

# all sources except main
set(SOURCE
  montecarlo.cpp
  lattice.cpp
  fileio.cpp)

# store sources for test
set(libsrc)
foreach (src IN LISTS SOURCE)
  list(APPEND libsrc "${CMAKE_CURRENT_LIST_DIR}/${src}")
endforeach ()
set(BASE_SOURCE ${libsrc} CACHE INTERNAL "Source files of the base program" FORCE)
unset(libsrc)

add_executable(ising ${SOURCE} main.cpp)
set_target_properties(ising PROPERTIES CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON)
target_compile_options(ising PUBLIC -Wall -Wextra -Wpedantic)
target_link_libraries(ising stdc++fs)

add_subdirectory(test)
add_test(NAME IsingTest COMMAND ising-test)
enable_testing()

find_package(yaml-cpp REQUIRED)
target_include_directories(ising PUBLIC ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(ising ${YAML_CPP_LIBRARIES})
