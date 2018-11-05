set(TEST_SOURCE
  lattice.cpp
  rng.cpp
  ising.cpp
  fileio.cpp
  test.cpp)

add_executable(ising-test EXCLUDE_FROM_ALL ${TEST_SOURCE} ${BASE_SOURCE})
set_target_properties(ising-test PROPERTIES CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON)
target_compile_options(ising-test PUBLIC -Wall -Wextra)

target_include_directories(ising-test PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(ising-test stdc++fs)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
    OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(ising-test PUBLIC -fsanitize=address -fno-omit-frame-pointer)
endif()
set_target_properties(ising-test PROPERTIES LINK_FLAGS "-fsanitize=address -fno-omit-frame-pointer")

target_compile_definitions(ising-test PUBLIC "-DISLE_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}")

find_package(yaml-cpp REQUIRED)
target_include_directories(ising-test PUBLIC ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(ising-test ${YAML_CPP_LIBRARIES})
