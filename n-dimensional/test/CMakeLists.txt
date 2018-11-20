set(TEST_SOURCE
  lattice.cpp
  rng.cpp
  ising.cpp
  fileio.cpp
  test.cpp)

add_executable(ising-test EXCLUDE_FROM_ALL ${TEST_SOURCE} ${BASE_SOURCE})
set_target_properties(ising-test PROPERTIES CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(ising PUBLIC ${GCC_CLANG_WARNINGS})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(ising PUBLIC ${GCC_CLANG_WARNINGS} ${GCC_EXTRA_WARNINGS})
endif ()

target_include_directories(ising-test PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(ising-test stdc++fs)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
    OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(ising-test PUBLIC -fsanitize=address -fno-omit-frame-pointer)
  set_target_properties(ising-test PROPERTIES LINK_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
endif()

target_compile_definitions(ising-test PUBLIC "-DISLE_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}")

find_package(yaml-cpp REQUIRED)
target_include_directories(ising-test PUBLIC ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(ising-test ${YAML_CPP_LIBRARIES})
