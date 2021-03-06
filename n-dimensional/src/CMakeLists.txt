# all sources except main
set(SOURCE
  montecarlo.cpp
  lattice.cpp
  fileio.cpp)

# store sources for other modules
set(isingsrc)
foreach (src IN LISTS SOURCE)
  list(APPEND isingsrc "${CMAKE_CURRENT_LIST_DIR}/${src}")
endforeach ()
set(BASE_SOURCE ${isingsrc} PARENT_SCOPE)
unset(isingsrc)

add_executable(ising ${SOURCE} main.cpp)
set_target_properties(ising PROPERTIES CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON)
target_link_libraries(ising stdc++fs)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(ising PUBLIC ${GCC_CLANG_WARNINGS})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(ising PUBLIC ${GCC_CLANG_WARNINGS} ${GCC_EXTRA_WARNINGS})
endif ()

find_package(yaml-cpp REQUIRED)
target_include_directories(ising PUBLIC ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(ising ${YAML_CPP_LIBRARIES})
