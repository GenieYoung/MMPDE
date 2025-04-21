cmake_minimum_required(VERSION 3.10...3.29)

project(MMPDE)

set(CMAKE_BUILD_TYPE Debug)

add_subdirectory(src)

add_subdirectory(example)

add_executable(main ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(main PUBLIC MMPDE_LIB)
