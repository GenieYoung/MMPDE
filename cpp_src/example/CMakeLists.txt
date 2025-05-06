file(GLOB EXAMPLE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

foreach(example ${EXAMPLE_LIST})
    get_filename_component(exec ${example} NAME_WE)
    add_executable(${exec} ${example})
    target_link_libraries(${exec} PUBLIC MMPDE_LIB)
endforeach()
