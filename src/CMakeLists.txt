file(GLOB MMPDE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(MMPDE_LIB STATIC ${MMPDE_SRC})

set(CMTL_SRC /mnt/d/yjl/code/CMTL/src)

target_include_directories(MMPDE_LIB PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMTL_SRC})