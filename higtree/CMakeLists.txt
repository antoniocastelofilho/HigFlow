CMAKE_minimum_required(VERSION 3.12)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
project(hig${dim}d)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_VIENNACL=1 -DUSE_SOR=1 -fopenmp")

add_library(${PROJECT_NAME} STATIC   
    ${HIGTREE_SRC_DIR}/allocator.c
    ${HIGTREE_SRC_DIR}/build-fringe.cpp
    ${HIGTREE_SRC_DIR}/coord.c
    ${HIGTREE_SRC_DIR}/domain.c
    ${HIGTREE_SRC_DIR}/global-stencil.c
    ${HIGTREE_SRC_DIR}/higtree.c 
    ${HIGTREE_SRC_DIR}/higtree-io.c
    ${HIGTREE_SRC_DIR}/higtree-iterator.c
    ${HIGTREE_SRC_DIR}/higtree-parallel.c 
    ${HIGTREE_SRC_DIR}/higtree-serialize.c
    ${HIGTREE_SRC_DIR}/isend-pool.c
    ${HIGTREE_SRC_DIR}/lbal.c 
    ${HIGTREE_SRC_DIR}/mapper.c
    ${HIGTREE_SRC_DIR}/mapper-syncer.c 
    ${HIGTREE_SRC_DIR}/pdomain.c
    ${HIGTREE_SRC_DIR}/point-cloud.c
    ${HIGTREE_SRC_DIR}/point-mapper.c
    ${HIGTREE_SRC_DIR}/rect.c
    ${HIGTREE_SRC_DIR}/rng.cpp
    ${HIGTREE_SRC_DIR}/solver.c 
    ${HIGTREE_SRC_DIR}/solver.camila.c
    ${HIGTREE_SRC_DIR}/solver-centralized.c 
    ${HIGTREE_SRC_DIR}/solver-debug-write.c
    ${HIGTREE_SRC_DIR}/solver-hypre.c
    ${HIGTREE_SRC_DIR}/solver-petsc.c
    ${HIGTREE_SRC_DIR}/solver-petsc.camila.c
    ${HIGTREE_SRC_DIR}/solver-sor.c
    ${HIGTREE_SRC_DIR}/solver-viennacl.cpp
    ${HIGTREE_SRC_DIR}/term-det.cpp
    ${HIGTREE_SRC_DIR}/utils.c
    ${HIGTREE_SRC_DIR}/wls.c
    ${HIGTREE_SRC_DIR}/uniqueid.c
)

target_include_directories(${PROJECT_NAME} PRIVATE ${HIGTREE_SRC_DIR})

target_link_libraries(${PROJECT_NAME} ${HIGFLOW_DEPENDENCIES})