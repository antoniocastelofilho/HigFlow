CMAKE_minimum_required(VERSION 3.12)

########################################
# compiler config
########################################
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wall -Wno-unused-result -Wno-unused-variable -DDIM=${dim}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -Wall -Wno-unused-result -Wno-unused-variable -DDIM=${dim}")

if(${debug})
    add_compile_options(-DDEBUG=1 -g -ggdb -Og)
    add_link_options(-DDEBUG=1 -g -ggdb -Og)
else()
    add_compile_options(-Ofast -DNDEBUG=1)
    add_link_options(-Ofast -DNDEBUG=1 -march=native -mtune=native -flto=8 -fno-fat-lto-objects -fno-strict-aliasing)
endif()

########################################
# dependencies 
########################################
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)

include_directories(${MPI_INCLUDE_PATH} ${HDF5_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS})
set(HIGFLOW_DEPENDENCIES ${HDF5_LIBRARIES} ${MPI_C_LIBRARIES})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# blas/lapack
find_package(BLAS)
list(APPEND MFSIM_DEPENDENCIES ${BLAS_LIBRARIES})

find_package(LAPACK)
list(APPEND MFSIM_DEPENDENCIES ${LAPACK_LIBRARIES})

# loadbalance
find_package(Zoltan)
include_directories(${ZOLTAN_INCLUDE_DIRS})
list(APPEND HIGFLOW_DEPENDENCIES ${ZOLTAN_LIBRARIES})

# petsc
find_package(PETSc)
include_directories(${PETSC_INCLUDE_DIRS})
list(APPEND HIGFLOW_DEPENDENCIES ${PETSC_LIBRARIES})

# threads and boost
find_package(Threads REQUIRED)
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND HIGFLOW_DEPENDENCIES ${Boost_LIBRARIES} Threads::Threads)

# openmp
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# viennacl
find_package(ViennaCL)
include_directories(${VIENNACL_INCLUDE_DIRS})
list(APPEND HIGFLOW_DEPENDENCIES ${VIENNACL_LIBRARIES})

# glib
FIND_PACKAGE(PkgConfig REQUIRED)

PKG_CHECK_MODULES(GLIB REQUIRED glib-2.0)

include_directories(${GLIB_INCLUDE_DIRS})
list(APPEND HIGFLOW_DEPENDENCIES ${GLIB_LIBRARIES})

# libfyaml

PKG_CHECK_MODULES(LIBFYAML REQUIRED libfyaml)

include_directories(${LIBFYAML_INCLUDE_DIRS})
list(APPEND HIGFLOW_DEPENDENCIES ${LIBFYAML_LIBRARIES})


PKG_CHECK_MODULES(NUMACTL REQUIRED numa)

include_directories(${NUMACTL_INCLUDE_DIRS})
list(APPEND HIGFLOW_DEPENDENCIES ${NUMACTL_LIBRARIES})

########################################
# paths
########################################

set(HIGTREE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higtree)
set(HIGTREE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higtree/src)

set(HIGFLOW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow)
set(HIGFLOW_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/src)

set(HIGFLOW_BMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_BMP)
set(HIGFLOW_GPTT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_Gptt)
set(HIGFLOW_KBKZ_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_KBKZ)
set(HIGFLOW_NEWT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_Newt)
set(HIGFLOW_OLDROYD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_Oldroyd)
set(HIGFLOW_VOF_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_VOF)
set(HIGFLOW_VOF_GPTT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_VOF_Gptt)
set(HIGFLOW_VOF_OLDROYD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/higflow/example2d_VOF_Oldroyd)

set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

########################################
# modules
########################################

# higtree
include(${HIGTREE_DIR}/CMakeLists.txt)

# higflow
include(${HIGFLOW_DIR}/CMakeLists.txt)

########################################
# examples
########################################

# example_BMP
include(${HIGFLOW_BMP_DIR}/CMakeLists.txt)

# example_Gptt
include(${HIGFLOW_GPTT_DIR}/CMakeLists.txt)

# example_KBKZ
include(${HIGFLOW_KBKZ_DIR}/CMakeLists.txt)

# example_Newt
include(${HIGFLOW_NEWT_DIR}/CMakeLists.txt)

# example_Oldroyd
include(${HIGFLOW_OLDROYD_DIR}/CMakeLists.txt)

# example_VOF
include(${HIGFLOW_VOF_DIR}/CMakeLists.txt)

# example_VOF_Gptt
include(${HIGFLOW_VOF_GPTT_DIR}/CMakeLists.txt)

# example_VOF_Oldroyd
include(${HIGFLOW_VOF_OLDROYD_DIR}/CMakeLists.txt)

########################################
# Install
########################################

if(NOT ${prefix} STREQUAL "")
    set(INSTALL_PATH ${prefix})
elseif(NOT ${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local")
    set(INSTALL_PATH ${CMAKE_INSTALL_PREFIX})
else()
    set(INSTALL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/install)
endif()

# includes
install(DIRECTORY ${HIGTREE_SRC_DIR}/ DESTINATION ${INSTALL_PATH}/include/higtree FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${HIGTREE_SRC_DIR}/ DESTINATION ${INSTALL_PATH}/include/higtree FILES_MATCHING PATTERN "*.hpp")

install(DIRECTORY ${HIGFLOW_SRC_DIR}/ DESTINATION ${INSTALL_PATH}/include/higflow FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${HIGFLOW_SRC_DIR}/ DESTINATION ${INSTALL_PATH}/include/higflow FILES_MATCHING PATTERN "*.hpp")

# examples 
install(FILES ${BUILD_DIR}/example_BMP PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_BMP)
install(DIRECTORY ${HIGFLOW_BMP_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_BMP/input)
install(DIRECTORY ${HIGFLOW_BMP_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_BMP/amrs)

install(FILES ${BUILD_DIR}/example_Gptt PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_Gptt)
install(DIRECTORY ${HIGFLOW_GPTT_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_Gptt/input)
install(DIRECTORY ${HIGFLOW_GPTT_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_Gptt/amrs)

install(FILES ${BUILD_DIR}/example_KBKZ PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_KBKZ)
install(DIRECTORY ${HIGFLOW_KBKZ_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_KBKZ/input)
install(DIRECTORY ${HIGFLOW_KBKZ_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_KBKZ/amrs)

install(FILES ${BUILD_DIR}/example_Newt PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_Newt)
install(DIRECTORY ${HIGFLOW_NEWT_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_Newt/input)
install(DIRECTORY ${HIGFLOW_NEWT_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_Newt/amrs)

install(FILES ${BUILD_DIR}/example_Oldroyd PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_Oldroyd)
install(DIRECTORY ${HIGFLOW_OLDROYD_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_Oldroyd/input)
install(DIRECTORY ${HIGFLOW_OLDROYD_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_Oldroyd/amrs)

install(FILES ${BUILD_DIR}/example_VOF PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_VOF)
install(DIRECTORY ${HIGFLOW_VOF_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_VOF/input)
install(DIRECTORY ${HIGFLOW_VOF_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_VOF/amrs)

install(FILES ${BUILD_DIR}/example_VOF_Gptt PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_VOF_Gptt)
install(DIRECTORY ${HIGFLOW_VOF_GPTT_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_VOF_Gptt/input)
install(DIRECTORY ${HIGFLOW_VOF_GPTT_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_VOF_Gptt/amrs)

install(FILES ${BUILD_DIR}/example_VOF_Oldroyd PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/example_VOF_Oldroyd)
install(DIRECTORY ${HIGFLOW_VOF_OLDROYD_DIR}/input/ DESTINATION ${INSTALL_PATH}/example_VOF_Oldroyd/input)
install(DIRECTORY ${HIGFLOW_VOF_OLDROYD_DIR}/amrs/ DESTINATION ${INSTALL_PATH}/example_VOF_Oldroyd/amrs)


install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/etc/higflow-env.sh PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/etc)

install(FILES ${BUILD_DIR}/libhig${dim}d.a PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/lib)
install(FILES ${BUILD_DIR}/libhigflow${dim}d.a PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION ${INSTALL_PATH}/lib)