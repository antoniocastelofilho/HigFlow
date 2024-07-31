# - Try to find blas
# Once done this will define
#  BLAS_FOUND - System has BLAS
#  BLAS_LIBRARIES - The libraries needed to use BLAS
#  BLAS_DEFINITIONS - Compiler switches required for using BLAS


set(BLAS_LMOD 0)
if(DEFINED ENV{BLAS_LIB_DIR})
    find_library(BLAS_LIBRARY blas HINTS "$ENV{BLAS_LIB_DIR}") 
    set(BLAS_LMOD 1)
endif()

set(BLAS_DOCKER 0)
if(DEFINED ENV{BLAS_LIB})
    find_library(BLAS_LIBRARY blas PATHS "$ENV{BLAS_LIB}") 
    set(BLAS_DOCKER 1)
endif()

set(BLAS_LAPACK_LMOD 0)
if(DEFINED ENV{LAPACK_LIB_DIRR})
    find_library(BLAS_LIBRARY blas HINTS "$ENV{LAPACK_LIB_DIR}") 
    set(BLAS_LAPACK_LMOD 1)
endif()

set(BLAS_LAPACK_DOCKER 0)
if(DEFINED ENV{LAPACK_LIB})
    find_library(BLAS_LIBRARY blas PATHS "$ENV{LAPACK_LIB}") 
    set(BLAS_LAPACK_DOCKER 1)
endif()

set(BLAS_ENV 0)
if(DEFINED ENV{BLAS})
    set(BLAS_LIBRARY "$ENV{BLAS}") 
    set(BLAS_ENV 1)
endif()

if(NOT ${BLAS_LMOD} AND NOT ${BLAS_DOCKER} AND NOT ${BLAS_LAPACK_LMOD} AND NOT ${BLAS_LAPACK_DOCKER} AND NOT ${BLAS_LAPACK_DOCKER} AND NOT ${BLAS_ENV})
    set(BLAS_PREFIX "${BLAS_PREFIX_DEFAULT}" CACHE STRING "BLAS install directory")
    if(BLAS_PREFIX)
        message(STATUS "BLAS_PREFIX ${BLAS_PREFIX}")
    endif()

    find_library(BLAS_LIBRARY blas PATHS "${BLAS_PREFIX}/lib")  
endif()

set(BLAS_LIBRARIES ${BLAS_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set BLAS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    BLAS
    DEFAULT_MSG
    BLAS_LIBRARY
)

mark_as_advanced(BLAS_LIBRARY)