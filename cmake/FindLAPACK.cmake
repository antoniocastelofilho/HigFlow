# - Try to find blas
# Once done this will define
#  LAPACK_FOUND - System has LAPACK
#  LAPACK_LIBRARIES - The libraries needed to use LAPACK
#  LAPACK_DEFINITIONS - Compiler switches required for using LAPACK


set(LAPACK_LMOD 0)
if(DEFINED ENV{LAPACK_LIB_DIR})
    find_library(LAPACK_LIBRARY lapack HINTS "$ENV{LAPACK_LIB_DIR}") 
    set(LAPACK_LMOD 1)
endif()

set(LAPACK_DOCKER 0)
if(DEFINED ENV{LAPACK_LIB})
    find_library(LAPACK_LIBRARY lapack PATHS "$ENV{LAPACK_LIB}") 
    set(LAPACK_DOCKER 1)
endif()

set(LAPACK_ENV 0)
if(DEFINED ENV{LAPACK})
    set(LAPACK_LIBRARY "$ENV{LAPACK}") 
    set(LAPACK_ENV 1)
endif()

if(NOT ${LAPACK_LMOD} AND NOT ${LAPACK_DOCKER} AND NOT ${LAPACK_ENV})
    set(LAPACK_PREFIX "${LAPACK_PREFIX_DEFAULT}" CACHE STRING "LAPACK install directory")
    if(LAPACK_PREFIX)
        message(STATUS "LAPACK_PREFIX ${LAPACK_PREFIX}")
    endif()

    find_library(LAPACK_LIBRARY lapack PATHS "${LAPACK_PREFIX}/lib")  
endif()

set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LAPACK_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    LAPACK
    DEFAULT_MSG
    LAPACK_LIBRARY
)

mark_as_advanced(LAPACK_LIBRARY)