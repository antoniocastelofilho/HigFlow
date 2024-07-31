# - Try to find zoltan
# Once done this will define
#  ZOLTAN_FOUND - System has ZOLTAN
#  ZOLTAN_INCLUDE_DIRS - The ZOLTAN include directories
#  ZOLTAN_LIBRARIES - The libraries needed to use ZOLTAN
#  ZOLTAN_DEFINITIONS - Compiler switches required for using ZOLTAN

set(ZOLTAN_LMOD 0)
if(DEFINED ENV{ZOLTAN_INC_DIR} AND DEFINED ENV{ZOLTAN_LIB_DIR})
    find_path(ZOLTAN_INCLUDE_DIR zoltan.h PATHS "$ENV{ZOLTAN_INC_DIR}")
    find_library(ZOLTAN_LIBRARY zoltan PATHS "$ENV{ZOLTAN_LIB_DIR}") 
    set(ZOLTAN_LMOD 1)
endif()

set(ZOLTAN_DOCKER 0)
if(DEFINED ENV{ZOLTAN_INC} AND DEFINED ENV{ZOLTAN_LIB})
    find_path(ZOLTAN_INCLUDE_DIR zoltan.h PATHS "$ENV{ZOLTAN_INC}")
    find_library(ZOLTAN_LIBRARY zoltan PATHS "$ENV{ZOLTAN_LIB}") 
    set(ZOLTAN_DOCKER 1)
endif()

if(NOT ${ZOLTAN_LMOD} AND NOT ${ZOLTAN_DOCKER})
    set(ZOLTAN_PREFIX "${ZOLTAN_PREFIX_DEFAULT}" CACHE STRING "Zoltan install directory")
    if(ZOLTAN_PREFIX)
        message(STATUS "ZOLTAN_PREFIX ${ZOLTAN_PREFIX}")
    endif()

    find_path(ZOLTAN_INCLUDE_DIR zoltan.h PATHS "${ZOLTAN_PREFIX}/include")
    find_library(ZOLTAN_LIBRARY zoltan PATHS "${ZOLTAN_PREFIX}/lib")  
endif()

set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY})
set(ZOLTAN_INCLUDE_DIRS ${ZOLTAN_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ZOLTAN_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    ZOLTAN
    DEFAULT_MSG
    ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR
)

mark_as_advanced(ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARY)