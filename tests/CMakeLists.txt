# Verification tests, registered with CTest.
#
# Standalone, so it needs nothing from the project's own CMake configuration:
#
#     cmake -S tests -B build-tests
#     ctest --test-dir build-tests --output-on-failure
#
# It is written this way on purpose. The project has two build systems that
# compile different source sets, and wiring tests into either of them is a
# change to how the project builds, which is a separate conversation. Once
# there is one canonical configuration, `add_subdirectory(tests)` from the root
# picks these up unchanged.

cmake_minimum_required(VERSION 3.12)
project(higflow_verification NONE)

enable_testing()

find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(VERIFY "${CMAKE_CURRENT_SOURCE_DIR}/verification/verify.py")

# numpy is the only third-party dependency, and a missing one should say so
# rather than surface as an import traceback in the middle of a test run.
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy"
    RESULT_VARIABLE HAVE_NUMPY
    OUTPUT_QUIET ERROR_QUIET)
if(NOT HAVE_NUMPY EQUAL 0)
    message(WARNING "numpy not found for ${Python3_EXECUTABLE}; the verification tests will not run")
endif()

# HIGFLOW_VTKS points at the output of a run that has already happened. Without
# it the whole-case test runs the case itself, which needs the container.
set(HIGFLOW_VTKS "" CACHE PATH "VTKS directory of an existing run of example2d_Newt")

if(HIGFLOW_VTKS)
    add_test(NAME channel_steady_state
             COMMAND ${Python3_EXECUTABLE} ${VERIFY} transient --vtks ${HIGFLOW_VTKS} --every 10)
    set_tests_properties(channel_steady_state PROPERTIES
        PASS_REGULAR_EXPRESSION "the run is at steady state"
        LABELS "verification")

    add_test(NAME channel_exact_solution
             COMMAND ${Python3_EXECUTABLE} ${VERIFY} check --vtks ${HIGFLOW_VTKS})
    set_tests_properties(channel_exact_solution PROPERTIES LABELS "verification")
else()
    # Runs the case and then checks it. Slower, and needs docker.
    add_test(NAME channel_exact_solution
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/run_verification.sh)
    set_tests_properties(channel_exact_solution PROPERTIES
        LABELS "verification;slow"
        TIMEOUT 1800)
endif()

# The convergence study needs four runs at different resolutions, so it is not
# registered by default. Give it the directories and it becomes a test:
#
#     cmake -S tests -B build-tests -DHIGFLOW_ORDER_RUNS="40=/a;80=/b;160=/c;320=/d"
set(HIGFLOW_ORDER_RUNS "" CACHE STRING "NX=DIR entries for the convergence study")

if(HIGFLOW_ORDER_RUNS)
    set(ORDER_ARGS "")
    foreach(entry ${HIGFLOW_ORDER_RUNS})
        list(APPEND ORDER_ARGS --run ${entry})
    endforeach()
    add_test(NAME channel_convergence_order
             COMMAND ${Python3_EXECUTABLE} ${VERIFY} order ${ORDER_ARGS}
                     --skip-edges 1 --expect 2.0 --order-tol 0.4)
    set_tests_properties(channel_convergence_order PROPERTIES
        LABELS "verification;slow")
endif()
