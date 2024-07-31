CMAKE_minimum_required(VERSION 3.12)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
project(example_Gptt)

add_executable(${PROJECT_NAME}     
    ${HIGFLOW_GPTT_DIR}/ns-example-2d.c
)

target_include_directories(${PROJECT_NAME} PRIVATE ${HIGFLOW_SRC_DIR} 
                                                        ${HIGTREE_SRC_DIR} 
                                                        ${HIGFLOW_GPTT_DIR})

target_link_libraries(${PROJECT_NAME}
                            higflow${dim}d
                            ${HIGFLOW_DEPENDENCIES})