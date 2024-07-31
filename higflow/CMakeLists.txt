CMAKE_minimum_required(VERSION 3.12)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
project(higflow${dim}d)

add_library(${PROJECT_NAME} STATIC
    ${HIGFLOW_SRC_DIR}/hig-flow-kernel.c
	${HIGFLOW_SRC_DIR}/hig-flow-io.c
	${HIGFLOW_SRC_DIR}/hig-flow-ic.c
	${HIGFLOW_SRC_DIR}/hig-flow-bc.c
	${HIGFLOW_SRC_DIR}/hig-flow-eval.c
	${HIGFLOW_SRC_DIR}/hig-flow-discret.c
	${HIGFLOW_SRC_DIR}/hig-flow-terms.c
	${HIGFLOW_SRC_DIR}/hig-flow-res.c
	${HIGFLOW_SRC_DIR}/hig-flow-step.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-generalized-newtonian.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-multiphase.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-multiphase-viscoelastic.c
	${HIGFLOW_SRC_DIR}/hig-flow-vof-plic.c
	${HIGFLOW_SRC_DIR}/hig-flow-vof-adap-hf.c
	${HIGFLOW_SRC_DIR}/hig-flow-vof-finite-difference-normal-curvature.c
    ${HIGFLOW_SRC_DIR}/hig-flow-vof-elvira.c
    ${HIGFLOW_SRC_DIR}/hig-flow-vof-9-cells.c
    ${HIGFLOW_SRC_DIR}/hig-flow-mittag-leffler.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-viscoelastic.c
    ${HIGFLOW_SRC_DIR}/hig-flow-step-viscoelastic-integral.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-electroosmotic.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-shear-thickening-suspension.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-viscoelastic-shear-banding.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-viscoelastic-variable-viscosity.c
	${HIGFLOW_SRC_DIR}/hig-flow-step-elastoviscoplastic.c
)

target_include_directories(${PROJECT_NAME} PRIVATE ${HIGFLOW_SRC_DIR} ${HIGTREE_SRC_DIR})

target_link_libraries(${PROJECT_NAME} ${HIGFLOW_DEPENDENCIES} hig${dim}d)