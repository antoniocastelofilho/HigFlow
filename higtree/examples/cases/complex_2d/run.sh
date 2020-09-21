#!/bin/bash

cfl=0.5
simtime=60
freq=1
reynolds=200
output=out/re200
method=monolytic-newton

#xterm -e gdb --ex "set pagination off" --ex run --args \

mpirun -n 4 \
../../read-amr-solve-ns-cfl-example-2d \
	10 \
		ct-d-0.amr \
		ct-d-1.amr \
		ct-d-2.amr \
		ct-d-3.amr \
		ct-d-4.amr \
		ct-d-5.amr \
		ct-d-6.amr \
		ct-d-7.amr \
		ct-d-8.amr \
		ct-d-9.amr \
	38 \
		ct-bc-0.amr \
		1 0 0 1 0 0 \
		ct-bc-1.amr \
		1 0 0 0 0 0 \
		ct-bc-2.amr \
		1 0 0 0 0 0 \
		ct-bc-3.amr \
		1 0 0 0 0 0 \
		ct-bc-4.amr \
		1 0 0 0 0 0 \
		ct-bc-5.amr \
		1 0 0 0 0 0 \
		ct-bc-6.amr \
		1 0 0 0 0 0 \
		ct-bc-7.amr \
		1 0 0 0 0 -1 \
		ct-bc-8.amr \
		1 0 0 0 0 0 \
		ct-bc-9.amr \
		1 0 0 0 0 0 \
		ct-bc-10.amr \
		1 0 0 0 0 0 \
		ct-bc-11.amr \
		1 0 0 0 0 0 \
		ct-bc-12.amr \
		1 0 0 0 0 0 \
		ct-bc-13.amr \
		1 0 0 0 0 0 \
		ct-bc-14.amr \
		1 0 0 0 0 0 \
		ct-bc-15.amr \
		1 0 0 0 0 0 \
		ct-bc-16.amr \
		0 0 1 0 1 0 \
		ct-bc-17.amr \
		1 0 0 0 0 0 \
		ct-bc-18.amr \
		1 0 0 0 0 0 \
		ct-bc-19.amr \
		1 0 0 0 0 0 \
		ct-bc-20.amr \
		0 0 1 0 1 0 \
		ct-bc-21.amr \
		1 0 0 0 0 0 \
		ct-bc-22.amr \
		1 0 0 0 0 0 \
		ct-bc-23.amr \
		1 0 0 0 0 0 \
		ct-bc-24.amr \
		1 0 0 0 0 0 \
		ct-bc-25.amr \
		1 0 0 0 0 0 \
		ct-bc-26.amr \
		1 0 0 0 0 0 \
		ct-bc-27.amr \
		1 0 0 0 0 0 \
		ct-bc-28.amr \
		1 0 0 0 0 0 \
		ct-bc-29.amr \
		1 0 0 0 0 0 \
		ct-bc-30.amr \
		1 0 0 0 0 0 \
		ct-bc-31.amr \
		1 0 0 0 0 0 \
		ct-bc-32.amr \
		1 0 0 0 0 0 \
		ct-bc-33.amr \
		1 0 0 0 0 0 \
		ct-bc-34.amr \
		1 0 0 0 0 0 \
		ct-bc-35.amr \
		1 0 0 0 0 0 \
		ct-bc-36.amr \
		1 0 0 0 0 0 \
		ct-bc-37.amr \
		1 0 0 0 0 0 \
	1 \
		ct-probe-0.amr \
	$cfl $reynolds $simtime $freq $method $output \
	-ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition self -fieldsplit_0_ksp_type bcgs -fieldsplit_0_pc_type hypre -fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type lsc -fieldsplit_1_lsc_pc_type hypre -fieldsplit_1_lsc_pc_hypre_boomeramg_cycle_type w # -ksp_monitor_singular_value
