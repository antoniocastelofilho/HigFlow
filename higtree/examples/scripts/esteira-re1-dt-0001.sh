#!/bin/bash -f

amrs=$1
steps=$2
res=$3
mul=$4
vtks=$5

./read-amr-solve-ns-example-2d \
	4 \
		$amrs/esteira-d-0.amr \
		$amrs/esteira-d-1.amr \
		$amrs/esteira-d-2.amr \
		$amrs/esteira-d-3.amr \
	8 \
		$amrs/esteira-bc-0.amr \
		1 0 0 11 0 0   \
		$amrs/esteira-bc-1.amr \
		1 0 0 0 0 0   \
		$amrs/esteira-bc-2.amr \
		0 0 1 0 1 0   \
		$amrs/esteira-bc-3.amr \
		1 0 0 0 0 0   \
		$amrs/esteira-bc-4.amr \
		1 0 0 0 0 0   \
		$amrs/esteira-bc-5.amr \
		1 0 0 0 0 0   \
		$amrs/esteira-bc-6.amr \
		1 0 0 0 0 0   \
		$amrs/esteira-bc-7.amr \
		1 0 0 0 0 0   \
	0.0001 1.0 $steps $res   \
	3 \
		$amrs/esteira-probe-0.amr 0 \
		$amrs/esteira-probe-1.amr 0 \
		$amrs/esteira-probe-2.amr 0 \
	$mul \
	$vtks \
	vtk -ksp_type bcgs -pc_type hypre

