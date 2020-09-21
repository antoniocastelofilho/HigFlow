#!/bin/bash -f

amrs=$1
steps=$2
res=$3
mul=$4
vtks=$5

./read-amr-solve-ns-example-2d \
	10 \
		$amrs/ct-d-0.amr \
		$amrs/ct-d-1.amr \
		$amrs/ct-d-2.amr \
		$amrs/ct-d-3.amr \
		$amrs/ct-d-4.amr \
		$amrs/ct-d-5.amr \
		$amrs/ct-d-6.amr \
		$amrs/ct-d-7.amr \
		$amrs/ct-d-8.amr \
		$amrs/ct-d-9.amr \
	38 \
		$amrs/ct-bc-0.amr \
		1 0 0 1 0 0 \
		$amrs/ct-bc-1.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-2.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-3.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-4.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-5.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-6.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-7.amr \
		1 0 0 0 0 -1 \
		$amrs/ct-bc-8.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-9.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-10.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-11.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-12.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-13.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-14.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-15.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-16.amr \
		0 0 1 0 1 0 \
		$amrs/ct-bc-17.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-18.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-19.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-20.amr \
		0 0 1 0 1 0 \
		$amrs/ct-bc-21.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-22.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-23.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-24.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-25.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-26.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-27.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-28.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-29.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-30.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-31.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-32.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-33.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-34.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-35.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-36.amr \
		1 0 0 0 0 0 \
		$amrs/ct-bc-37.amr \
		1 0 0 0 0 0 \
	0.001 1.0 $steps $res   \
	1 \
		$amrs/ct-probe-0.amr 0 \
	$mul \
	$vtks \
	vtk -ksp_type bcgs -pc_type hypre

