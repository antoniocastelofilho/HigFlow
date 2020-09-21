#!/bin/bash -f

amrs=$1
steps=$2
res=$3
mul=$4
vtks=$5
loadvel=$6
savevel=$7

./read-amr-solve-ns-restart-example-2d \
	1 \
		$loadvel \
	1 \
		$savevel \
	1 \
		$amrs/ch-d-0.amr \
	4 \
		$amrs/ch-bc-0.amr \
		1 0 0 1 0 0 \
		$amrs/ch-bc-1.amr \
		1 0 0 0 0 0 \
		$amrs/ch-bc-2.amr \
		0 0 1 0 1 0 \
		$amrs/ch-bc-3.amr \
		1 0 0 0 0 0 \
	0.001 100.0 $steps $res   \
	2 \
		$amrs/ch-probe-0.amr 0 \
		$amrs/ch-probe-1.amr 0 \
	$mul \
	$vtks \
	-ksp_type bcgs -pc_type hypre

