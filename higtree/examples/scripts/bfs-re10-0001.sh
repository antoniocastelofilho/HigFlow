#!/bin/bash -f

amrs=$1
steps=$2
res=$3
mul=$4
vtks=$5

./read-amr-solve-ns-example-2d \
	2 \
		$amrs/bfs-test-1-2d.amr \
		$amrs/bfs-test-2-2d.amr \
	6 \
		$amrs/bfs-test-bc-1-2d.amr \
		1 0 0 2 0 0   \
		$amrs/bfs-test-bc-2-2d.amr \
		1 0 0 0 0 0   \
		$amrs/bfs-test-bc-3-2d.amr \
		0 0 1 0 1 0   \
		$amrs/bfs-test-bc-4-2d.amr \
		1 0 0 0 0 0   \
		$amrs/bfs-test-bc-5-2d.amr \
		1 0 0 0 0 0   \
		$amrs/bfs-test-bc-6-2d.amr \
		1 0 0 0 0 0   \
	0.0001 10.0 $steps $res   \
	4 \
		$amrs/bfs-test-probe-3-2d.amr 0 \
		$amrs/bfs-test-probe-1-2d.amr 0 \
		$amrs/bfs-test-probe-2-2d.amr 0 \
		$amrs/bfs-test-probe-4-2d.amr 0 \
	$mul \
	$vtks \
	vtk -ksp_type bcgs -pc_type hypre

