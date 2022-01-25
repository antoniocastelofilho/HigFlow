#!/bin/bash -f

amrs=amrs/bfs-test-refined
dt=$1
Re=$2
steps=$3
res=$4
mul=$5

./read-amr-solve-channel-ns-new-partition-example-2d \
	2 \
		$amrs/bfs-test-1-2d.amr \
		$amrs/bfs-test-2-2d.amr \
	6 \
		$amrs/bfs-test-bc-1-2d.amr \
		1 0 0 1 0 0   \
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
	$dt $Re $steps $res   \
	4 \
		$amrs/bfs-test-probe-3-2d.amr 0 \
		$amrs/bfs-test-probe-1-2d.amr 0 \
		$amrs/bfs-test-probe-2-2d.amr 0 \
		$amrs/bfs-test-probe-4-2d.amr 0 \
	$mul \
	-ksp_type bcgs -pc_type hypre

