#!/bin/bash -f

steps=$1
res=$2

./read-amr-solve-channel-ns-new-partition-example-2d \
	2 \
		amrs/test-bfs-1-2d.amr \
		amrs/test-bfs-2-2d.amr \
	7 \
		amrs/test-bfs-bc-1-2d.amr \
		1 0 0 10 0 0   \
		amrs/test-bfs-bc-2-2d.amr \
		1 0 0 0 0 0   \
		amrs/test-bfs-bc-3-2d.amr \
		1 0 0 0 0 0   \
		amrs/test-bfs-bc-4-2d.amr \
		0 0 0 0 0 0   \
		amrs/test-bfs-bc-5-2d.amr \
		1 0 0 0 0 0   \
		amrs/test-bfs-bc-6-2d.amr \
		1 0 0 0 0 0   \
		amrs/test-bfs-bc-7-2d.amr \
		1 0 0 0 0 0   \
	0.00001 10.0 $steps $res   \
	3 \
		amrs/test-bfs-probe-1-2d.amr 0 \
		amrs/test-bfs-probe-2-2d.amr 0 \
		amrs/test-bfs-probe-3-2d.amr 0 \
	5.0
