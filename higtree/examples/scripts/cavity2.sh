#!/bin/bash -f

amrs=arms
dt=$1
Re=$2
steps=$3
res=$4

./read-amr-solve-channel-ns-new-partition-example-2d \
	1 \
		amrs/test-ch-1-2d.amr \
	4 \
		amrs/test-ch-bc-1-2d.amr \
		1 0 0 1 0 0   \
		amrs/test-ch-bc-2-2d.amr \
		1 0 0 0 0 0   \
		amrs/test-ch-bc-3-2d.amr \
		0 0 1 0 0 0   \
		amrs/test-ch-bc-4-2d.amr \
		1 0 0 0 0 0   \
	$dt $Re $steps $res   \
	2 \
		amrs/test-ch-probe-2-2d.amr 0 \
		amrs/test-ch-probe-1-2d.amr 0 \
	20.0
