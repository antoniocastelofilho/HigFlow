#!/bin/bash -f

amrs=$1

echo ./read-amr-solve-channel-ns-new-partition-example-2d 6 ${amrs}/test-1-2d.amr ${amrs}/test-2-2d.amr ${amrs}/test-3-2d.amr ${amrs}/test-4-2d.amr ${amrs}/test -5-2d.amr ${amrs}/test-6-2d.amr 3 ${amrs}/test-bc-1-2d.amr ${amrs}/test-bc-2-2d.amr ${amrs}/test-bc-3-2d.amr 1 0 1 0 0 0   1 0 1 0 0 0   1 0 0 0 1 0   0.00 1 10.0 10 1    1 ${amrs}/test-probe-1-2d.amr 0

./read-amr-solve-channel-ns-new-partition-example-2d \
	6 \
		${amrs}/test-1-2d.amr \
		${amrs}/test-2-2d.amr \
		${amrs}/test-3-2d.amr \
		${amrs}/test-4-2d.amr \
		${amrs}/test-5-2d.amr \
		${amrs}/test-6-2d.amr \
	3 \
		${amrs}/test-bc-1-2d.amr \
		${amrs}/test-bc-2-2d.amr \
		${amrs}/test-bc-3-2d.amr \
		1 1 1 0 0 0   \
		1 0 1 0 0 0   \
		1 0 0 0 1 0   \
	0.001 10.0 10 1    \
	1 \
		${amrs}/test-probe-1-2d.amr 0

echo ${amrs}
