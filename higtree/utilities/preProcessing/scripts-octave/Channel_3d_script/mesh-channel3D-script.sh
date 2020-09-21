#!/bin/bash -f

amrs=$1
steps=$2
res=$3
mul=$4
vtks=$5

	./read-amr-solve-ns-example-3d \
		16 \
			amrs/Channel_3d_d/mesh-channel3D-d-0.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-1.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-2.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-3.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-4.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-5.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-6.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-7.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-8.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-9.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-10.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-11.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-12.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-13.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-14.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-15.amr \
		81 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-0.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-1.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-2.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-3.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-4.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-5.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-6.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-7.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-8.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-9.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-10.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-11.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-12.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-13.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-14.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-15.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-16.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-17.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-18.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-19.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-20.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-21.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-22.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-23.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-24.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-25.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-26.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-27.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-28.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-29.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-30.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-31.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-32.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-33.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-34.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-35.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-36.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-37.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-38.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-39.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-40.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-41.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-42.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-43.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-44.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-45.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-46.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-47.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-48.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-49.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-50.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-51.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-52.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-53.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-54.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-55.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-56.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-57.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-58.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-59.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-60.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-61.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-62.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-63.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-64.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-65.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-66.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-67.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-68.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-69.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-70.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-71.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-72.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-73.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-74.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-75.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-76.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-77.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-78.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-79.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-80.amr \
			1 0 0 0 0 0 0 0 \
		0.010000 10.000000 $steps $res \
		0 \
		$mul \
		$vtks \
		-ksp_type bcgs -pc_type hypre
