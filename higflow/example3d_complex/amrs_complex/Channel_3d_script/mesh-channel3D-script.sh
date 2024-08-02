#!/bin/bash -f

amrs=$1
steps=$2
res=$3
mul=$4
vtks=$5

	./read-amr-solve-ns-example-3d \
		33 \
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
			amrs/Channel_3d_d/mesh-channel3D-d-16.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-17.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-18.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-19.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-20.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-21.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-22.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-23.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-24.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-25.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-26.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-27.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-28.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-29.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-30.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-31.amr \
			amrs/Channel_3d_d/mesh-channel3D-d-32.amr \
		169 \
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
			amrs/Channel_3d_bc/mesh-channel3D-bc-81.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-82.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-83.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-84.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-85.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-86.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-87.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-88.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-89.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-90.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-91.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-92.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-93.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-94.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-95.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-96.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-97.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-98.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-99.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-100.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-101.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-102.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-103.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-104.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-105.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-106.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-107.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-108.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-109.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-110.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-111.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-112.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-113.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-114.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-115.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-116.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-117.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-118.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-119.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-120.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-121.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-122.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-123.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-124.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-125.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-126.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-127.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-128.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-129.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-130.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-131.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-132.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-133.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-134.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-135.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-136.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-137.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-138.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-139.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-140.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-141.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-142.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-143.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-144.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-145.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-146.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-147.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-148.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-149.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-150.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-151.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-152.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-153.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-154.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-155.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-156.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-157.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-158.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-159.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-160.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-161.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-162.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-163.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-164.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-165.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-166.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-167.amr \
			1 0 0 0 0 0 0 0 \
			amrs/Channel_3d_bc/mesh-channel3D-bc-168.amr \
			1 0 0 0 0 0 0 0 \
		0.010000 10.000000 $steps $res \
		0 \
		$mul \
		$vtks \
		-ksp_type bcgs -pc_type hypre
