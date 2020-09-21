#/bin/sh

ulimit -s unlimited

#xterm -e gdb --ex "set pagination off" --ex run --args \
#valgrind --main-stacksize=304934672 \

mpirun -n 1 \
../../read-amr-solve-ns-cfl-example-3d \
13 \
d_000.amr \
d_001.amr \
d_002.amr \
d_003.amr \
d_004.amr \
d_005.amr \
d_006.amr \
d_007.amr \
d_008.amr \
d_009.amr \
d_010.amr \
d_011.amr \
d_012.amr \
54 \
bc_000.amr 1 0 0 0 0 0 0 0 \
bc_001.amr 1 0 0 0 0 0 0 0 \
bc_002.amr 1 0 0 0 0 0 0 0 \
bc_003.amr 1 0 0 0 0 0 0 0 \
bc_004.amr 1 0 0 0 0 0 0 0 \
bc_005.amr 1 0 0 0 0 0 0 0 \
bc_006.amr 1 0 0 0 0 0 0 0 \
bc_007.amr 1 0 0 0 0 0 0 0 \
bc_008.amr 1 0 0 0 0 0 0 0 \
bc_009.amr 1 0 0 0 0 0 0 0 \
bc_010.amr 1 0 0 0 0 0 0 0 \
bc_011.amr 1 0 0 0 0 0 0 0 \
bc_012.amr 1 0 0 0 0 0 0 0 \
bc_013.amr 1 0 0 0 0 0 0 0 \
bc_014.amr 1 0 0 0 0 0 0 0 \
bc_015.amr 1 0 0 0 0 0 0 0 \
bc_016.amr 1 0 0 0 0 0 0 0 \
bc_017.amr 1 0 0 0 0 0 0 0 \
bc_018.amr 1 0 0 0 0 0 0 0 \
bc_019.amr 1 0 0 0 0 0 0 0 \
bc_020.amr 1 0 0 0 0 0 0 0 \
bc_021.amr 1 0 0 0 0 0 0 0 \
bc_022.amr 0 0 0 1 0 0 0 0 \
bc_023.amr 1 0 0 0 0 0 0 0 \
bc_024.amr 1 0 0 0 0 0 0 0 \
bc_025.amr 1 0 0 0 0 0 0 0 \
bc_026.amr 1 0 0 0 0 0 0 0 \
bc_027.amr 1 0 0 0 0 0 0 0 \
bc_028.amr 1 0 0 0 0 0 0 0 \
bc_029.amr 1 0 0 0 0 0 0 0 \
bc_030.amr 1 0 0 0 0 0 0 0 \
bc_031.amr 1 0 0 0 0 0 0 0 \
bc_032.amr 1 0 0 0 0 0 0 0 \
bc_033.amr 1 0 0 0 0 0 0 0 \
bc_034.amr 1 0 0 0 0 0 0 0 \
bc_035.amr 1 0 0 0 0 0 0 0 \
bc_036.amr 1 0 0 0 0 0 0 0 \
bc_037.amr 1 0 0 0 0 0 0 0 \
bc_038.amr 1 0 0 0 0 0 0 0 \
bc_039.amr 1 0 0 0 0 0 0 0 \
bc_040.amr 1 0 0 0 0 0 0 0 \
bc_041.amr 1 0 0 0 0 0 0 0 \
bc_042.amr 1 0 0 0 0 0 0 0 \
bc_043.amr 1 0 0 0 0 0 0 0 \
bc_044.amr 1 0 0 0 0 0 0 0 \
bc_045.amr 1 0 0 0 0 0 0 0 \
bc_046.amr 1 0 0 0 0 0 0 0 \
bc_047.amr 1 0 0 0 0 0 0 0 \
bc_048.amr 1 0 0 0 0 0 0 0 \
bc_049.amr 1 0 0 0 0 0 0 0 \
bc_050.amr 1 0 0 0 0 0 0 0 \
bc_051.amr 1 0 0 1 0 0 0 0 \
bc_052.amr 1 0 0 0 0 0 0 0 \
bc_053.amr 1 0 0 0 0 0 0 0 \
0 \
0.20 200 0.04 125 out/out
