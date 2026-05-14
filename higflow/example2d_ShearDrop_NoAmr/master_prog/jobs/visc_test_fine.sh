#!/bin/bash

cd ~/HigFlow/
export HIGTREE_DIR=$(pwd)/higtree
export HIGFLOW_DIR=$(pwd)/higflow

### your script goes here

cd $HIGTREE_DIR
make clean
make DIM=2
cd $HIGFLOW_DIR
make clean
make DIM=2
cd example2d_ElectroOsmotic
make build_run IN="visc" MESH=fine ENAME=test NP=3
