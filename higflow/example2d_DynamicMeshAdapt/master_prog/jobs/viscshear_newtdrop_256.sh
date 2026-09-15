#!/bin/bash

cd ~/HigFlow/
export PYTHONPATH=$(pwd)/bibliotecas
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
make build_run IN="mult visc-newt dt=0.0005 Re=4.8 Ca=2.4 De=0.4 beta=0.5 tf=10.0" MESH=short_256 NP=40 ENAME=pdir_IF_ss_73
