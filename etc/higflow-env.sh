#!/bin/bash
HIGFLOW_DIR="$(dirname "$PWD")"
HIGFLOW_BIN=$HIGFLOW_DIR/bin
HIGFLOW_INC=$HIGFLOW_DIR/include
HIGFLOW_LIB=$HIGFLOW_DIR/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HIGFLOW_LIB
export PATH=$PATH:$HIGFLOW_BIN
export INCLUDE=$INCLUDE:$HIGFLOW_INC
export HIGFLOW_DIR
export HIGTREE_DIR