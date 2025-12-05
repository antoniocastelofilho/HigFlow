#!/bin/bash

# for the given  folder, do "./read-amr-write-vtk file.amr file.vtk" for each .amr file that has "-d" in their name in them
# receive folder from user
folder=$1
cd $folder
cd domain
for file in $(ls *-d*.amr); do
    echo "Processing $file"
    ./../../read-amr-write-vtk $file ${file%.amr}.vtk
done
cd ../..
