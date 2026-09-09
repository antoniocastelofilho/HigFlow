#!/bin/bash

# for all folders in the current folder, do "./read-amr-write-vtk file.amr file.vtk" for each .amr file that has "-d" in their name in them
for folder in $(ls -d */); do
    ./2vtk.sh $folder
done
