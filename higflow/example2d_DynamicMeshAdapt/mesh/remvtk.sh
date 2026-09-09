#!/bin/bash

# for the given folder, remove the vtk files in it
# receive folder from user
folder=$1
cd $folder
find -type f -name "*.vtk" -delete
echo "Removed all vtk files in $folder"
cd ..
