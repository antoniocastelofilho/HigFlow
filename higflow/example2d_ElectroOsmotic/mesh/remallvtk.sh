#!/bin/bash

# if remvtk.sh is not here, error
if [ ! -f remvtk.sh ]; then
    echo "remvtk.sh is not here"
    echo "Please run this script in the /mesh folder unless you want to delete all your precious .vtk files"
    exit 1
fi

# for for all folders, remove the vtk files in them
for folder in $(ls -d */); do
    ./remvtk.sh $folder
done
