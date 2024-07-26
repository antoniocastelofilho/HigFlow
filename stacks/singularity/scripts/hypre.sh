#!/bin/bash

echo "*******************************************************"
echo ""
echo "              HYPRE "
echo ""
echo "*******************************************************"

if [ ! -f /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE.so ]
then 
    ln -s /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE_krylov.so /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE.so 
elif [ ! -f /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE_krylov.so ]
then 
    ln -s /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE.so /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE_krylov.so 
else
    echo ""
fi


