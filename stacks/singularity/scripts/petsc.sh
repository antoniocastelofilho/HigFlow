#!/bin/bash

echo "*******************************************************"
echo ""
echo "              PETSC "
echo ""
echo "*******************************************************"

cd /bibliotecas
if [ -d petsc-3.14.0 ]; then rm -rf petsc-3.14.0; fi
tar -vzxf petsc-3.14.0.tar.gz
cd petsc-3.14.0
if ! ./configure --prefix=/opt/petsc-3.14.0-openmnpi-hypre-hdf5 \
                    --PETSC_ARCH=x86_64 \
                    --download-openmpi \
                    --download-hdf5 \
                    --download-hypre \
                    --download-fblaslapack \
                    --with-debubbing=yes \
                    --with-cc=gcc \
                    --with-cxx=g++ \
                    --with-fc=gfortran
then
    cp configure.log /pacotes/.
    exit 1
fi
if ! make PETSC_DIR=`pwd` PETSC_ARCH=x86_64 all 
then
    exit 1
fi
if ! make PETSC_DIR=`pwd` PETSC_ARCH=x86_64 install
then
    exit 1
fi
if ! make PETSC_DIR=/opt/petsc-3.14.0-openmnpi-hypre-hdf5 PETSC_ARCH="" check
then
    exit 1
fi

export PETSC_INSTALL_DIR="/opt/petsc-3.14.0-openmnpi-hypre-hdf5"
export PETSC_DIR="${PETSC_INSTALL_DIR}"
export PETSC_LIB="${PETSC_INSTALL_DIR}/lib"
export PETSC_INC="${PETSC_INSTALL_DIR}/include"

export LIBRARY_PATH="${PETSC_LIB}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PETSC_LIB}:${LD_LIBRARY_PATH}"
    
export INCLUDE="${PETSC_INC}:${INCLUDE}"
export CPLUS_INCLUDE_PATH="${PETSC_INC}:${CPLUS_INCLUDE_PATH}"
export CPATH="${PETSC_INC}:${CPATH}"
export C_INCLUDE_PATH="${PETSC_INC}:${C_INCLUDE_PATH}"

export PKG_CONFIG_PATH="${PETSC_LIB}/pkgconfig:${PKG_CONFIG_PATH}"
