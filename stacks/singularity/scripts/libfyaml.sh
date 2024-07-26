#!/bin/bash

echo "*******************************************************"
echo ""
echo "              LIBFYAML "
echo ""
echo "*******************************************************"

cd /bibliotecas
if [ -d libfyaml-master ]; then rm -rf libfyaml-master; fi
unzip libfyaml-master.zip
cd libfyaml-master
if ! ./bootstrap.sh 
then
    exit 1
fi
if ! ./configure -prefix=/opt/libfyaml
then
    exit 1
fi
if ! make
then
    exit 1
fi 
if ! make check 
then
    exit 1
fi
if ! make install 
then
    exit 1
fi

export LIBFYAML_INSTALL_DIR="/opt/libfyaml"
export LIBFYAML_DIR="${LIBFYAML_INSTALL_DIR}"
export LIBFYAML_LIB="${LIBFYAML_INSTALL_DIR}/lib"
export LIBFYAML_INC="${LIBFYAML_INSTALL_DIR}/include"

export LIBRARY_PATH="${LIBFYAML_LIB}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${LIBFYAML_LIB}:${LD_LIBRARY_PATH}"
    
export INCLUDE="${LIBFYAML_INC}:${INCLUDE}"
export CPLUS_INCLUDE_PATH="${LIBFYAML_INC}:${CPLUS_INCLUDE_PATH}"
export CPATH="${LIBFYAML_INC}:${CPATH}"
export C_INCLUDE_PATH="${LIBFYAML_INC}:${C_INCLUDE_PATH}"

export PKG_CONFIG_PATH="${LIBFYAML_LIB}/pkgconfig:${PKG_CONFIG_PATH}"