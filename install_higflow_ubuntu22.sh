echo 'Install:= Update ...'
sudo apt update 

echo 'Install:= Upgrade ...'
sudo apt -y upgrade

echo 'Install:= Glib ...' 
sudo apt -y install libglib2.0-dev 

echo 'Install:= boost ...'
sudo apt -y install libboost-all-dev

echo 'Install:= valgrind ...' 
sudo apt -y install valgrind

echo 'Install:= cmake ...' 
sudo apt -y install cmake

echo 'Install:= gfortran ...' 
sudo apt -y install gfortran

echo 'Install:= hypre ...' 
sudo apt -y install libhypre-dev

echo 'Install:= hdf5 ...' 
sudo apt -y install libhdf5-openmpi-dev

echo 'Install:= mpi ...' 
sudo apt -y install openmpi-bin sudo apt -y install mpi

echo 'Install:= libopenmpi ...' 
sudo apt -y install libopenmpi-dev

echo 'Install:= trilinus ...' 
sudo apt -y install libtrilinos-zoltan-dev

echo 'Install:= mpich ...'
sudo apt -y install mpich

echo 'Install:= libfyaml dependence ...' 
sudo apt -y install gcc autoconf automake libtool git make libltdl-dev pkg-config
sudo apt -y install libyaml-dev
sudo apt -y install check
sudo apt -y install python3 python3-pip python3-setuptools
pip3 install wheel sphinx git+http://github.com/return42/linuxdoc.git sphinx\_rtd\_theme sphinx-markdown-builder

echo 'Install:= libfyaml lib ...' 

cd bibliotecas/
unzip libfyaml-master.zip
sudo chmod 777 libfyaml-master
cd libfyaml-master/ ./bootstrap.sh ./configure make make check sudo make install 
cd ../..

echo 'Install:= petsc-3_14_0 ...' 
cd bibliotecas/
tar -vzxf petsc-3.14.0.tar.gz
cd petsc-3.14.0/ 
sudo ./configure --prefix=/opt/petsc-3.14.0-openmnpi-hypre-hdf5 --PETSC_ARCH=x86_64 --download-openmpi --download-hdf5 --download-hypre --download-fblaslapack --with-debubbing=yes --with-cc=gcc --with-cxx=g++ --with-fc=gfortran
sudo make PETSC_DIR=$PWD PETSC_ARCH=x86_64 all 
sudo make PETSC_DIR=$PWD PETSC_ARCH=x86_64 install 
sudo make PETSC_DIR=/opt/petsc-3.14.0-openmnpi-hypre-hdf5 PETSC_ARCH="" check
cd ../../

echo 'Configure hypre...'
sudo ln -s /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE.so /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE_krylov.so
export PKG_CONFIG_PATH=$PWD/bibliotecas/hypre/

echo 'Configure libfyaml ...'
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig/libfyaml.pc

echo 'Re-install Libfyaml...'
sudo apt -y install libyaml-dev
sudo apt autoremove

echo 'LDCONFIG for the libraries...'
sudo /sbin/ldconfig -v
