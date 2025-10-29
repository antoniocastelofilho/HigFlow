# Dockerfile to build Ubuntu22.04x64 + OpenFOAM-9 + Python 3.10

FROM ubuntu_petsc3.14:v01

RUN                                                                            \
cd bibliotecas/                                                              &&\
unzip libfyaml-master.zip                                                    &&\
chmod 777 libfyaml-master                                                    &&\
cd libfyaml-master/                                                          &&\
./bootstrap.sh && ./configure && make && sudo make install                   &&\
sudo apt-get update && sudo apt -y install libyaml-dev                       &&\
cd /home/hig_user/HigFlow                                                    &&\
source varsrc                                                                &&\
(cd higtree/ && make clean && make DIM=2 && make DIM=3)                      &&\
sudo /sbin/ldconfig -v                                                       &&\
echo "source /home/hig_user/HigFlow/varsrc" >> /home/hig_user/.bashrc
# (cd higflow/; make clean && make DIM=2);                                       \
# cd ../../                                                                    &&\
# echo 'Configure hypre...'                                                    &&\
# sudo ln -s /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE.so             \
# /opt/petsc-3.14.0-openmnpi-hypre-hdf5/lib/libHYPRE_krylov.so                 &&\
# export PKG_CONFIG_PATH=$PWD/bibliotecas/hypre/                               &&\
# sudo /sbin/ldconfig -v
