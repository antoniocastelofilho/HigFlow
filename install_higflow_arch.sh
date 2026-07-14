echo 'Pacotes nos repos core e extra'
sudo pacman -Syu --noconfirm --needed glib2 boost valgrind cmake gcc-fortran hdf5-openmpi openmpi autoconf automake libtool git make pkg-config libyaml

echo 'Pacotes na AUR'
yay -S --noconfirm  --needed hypre zoltan mpich 

python -m venv venv
source venv/bin/activate
pip3 install wheel sphinx sphinx_rtd_theme sphinx-markdown-builder xdrlib3
deactivate

echo 'libfyaml lib'
cd bibliotecas/
unzip libfyaml-master.zip
chmod +x libfyaml-master
cd libfyaml-master/
./bootstrap.sh
./configure
make
make check
sudo make install

cd ../..
echo 'Install petsc-3.23.6' 
cd bibliotecas/
wget https://aur.archlinux.org/cgit/aur.git/snapshot/petsc.tar.gz
tar -vzxf petsc.tar.gz
cd petsc/
makepkg -si

cd ../..
echo "Compila HigTree"
cd higtree/

echo "Compila HigFlow"
make DIM=2 && make DIM=3
cd ../higflow/
make
