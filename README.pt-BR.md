# HigFlow
Após realiza o clone deste repositório será necessário instalar algumas dependências. 

Para realizar a instalação existe um script denominado 'install_higflow_ubuntu22' que contém todos os passos para instalar as dependências.
Basta dar permissão de execução ao arquivo fazendo:

* chmod +x install_higflow_ubuntu22

Excecute o arquivo:

* ./install_higflow_ubuntu22

Após a instalação terminar (pode levar bastante tempo), você pode testar se o código está funcionando fazendo:

* ./TUM_2D

(Será necessário dar permissão de execução para este arquivo também!)

## Pre-definições:

O arquivo 'varsrc' contém os caminhos de 'HIGTREE_DIR', 'HIGFLOW_DIR' e 'PETSC_DIR' (necessários para o funcionamento do código). Para cada terminal que abrir você deverá carregar esse arquivo fazendo:

* source varsrc

Você pode também editar o '.bashrc' para que carregue automaticamente o arquivo 'varsrc', adicionando uma linha com o comando no .bashrc com o seguinte comando:

* source /caminho/ate/o/arquivo/varsrc


Após re./configure --prefix=/opt/petsc-3.14.0-openmnpi-hypre-hdf5 --PETSC_ARCH=x86_64 --download-openmpi --download-hdf5 --download-hypre --download-fblaslapack --with-debubbing=yes --with-cc=gcc --with-cxx=g++ --with-fc=gfortran ;
sleep 5alizar um dos passos anterior você já pode utilizar o sistema!

## Usando o sistema HigFlow
Inicialmente é preciso ter um terminal aberto no diretório do sistema HigFlow. Navegue via terminal para o diretório 'higtree' e compile, fazendo:

* make clean && make DIM=2 && make DIM=3

Retorne ao diretório 'higflow' (no mesmo terminal), configure as variáveis de ambiente

* source ../etc/higflow-env-mak.sh

então, compile o código fazendo:

* make clean && make DIM=2

Posteriormente basta escolher qual exemplo (já adicionado) deseja estudar e executar (dentro do diretório):

* make clean && make && make run


## cmake with lmod (mflab'ss stack 13 - usp version)

Para compilar, cria a pasta build na raiz do projeto

Carrega os módulos

```
ml gnu/13.2.0 cmake glib boost libfyaml glib viennacl openmpi hdf5 petsc zoltan
```

Compila com a dimensão que for usar o higtree e higflow

```
cmake .. -Wno-dev -Ddim=2
make
make install
```

Opções do cmake

```
-Dprefix=PATH   Indica um path para instalar as libs e binários
-Ddebug=1       Ativa a compilação em modo debug
```

Para rodar vai na pasta onde instalou o código 

```
cd PASTA_ONDE_INSTALOU
```

carrega os módulos (se abrindo em outro terminal)

```
ml gnu/13.2.0 cmake glib boost libfyaml glib viennacl openmpi hdf5 petsc zoltan
```

Configura as variáveis de ambiente e rode

```
source ../etc/higflow-env.sh
mpirun ...
```
