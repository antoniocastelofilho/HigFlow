<!--
Versao em portugues do README. O README.md em ingles e o documento principal.
Este texto vem da branch PC_Daniel_Mesh, de Pedro Coimbra, que ja havia
reescrito varias secoes e corrigido o trecho danificado.
-->

# HigFlow
Após realiza o clone deste repositório será necessário instalar algumas dependências. 

Para realizar a instalação existe um script denominado 'install_higflow_ubuntu22' que contém todos os passos para instalar as dependências.

Excecute o arquivo:

```bash
./install_higflow_ubuntu22
```

Após a instalação terminar (pode levar bastante tempo), você pode testar se o código está funcionando fazendo:

```bash
./TUM_2D
```

## Pre-definições:

O arquivo 'varsrc' contém os as variáveis de ambiente 'HIGTREE_DIR', 'HIGFLOW_DIR' e 'PETSC_DIR' que contém os caminhos das bibliotecas do HigFlow e da HiGTree  (necessários para o funcionamento do código). Para cada terminal que abrir você deverá carregar esse arquivo fazendo:

```bash
source varsrc
```

Você pode também editar o '.bashrc' para que carregue automaticamente o arquivo 'varsrc', adicionando uma linha com o comando no .bashrc com o seguinte comando:

```bash
source /caminho/ate/o/arquivo/varsrc
```
Após carregar as variáveis de ambiente, já é possível utlizar as bibliotecas.

## Usando o sistema HigFlow
Inicialmente é preciso ter um terminal aberto no diretório do sistema HigFlow. Navegue via terminal para o diretório 'higtree' e compile, fazendo:

```bash
make clean && make DIM=2 && make DIM=3
```

Retorne ao diretório 'higflow' (no mesmo terminal), configure as variáveis de ambiente

```bash
source ../etc/higflow-env-mak.sh
```

então, compile o código fazendo:

```bash
make clean && make DIM=2
```

Posteriormente basta escolher qual exemplo (já adicionado) deseja estudar e executar (dentro do diretório):

```bash
make clean && make && make run
```


## cmake with lmod (mflab'ss stack 13 - usp version)

Para compilar, cria a pasta build na raiz do projeto

Carrega os módulos

```bash
ml gnu/13.2.0 cmake glib boost libfyaml glib viennacl openmpi hdf5 petsc zoltan
```
Compila com a dimensão que for usar o higtree e higflow

```bash
cmake .. -Wno-dev -Ddim=2
make
make install
```

Opções do cmake

```bash
-Dprefix=PATH   Indica um path para instalar as libs e binários
-Ddebug=1       Ativa a compilação em modo debug
```

Para rodar vai na pasta onde instalou o código 
```
cd PASTA_ONDE_INSTALOU
```

carrega os módulos (se abrindo em outro terminal)

```bash
ml gnu/13.2.0 cmake glib boost libfyaml glib viennacl openmpi hdf5 petsc zoltan
```

Configura as variáveis de ambiente e rode

```bash
source ../etc/higflow-env.sh
mpirun ...
```

# Docker

Como configurar docker?

```bash
docker build . -t higflow:v01 -f conteiner/Dockerfile
```
