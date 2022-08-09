# HigFlow
Após realiza o clone deste repositório será necessário instalar algumas dependências. 

Para realizar a instalação existe um script denominado 'install_higflow_ubuntu22' que contém todos os passos para instalar as dependências.
Basta dar permissão de execução ao arquivo fazendo:

* chmod +x install_higflow_ubuntu22

Excecute o arquivo em modo de root:

* sudo ./install_higflow_ubuntu22

Após a instalação terminar (pode levar bastante tempo), você pode testar se o código está funcionando fazendo:

* ./TUM_2D

(Será necessário dar permissão de execução para este arquivo também!)

## Pre-definições:

O arquivo 'varsrc' contém os caminhos de 'HIGTREE_DIR', 'HIGFLOW_DIR' e 'PETSC_DIR' (necessários para o funcionamento do código). Para cada terminal que abrir você deverá carregar esse arquivo fazendo:

* source varsrc

Você pode também editar o '.bashrc' para que carregue automaticamente o arquivo 'varsrc', adicionando uma linha com o comando no .bashrc com o seguinte comando:

* source /caminho/ate/o/arquivo/varsrc


Após realizar um dos passos anterior você já pode utilizar o sistema!

## Usando o sistema HigFlow

