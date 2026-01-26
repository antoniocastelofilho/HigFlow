# HigTree

## Formato `amr`

* Adaptative Mesh Refinement

Formato para input da malha e da fronteira.

Exemplo `example2d_Newt/output/ch-bc-0.amr`:

Nesse caso temos uma higtree degenerada, ou seja um domínio em apenas uma
dimensão para o contorno do domínio.

```
0.0 0.0 -1.0 1.0
1
0.0 0.05 1
1 1 1 40
```

1. Primeira linha:
  * Boundinbox global
  ```
  x1  x2   y1  y2
  0.0 0.0 -1.0 1.0
  ```

2. Segunda linha:
  * Número de higtree
   ```
   1
   ```
3. Terceira linha
  * Início das especificações de uma das higtrees
  1) Dimensões das células
  ```
  x   y    n
  0.0 0.05 1
  ```
  onde `n` é o número de refinamentos da higtree.

4.  Quarta linha
  * Quais pedaços da higtree vamos preencher com malha.
  ```
  1 1 1 40
  ```

# Arquivos

# `higtree.c`

# Documentação do Código da Biblioteca HiG-Tree

### 1. Visão Geral

A biblioteca fornece um conjunto de ferramentas para realizar simulações
numéricas em domínios complexos, que são discretizados usando uma estrutura de
dados de malha cartesiana hierárquica chamada **HiG-Tree**. Ela foi projetada
com o paralelismo em mente, utilizando o padrão MPI para distribuição de dados
e computação.

**Principais Características:**
* **Malha Hierárquica (HiG-Tree):** Permite o refinamento adaptativo da malha,
concentrando a resolução computacional apenas onde é necessário, economizando
memória e tempo de processamento.
* **Decomposição de Domínio e Balanceamento de Carga:** Utiliza a biblioteca
Zoltan para particionar a malha de forma eficiente entre múltiplos processos,
garantindo que a carga de trabalho seja distribuída de maneira equilibrada.
* **Cálculo de Stencil de Alto Nível:** Abstrai o complexo processo de cálculo
de stencils para métodos de diferenças finitas.
Utiliza o método de **Mínimos Quadrados Ponderados (WLS)** para
gerar stencils de alta ordem em pontos arbitrários, mesmo em malhas não
uniformes e perto de fronteiras.
* **Interface Abstrata de Solver Linear:** Fornece uma interface unificada para
interagir com vários solucionadores de sistemas lineares de alto desempenho,
como PETSc e HYPRE, permitindo flexibilidade na escolha do algoritmo de
solução.
* **Entrada/Saída Paralela:** Inclui funcionalidades robustas para salvar e
carregar dados de simulação em paralelo usando o formato HDF5, além de gerar
arquivos de visualização em formatos como VTK e XDMF para análise de dados com
ferramentas como ParaView e VisIt.
* **Gerenciamento de Domínio e Condições de Contorno:** Oferece estruturas
claras para definir o domínio computacional (composto por múltiplas HiG-Trees),
incluindo as "franjas" (células fantasmas) para comunicação entre processos e a
aplicação de condições de contorno como Dirichlet e Neumann.

---
## Módulos e Arquivos Associados

### Módulos Principais
---
#### **Malha Hierárquica: `higtree.h`, `higtree.c`**
A HiG-Tree é a estrutura de dados fundamental que representa a malha. Ela é uma
árvore onde cada nó (`hig_cell`) representa uma região retangular do espaço
(célula).

* **`hig_cell`**: Uma célula que pode ser uma folha (a menor resolução em uma
região) ou um nó interno que é subdividido em células filhas.
Cada célula armazena seus limites (`lowpoint`, `highpoint`), seu
pai, e ponteiros para seus filhos.
* **`hig_facet`**: Representa uma face (ou aresta em 2D, lado em 1D) de uma
`hig_cell`.
* **Refinamento**: Uma célula pode ser refinada usando `hig_refine_uniform`,
que a subdivide em um grid regular de células filhas.
* **Consultas**: A biblioteca oferece funções para consultar propriedades da
malha, como encontrar a célula que contém um determinado ponto
(`hig_get_cell_with_point`) , obter o centro de uma célula (`hig_get_center`) ,
ou encontrar células vizinhas (`hig_get_neighbour`).
---
#### **Iteradores da Malha: `higtree-iterator.h`, `higtree-iterator.c`, `higtree-iterator-internal.h`**
Para navegar pela complexa estrutura da HiG-Tree, a biblioteca oferece um
sistema de iteradores flexível e poderoso.

* **`higcit_celliterator`**: Um iterador para percorrer células (`hig_cell`).
* **`higfit_facetiterator`**: Um iterador para percorrer faces (`hig_facet`).
* **Tipos de Iteradores**: É possível criar iteradores para diferentes propósitos:
    * `higcit_create_all_leaves`: Itera sobre todas as células-folha da árvore.
    * `higcit_create_bounding_box`: Itera sobre as folhas dentro de uma caixa de delimitação.
    * `higcit_create_neighbours`: Itera sobre as células vizinhas a uma dada célula.
    * `higcit_create_concat`: Combina múltiplos iteradores em um só.

---
#### **Domínio de Simulação: `domain.h`, `domain.c`**
Esta estrutura organiza as HiG-Trees e as condições de contorno associadas.

* **`sim_domain`**: Um contêiner que agrupa um conjunto de HiG-Trees que
compõem o domínio computacional. Ele distingue entre árvores que fazem parte do
domínio local e árvores que representam a "franja" (células fantasmas de
processos vizinhos).
* **`sim_boundary`**: Representa uma condição de contorno (CC). Cada CC é
associada a uma HiG-Tree que define sua geometria e armazena os valores da
condição (ex: valores de Dirichlet ou gradientes de Neumann). Funções como
`sb_create` e `sd_add_boundary` permitem a criação e associação de CCs ao
domínio.
* **`sim_facet_domain`**: Uma especialização do `sim_domain` projetada para
trabalhar com dados centrados nas faces das células, como fluxos ou
velocidades.

---
#### **Stencils: `domain.h`, `global-stencil.h`, `global-stencil.c`**
O stencil é uma representação da dependência de um ponto em relação aos seus
vizinhos, fundamental para a discretização de equações diferenciais.

* **`sim_stencil`**: Uma estrutura que armazena os coeficientes do stencil como
uma lista de pares `(ID_vizinho, peso)` e um termo independente (`rhs`). É
definida em `domain.h`.
* **`global_stencil`**: Uma estrutura similar ao `sim_stencil`, mas que utiliza
IDs globais para representar o stencil no domínio particionado.
* **Cálculo de Stencil**: As funções `sd_get_stencil` (para dados em células) e
`sfd_get_stencil` (para dados em faces) são o coração do motor numérico. Elas
calculam os pesos de interpolação para um ponto arbitrário
`x`, encontrando as células/faces vizinhas mais relevantes e usando WLS para
determinar os pesos.
---
#### **Mínimos Quadrados Ponderados (WLS): `wls.h`, `wls.c`**
O WLS é o método matemático usado para garantir a precisão da interpolação em
malhas não uniformes.

* **`wls_interpolator`**: Implementa o algoritmo WLS. Ele constrói um pequeno
sistema linear localmente, baseado em uma base polinomial e nos pontos
vizinhos, e o resolve para encontrar os pesos do stencil. A função
`wls_set_samples_and_calc` é tipicamente chamada internamente pelas rotinas de
`get_stencil`.

---
### Módulos de Paralelismo
---
#### **Balanceamento de Carga: `lbal.h`, `lbal.c`**
O `load_balancer` é responsável por particionar o domínio global entre os
processos MPI.

* **`load_balancer`**: Orquestra o particionamento da malha.
* **Processo de Particionamento**:
    1.  Todos os processos criam um contexto `load_balancer` com `lb_create`.
    2.  As árvores a serem particionadas são adicionadas com
        `lb_add_input_tree`.
    3.  A função coletiva `lb_calc_partition` é chamada para calcular a
        distribuição.
    4.  Cada processo obtém suas árvores de domínio local e as árvores de
        franja com `lb_get_local_tree`.

---
#### **Domínio Particionado: `pdomain.h`, `pdomain.c`**
Após o balanceamento de carga, o domínio de cada processo é representado por
estruturas `psim_*`.

* **`partition_graph`**: Armazena o resultado do particionamento, descrevendo
quais processos são vizinhos e as regiões de fronteira entre eles.
* **`psim_domain`**: Representa a visão de um processo sobre o domínio
distribuído. Ele contém o `sim_domain` local (com domínio e
franjas) e o `partition_graph` para entender a topologia da vizinhança.
* **`distributed_property`**: Uma abstração para um vetor de dados distribuído
entre os processos (ex: pressão, temperatura). A função
`dp_sync` é usada para atualizar os valores nas células de franja,
comunicando-se com os processos vizinhos.

---
#### **Comunicação e Sincronização: `higtree-parallel.h`, `higtree-parallel.c`,
`mapper-syncer.h`, `mapper-syncer.c`, `isend-pool.h`, `isend-pool.c`,
`term-det.h`, `build-fringe.h`**

Estes módulos fornecem as primitivas para a comunicação MPI.

* **`higtree-parallel`**: Funções para enviar e receber HiG-Trees inteiras ou
apenas suas estruturas de refinamento entre processos (e.g.,
`higp_send_uniform_tree_to_node`).
* **`mapper-syncer`**: Sincroniza os `mp_mapper` entre processos vizinhos para
garantir que os IDs de células/faces nas franjas sejam consistentes.
* **`isend-pool`**: Uma abstração para gerenciar múltiplas chamadas `MPI_Isend`
não bloqueantes, permitindo esperar por todas elas de uma só vez com
`isend_pool_wait_all`.
* **`term-det.h`**: Implementa um algoritmo de detecção de término para
coordenar atividades assíncronas entre processos.
* **`build-fringe.h`**: Contém a lógica para construir as regiões de franja
(células fantasmas) após o particionamento do domínio.

---
### Solvers Lineares
---
#### **Interface Abstrata: `solver.h`, `solver.c`**
Abstrai a interação com diferentes bibliotecas de resolução de sistemas lineares. Os arquivos `solver.camila.c` e `solver-petsc.camila.c` são variantes específicas de desenvolvimento.

* [cite_start]**Interface `solver`**: Define uma API genérica para configurar e resolver um sistema `Ax = b`[cite: 2982].
* **Funções Principais**:
    * [cite_start]`slv_create`: Cria uma instância de um solver, escolhendo um backend[cite: 1730, 1788, 2987].
    * [cite_start]`slv_set_Ai` / `slv_set_Aij`: Define os coeficientes da matriz `A`[cite: 1697, 1756, 1700, 1759].
    * [cite_start]`slv_assemble`: Finaliza a montagem da matriz e do vetor[cite: 1725, 1784].
    * [cite_start]`slv_solve`: Resolve o sistema[cite: 1726, 1785].

#### **Backends de Solvers**
* **`solver-petsc.h`, `solver-petsc.c`**: Implementação usando a biblioteca PETSc. Suporta matrizes esparsas e uma vasta gama de solvers iterativos e precondicionadores.
* **`solver-hypre.c`**: Implementação usando a biblioteca HYPRE, focada em precondicionadores de multigrid algébrico (AMG).
* **`solver-centralized.h`, `solver-centralized.c`**: Um framework para solvers que não são nativamente paralelos. Ele centraliza a matriz em um único processo (rank 0), resolve o sistema e distribui a solução de volta.
* **`solver-sor.c`**: Uma implementação do solver SOR (Successive Over-Relaxation) que utiliza o framework centralizado.
* **`solver-debug-write.c`**: Um "solver" para depuração que, em vez de resolver o sistema, escreve a matriz e o vetor RHS em arquivos para análise externa.

---
### Módulos de Entrada/Saída e Serialização
---
#### **I/O de Dados e Visualização: `higtree-io.h`, `higtree-io.c`**
* **Formatos Simples**: Funções como `higio_write_to_file` e `higio_read_from_file` salvam e carregam a estrutura da malha em um formato de texto simples.
* **Visualização com VTK**: `higio_print_in_vtk` exporta a malha no formato legado VTK, que pode ser lido por muitas ferramentas de visualização.
* **I/O Paralelo com HDF5**:
    * `higio_create_hdf5`: Cria um arquivo HDF5 para escrita.
    * `higio_write_tree_hdf5`: Salva a estrutura da malha em HDF5.
    * `higio_write_cell_property_hdf5` e `higio_write_facet_property_hdf5`: Salvam os dados de `distributed_property` em paralelo.
* **Visualização Paralela com XDMF**:
    * `xdmf_output`: Um contexto que gerencia a criação de arquivos `.xmf`, que descrevem como os dados HDF5 devem ser interpretados por ferramentas como o ParaView. O fluxo inclui `xdmf_init` , `xdmf_register_cell_property`  e `xdmf_write_timestep`.

---
#### **Serialização: `higtree-serialize.h`, `higtree-serialize.c`**
Fornece funções para converter a estrutura de dados da HiG-Tree em um formato
de array plano (`hig_serial_tree`) e vice-versa. Isso é essencial para enviar
árvores ou sub-árvores entre processos MPI.

* `hs_serialize`: Serializa uma HiG-Tree em um buffer.
* `hs_deserialize`: Reconstrói uma HiG-Tree a partir de um buffer serializado.

---
### Módulos Utilitários
---
* **`allocator.h`, `allocator.c`**: Um alocador de memória simples que permite a desalocação em massa de todos os buffers alocados de uma só vez (`allocator_destroy`).
* **`coord.h`, `coord.c`**: Contém a definição do tipo `Point` e um grande conjunto de macros para operações vetoriais eficientes (e.g., `POINT_ADD`, `POINT_SUB`, `POINT_MULT_SCALAR`).
* **`mapper.h`, `mapper.c`**: Fornece a estrutura `mp_mapper`, um hash map para mapear IDs únicos de células/faces para índices de vetores locais (`mp_lookup`, `mp_assign`).
* **`point-cloud.h`, `point-cloud.c`**: Uma estrutura de dados para gerenciar uma coleção de pontos.
* **`point-mapper.h`, `point-mapper.c`**: Um hash map especializado para mapear coordenadas de pontos (`Point`) a valores, tratando de problemas de precisão de ponto flutuante.
* **`rect.h`, `rect.c`**: Funções para operações com retângulos
(`Rect`), como `rect_intersect` e `rect_contains`.
* **`uniqueid.h`, `uniqueid.c`**: Gera IDs únicos para os elementos
da malha.
* **`utils.h`, `utils.c`**: Funções e macros de utilidade geral, como alocação
de memória com verificação de erro (`ALLOC`) , constantes de precisão
(`EPSDELTA`, `EPSMACH`), e a inicialização da biblioteca
(`higtree_initialize`).
* **`dim.h`, `types.h`, `Debug-c.h`, `rng.h`, `poisson.h`**: Arquivos de
cabeçalho que definem, respectivamente, a dimensionalidade do problema (`DIM`),
tipos de dados básicos (`real`), macros de depuração, um gerador de números
aleatórios (`rng.h`), e a interface para um possível solver de Poisson
(`poisson.h`).

# Documentação do Código da Biblioteca HiG-Tree

Esta é uma documentação detalhada da biblioteca de software fornecida, que
parece ser um framework para simulações numéricas baseadas em malhas
cartesianas hierárquicas. A biblioteca é escrita em C e utiliza MPI para
paralelismo.

### 1. Visão Geral
A biblioteca fornece um conjunto de ferramentas para realizar simulações numéricas em domínios complexos, que são discretizados usando uma estrutura de dados de malha cartesiana hierárquica chamada **HiG-Tree**. Ela foi projetada com o paralelismo em mente, utilizando o padrão MPI para distribuição de dados e computação.

**Principais Características:**
* **Malha Hierárquica (HiG-Tree):** Permite o refinamento adaptativo da malha, concentrando a resolução computacional apenas onde é necessário, economizando memória e tempo de processamento.
* **Decomposição de Domínio e Balanceamento de Carga:** Utiliza a biblioteca Zoltan para particionar a malha de forma eficiente entre múltiplos processos, garantindo que a carga de trabalho seja distribuída de maneira equilibrada.
* **Cálculo de Stencil de Alto Nível:** Abstrai o complexo processo de cálculo de stencils para métodos de diferenças finitas ou volumes finitos. Utiliza o método de **Mínimos Quadrados Ponderados (WLS)** para gerar stencils de alta ordem em pontos arbitrários, mesmo em malhas não uniformes e perto de fronteiras.
* **Interface Abstrata de Solver Linear:** Fornece uma interface unificada para interagir com vários solucionadores de sistemas lineares de alto desempenho, como PETSc e HYPRE, permitindo flexibilidade na escolha do algoritmo de solução.
* **Entrada/Saída Paralela:** Inclui funcionalidades robustas para salvar e carregar dados de simulação em paralelo usando o formato HDF5, além de gerar arquivos de visualização em formatos como VTK e XDMF para análise de dados com ferramentas como ParaView e VisIt.
* **Gerenciamento de Domínio e Condições de Contorno:** Oferece estruturas claras para definir o domínio computacional (composto por múltiplas HiG-Trees), incluindo as "franjas" (células fantasmas) para comunicação entre processos e a aplicação de condições de contorno como Dirichlet e Neumann.

---

### 2. Estruturas de Dados Centrais

#### 2.1. HiG-Tree (`higtree.h`): A Malha Hierárquica
A HiG-Tree é a estrutura de dados fundamental que representa a malha. Ela é uma árvore onde cada nó (`hig_cell`) representa uma região retangular do espaço (célula).

* **`hig_cell`**: Uma célula que pode ser uma folha (a menor resolução em uma região) ou um nó interno que é subdividido em células filhas. Cada célula armazena seus limites (`lowpoint`, `highpoint`), seu pai, e ponteiros para seus filhos.
* **`hig_facet`**: Representa uma face (ou aresta em 2D, lado em 1D) de uma `hig_cell`.
* **Refinamento**: Uma célula pode ser refinada usando `hig_refine_uniform`, que a subdivide em um grid regular de células filhas.
* **Consultas**: A biblioteca oferece funções para consultar propriedades da malha, como encontrar a célula que contém um determinado ponto (`hig_get_cell_with_point`), obter o centro de uma célula (`hig_get_center`), ou encontrar células vizinhas (`hig_get_neighbour`).

#### 2.2. Domínio de Simulação (`domain.h`): Gerenciamento da Malha e Condições de Contorno
Esta estrutura organiza as HiG-Trees e as condições de contorno associadas.

* **`sim_domain`**: Um contêiner que agrupa um conjunto de HiG-Trees que compõem o domínio computacional. Ele distingue entre árvores que fazem parte do domínio local e árvores que representam a "franja" (células fantasmas de processos vizinhos).
* **`sim_boundary`**: Representa uma condição de contorno (CC). Cada CC é associada a uma HiG-Tree que define sua geometria e armazena os valores da condição (ex: valores de Dirichlet ou gradientes de Neumann). Funções como `sb_create` e `sd_add_boundary` permitem a criação e associação de CCs ao domínio.
* **`sim_facet_domain`**: Uma especialização do `sim_domain` projetada para trabalhar com dados centrados nas faces das células, como fluxos ou velocidades.

---

### 3. Paralelismo e Decomposição de Domínio

A biblioteca é projetada para rodar em clusters de computadores, gerenciando a distribuição da malha e a comunicação entre processos.

#### 3.1. Balanceador de Carga (`lbal.h`)
O `load_balancer` é responsável por particionar o domínio global entre os processos MPI.

* **Processo de Particionamento**:
    1.  O processo mestre (rank 0) normalmente cria a malha global.
    2.  Todos os processos criam um contexto `load_balancer`.
    3.  As árvores a serem particionadas são adicionadas ao balanceador com `lb_add_input_tree`.
    4.  A função coletiva `lb_calc_partition` é chamada. Internamente, ela usa a biblioteca Zoltan para calcular uma distribuição de células que minimiza a comunicação e equilibra a carga de trabalho.
    5.  Após a partição, cada processo pode obter suas árvores de domínio local e as árvores de franja (`lb_get_local_tree`).

#### 3.2. Domínio Particionado (`pdomain.h`)
Após o balanceamento de carga, o domínio de cada processo é representado por estruturas `psim_*`.

* **`partition_graph`**: Armazena o resultado do particionamento, descrevendo quais processos são vizinhos e as regiões de fronteira entre eles.
* **`psim_domain`**: Representa a visão de um processo sobre o domínio distribuído. Ele contém o `sim_domain` local (com domínio e franjas) e o `partition_graph` para entender a topologia da vizinhança.
* **`distributed_property`**: Uma abstração crucial para dados de simulação (como pressão, temperatura, etc.). É um vetor distribuído que armazena os valores para as células locais e as células de franja. A função `dp_sync` é usada para atualizar os valores nas células de franja, comunicando-se com os processos vizinhos.

---

### 4. Métodos Numéricos

#### 4.1. Stencils (`domain.h`, `global-stencil.h`)
O stencil é uma representação da dependência de um ponto em relação aos seus vizinhos, fundamental para a discretização de equações diferenciais.

* **`sim_stencil`**: Uma estrutura que armazena os coeficientes do stencil como uma lista de pares `(ID_vizinho, peso)` e um termo independente (`rhs`).
* **Cálculo de Stencil**: As funções `sd_get_stencil` (para dados em células) e `sfd_get_stencil` (para dados em faces) são o coração do motor numérico. Elas calculam os pesos de interpolação para um ponto arbitrário `x`, encontrando as células/faces vizinhas mais relevantes e usando WLS para determinar os pesos. Elas tratam automaticamente de casos complexos, como pontos próximos a fronteiras ou fora do domínio.

#### 4.2. Mínimos Quadrados Ponderados (WLS) (`wls.h`)
O WLS é o método matemático usado para garantir a precisão da interpolação em malhas não uniformes.

* **`wls_interpolator`**: Implementa o algoritmo WLS. Ele constrói um pequeno sistema linear localmente, baseado em uma base polinomial e nos pontos vizinhos, e o resolve para encontrar os pesos do stencil. A função `wls_set_samples_and_calc` é tipicamente chamada internamente pelas rotinas de `get_stencil`.

---

### 5. Solvers Lineares (`solver.h`)
Após a discretização, a simulação geralmente resulta em um grande sistema de equações lineares. A biblioteca abstrai a resolução deste sistema.

* **Interface `solver`**: Define uma API genérica para configurar e resolver um sistema `Ax = b`.
* **Funções Principais**:
    * `slv_create`: Cria uma instância de um solver, escolhendo um backend (e.g., PETSc).
    * `slv_set_Ai` / `slv_set_Aij`: Define os coeficientes da matriz `A`.
    * `slv_set_bi`: Define os valores do vetor `b` (lado direito).
    * `slv_assemble`: Finaliza a montagem da matriz e do vetor.
    * `slv_solve`: Resolve o sistema.
    * `slv_get_x`: Obtém o vetor solução `x`.
* **Backends**: A biblioteca pode ser compilada para usar diferentes backends, como PETSc (`solver-petsc.c`), HYPRE (`solver-hypre.c`), ou um solver SOR (Successive Over-Relaxation) interno (`solver-sor.c`). A escolha pode ser feita em tempo de execução através de variáveis de ambiente.

---

### 6. Iteradores (`higtree-iterator.h`)
Para navegar pela complexa estrutura da HiG-Tree, a biblioteca oferece um sistema de iteradores flexível e poderoso.

* **`higcit_celliterator`**: Um iterador para percorrer células (`hig_cell`).
* **`higfit_facetiterator`**: Um iterador para percorrer faces (`hig_facet`).
* **Tipos de Iteradores**: É possível criar iteradores para diferentes propósitos:
    * `higcit_create_all_leaves`: Itera sobre todas as células-folha da árvore.
    * `higcit_create_bounding_box`: Itera sobre as folhas dentro de uma caixa de delimitação.
    * `higcit_create_neighbours`: Itera sobre as células vizinhas a uma dada célula.
    * `higcit_create_concat`: Combina múltiplos iteradores em um só.

---

### 7. Entrada e Saída (`higtree-io.h`)
O módulo de I/O permite a persistência e visualização dos dados da simulação.

* **Formatos Simples**: Funções como `higio_write_to_file` e `higio_read_from_file` salvam e carregam a estrutura da malha em um formato de texto simples.
* **Visualização com VTK**: `higio_print_in_vtk` exporta a malha no formato legado VTK, que pode ser lido por muitas ferramentas de visualização.
* **I/O Paralelo com HDF5**:
    * A biblioteca fornece funções para criar e abrir arquivos HDF5 (`higio_create_hdf5`).
    * `higio_write_tree_hdf5`: Salva a estrutura da malha em HDF5.
    * `higio_write_cell_property_hdf5` e `higio_write_facet_property_hdf5`: Salvam os dados de `distributed_property` em paralelo.
* **Visualização Paralela com XDMF**:
    * O contexto `xdmf_output` gerencia a criação de arquivos `.xmf`, que são arquivos XML que descrevem como os dados armazenados nos arquivos HDF5 devem ser interpretados por ferramentas como o ParaView.
    * Isso permite visualizar resultados de simulações muito grandes que foram salvas em múltiplos arquivos HDF5 (um por processo). O fluxo de trabalho inclui `xdmf_init`, `xdmf_register_cell_property` e `xdmf_write_timestep`.

---

### 8. Módulos Utilitários

* **`allocator.h`**: Um alocador de memória simples que permite a desalocação em massa de todos os buffers alocados de uma só vez (`allocator_destroy`).
* **`coord.h`**: Contém a definição do tipo `Point` e um grande conjunto de macros para operações vetoriais eficientes (e.g., `POINT_ADD`, `POINT_SUB`, `POINT_MULT_SCALAR`).
* **`mapper.h`**: Fornece a estrutura `mp_mapper`, um hash map para mapear IDs únicos de células/faces para índices de vetores locais (`mp_lookup`, `mp_assign`).
* **`point-mapper.h`**: Um hash map especializado para mapear coordenadas de pontos (`Point`) a valores, tratando de problemas de precisão de ponto flutuante.
* **`rect.h`**: Funções para operações com retângulos (`Rect`), como `rect_intersect` e `rect_contains`.
* **`utils.h`**: Contém macros e funções de utilidade geral, como alocação de memória com verificação de erro (`ALLOC`), constantes de precisão (`EPSDELTA`, `EPSMACH`), e a inicialização da biblioteca (`higtree_initialize`).

---
### 9. Apêndice: Referência Rápida da API (Funções Principais por Módulo)

| Módulo (`.h`) | Estrutura/Função Principal | Descrição |
| :--- | :--- | :--- |
| **`higtree.h`** | `hig_cell* hig_create_root(...)` | Cria a raiz de uma nova malha hierárquica. |
| | `hig_refine_uniform(...)` | Refina uma célula em um grid cartesiano uniforme. |
| | `hig_get_cell_with_point(...)`| Encontra a célula-folha que contém um dado ponto. |
| | `hig_destroy(...)` | Libera a memória de uma árvore e seus descendentes. |
| **`domain.h`** | `sim_domain* sd_create(...)` | Cria um novo domínio de simulação. |
| | `sd_add_higtree(...)` | Adiciona uma HiG-Tree ao domínio local. |
| | `sd_add_boundary(...)` | Adiciona uma condição de contorno ao domínio. |
| | `sim_stencil* stn_create()` | Cria uma nova estrutura de stencil. |
| | `sd_get_stencil(...)` | Calcula os pesos do stencil para um ponto no domínio. |
| **`lbal.h`** | `load_balancer* lb_create(...)` | Cria um contexto para balanceamento de carga. |
| | `lb_calc_partition(...)` | Executa o particionamento do domínio entre os processos. |
| | `hig_cell* lb_get_local_tree(...)`| Obtém uma das árvores locais resultantes após o particionamento. |
| **`pdomain.h`** | `psim_domain* psd_create(...)` | Cria um domínio particionado a partir de um domínio local e um grafo de partição. |
| | `distributed_property* psd_create_property(...)`| Cria um vetor de dados distribuído associado ao domínio. |
| | `dp_sync(...)` | Sincroniza os valores das franjas (células fantasmas) de uma propriedade distribuída. |
| | `int psd_lid_to_gid(...)` | Converte um ID local para um ID global único no sistema. |
| **`solver.h`** | `solver* slv_create(...)` | Cria uma instância de um solver linear. |
| | `slv_set_Ai(...)` / `slv_set_bi(...)`| Define os valores da matriz e do vetor do lado direito. |
| | `slv_solve(...)` | Resolve o sistema linear `Ax=b`. |
| **`higtree-io.h`**| `xdmf_output* xdmf_init(...)` | Inicializa o contexto para escrita de arquivos de visualização XDMF/HDF5. |
| | `xdmf_register_cell_property(...)`| Registra uma propriedade para ser salva. |
| | `xdmf_write_timestep(...)` | Escreve os dados do passo de tempo atual para os arquivos. |
