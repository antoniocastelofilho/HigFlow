# O t8code na suíte da HiGTree

Este diretório contém o **adaptador de malha do t8code** e a sonda de
viabilidade. Nada aqui é construído por omissão: o t8code não é dependência da
HiGTree, e sem `--t8code` a suíte roda exatamente como antes.

Com o adaptador, as 18 cláusulas do contrato de Mesh
(`higtree/src/hig-mesh-contract.h`) passam a ser verificadas para **duas**
implementações de malha, não uma.

```
sem --t8code    104 casos, contrato 18/18
com --t8code    142 casos, contrato 18/18
```

## Construir o t8code

Não está empacotado em nenhuma distribuição — tem de vir do fonte. O
`libp4est-dev` existe no Debian/Ubuntu (2.3.6) mas não basta; o CMake do t8code
busca `sc` e `p4est` por `FetchContent`, então basta ter rede na configuração.

```bash
cd bibliotecas
git clone https://github.com/DLR-AMR/t8code.git
cd t8code
git checkout v4.0.0-26.09
cmake -B build -S . \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/install" \
  -DT8CODE_ENABLE_MPI=ON \
  -DT8CODE_BUILD_AS_SHARED_LIBRARY=ON \
  -DT8CODE_BUILD_TESTS=OFF -DT8CODE_BUILD_EXAMPLES=OFF \
  -DT8CODE_BUILD_BENCHMARKS=OFF -DT8CODE_BUILD_TUTORIALS=OFF
cmake --build build -j4
cmake --install build
```

**Não clone com `--depth 1`.** O CMake deriva a versão de `git describe`, e sem
as tags ele para com `VERSION ".." format invalid` — uma mensagem que não
menciona tag nenhuma. Se já clonou raso, `git fetch --unshallow --tags` resolve.

A versão em uso aqui é a tag `v4.0.0-26.09`; o `git describe` a reporta como
`v4.0.0-26.08-manually_fix_version`, que é o nome que o CMake grava na
biblioteca. Não é erro.

O diretório `bibliotecas/` é ignorado pelo git, então a biblioteca não entra na
árvore.

## Rodar a suíte com ele

```bash
set -a; . ./varsrc; set +a
python3 ci/run_higtree_tests.py --t8code bibliotecas/t8code/install
```

Carregue o `varsrc` **da raiz do repositório**: ele usa `$(pwd)`, e carregado de
outro diretório aponta `PETSC_DIR` para um `bibliotecas/` que não existe ali. O
sintoma é `undefined reference to PetscFinalize` ao ligar, que aponta para o
PETSc e não para o diretório errado.

Para compilar um teste à mão, sem o driver, exporte também `HIGTREE_DIR` — o
`varsrc` não o define, e sem ele `LIBPATH` fica vazio:

```bash
export HIGTREE_DIR="$PWD/higtree"
make -C higtree/tests DIM=2 T8CODE="$PWD/bibliotecas/t8code/install"
```

## Três restrições de integração

Cada uma custou um ciclo de diagnóstico, e nenhuma se anuncia com clareza.

1. **O t8code exige C++20**; a HiGTree constrói com `-std=gnu++17`. Por isso
   cada unidade aqui compila em C++20 separada e só assinaturas `extern "C"`
   atravessam para o lado C.

2. **O consumidor tem de definir `-DT8_ENABLE_MPI=1 -DT8_ENABLE_MPIIO=1`.** Sem
   isso o `t8.h` para com `MPI configured differently in t8code and libsc` — que
   aponta para o libsc, e não para a definição que falta. O Makefile já as passa.

3. **`t8_cmesh_init()` é obrigatório antes de `t8_cmesh_new_hypercube`.** Sem
   ele o segfault sai três quadros abaixo, dentro de `t8_cmesh_set_tree_class`,
   sem mencionar inicialização nenhuma.

## Os módulos, e por que a costura muda de forma

Uma segunda implementação de malha não herda de uma classe: ela **produz** as
estruturas que as consultas leem, ou **responde** às consultas. Qual das duas
depende do que se está julgando — e essa é a decisão mais fácil de errar aqui.

Uma costura de **produtor** (o t8code entrega uma malha e o MTree responde tudo)
serve para C11, que é sobre o que a malha *representa*. Não serve para nenhuma
cláusula sobre quem *responde*: o produtor materializa uma árvore `hig_cell`, e
aí quem responde é o MTree. Verde assim não diz nada sobre o t8code.

| módulo | costura | cláusulas |
|---|---|---|
| `t8-mesh-producer` | produtor de malha | C11 |
| `t8-point-locator` | localizador | C7, C8, C9 |
| `t8-stencil-support` | suporte de estêncil | C10 |
| `t8-boundary-faces` | faces de contorno | C12 |
| `t8-point-class` | classificador de ponto | C14 |
| `t8-partition` | partição com ghost | C13, P1–P4 |
| `t8-forest-cache` | floresta comum aos demais | — |
| `probe-c11.c` | sonda de viabilidade, fora da suíte | — |

O `t8-mesh-producer` tem **dois caminhos**, e a diferença entre eles é o ponto
central da integração:

- `t8_produz_malha_nao_graduada` materializa a floresta como octree de
  ponteiros — `O(folhas)` por produção. É tradução, e é o custo que se queria
  evitar.
- `t8_preenche_instantaneo` preenche o arranjo plano de
  `higtree/src/hig-mesh-snapshot.h` **direto da floresta**, sem construir árvore
  alguma. É a fronteira das consultas, movida.

## Armadilhas do t8code que os testes já registram

- `t8_forest_element_points_inside` devolve **true para dois vizinhos** quando o
  ponto está na face comum. A primitiva não desempata; a convenção da cláusula
  C8 é imposta sobre ela.
- O índice que `t8_forest_leaf_face_neighbors` devolve é da **floresta local**,
  não da árvore. Tratá-lo como índice de árvore dá um elemento existente, só que
  outro — erro silencioso. Use `t8_forest_get_leaf_element`.
- Sem `do_face_ghost = 1`, "face sem vizinho" deixa de distinguir contorno do
  domínio de fronteira de partição.

## O que ainda não existe

**O t8code em produção.** Os 777 sítios de `higflow/src` seguem chamando as
funções antigas, e nenhum `higflow_solver` produz um instantâneo. Migrá-los
exige antes produzir e guardar um instantâneo por domínio, e invalidá-lo a cada
adaptação de malha — isso é projeto, não varredura.

**Uma decisão de projeto em aberto.** Na interface entre árvores de cmesh, os
dois backends classificam o ponto de forma diferente: a HiGTree diz
`ON_BOUNDARY` (é limite de caixa), o t8code diz dentro (a face tem vizinho). O
caso `criterios_divergem_na_interface_entre_arvores` **registra** a divergência
em vez de eleger um vencedor. Qual dos dois é o desejado muda o que o fechamento
de contorno faz numa interface interna.
