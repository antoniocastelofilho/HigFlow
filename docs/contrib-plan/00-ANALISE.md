# HigFlow - Dossiê de Análise Técnica

> Documento de trabalho interno da contribuição de **Juniormar Organista**.
> Escrito em português por ser material de decisão pessoal. **Não integra nenhum
> Pull Request para o upstream** - todo artefato que sobe para o upstream é em inglês.
>
> Base da análise: commit `f5eb580` (`master`, sincronizada com `upstream/master`).
> Data da análise: 2026-08-20.

---

## 1. O que é o projeto

HigFlow é um simulador de escoamentos incompressíveis desenvolvido no ICMC/USP,
construído sobre duas camadas:

| Camada | Papel | Tamanho |
|---|---|---|
| **higtree** | Malha AMR hierárquica (árvore de células), decomposição de domínio MPI, interpolação por mínimos quadrados móveis (WLS), interface para solvers lineares (PETSc, HYPRE, ViennaCL, SOR) | ~52.000 linhas |
| **higflow** | Solver Navier–Stokes por método de projeção, modelos constitutivos reológicos, VOF multifásico, eletro-osmótica, I/O e VTK | ~124.000 linhas |

A discretização é por diferenças finitas em malha deslocada (*staggered*): pressão e
propriedades escalares no centro da célula, velocidades nas faces.

### 1.1 Modelos físicos implementados

| Família | Modelos |
|---|---|
| Newtoniano | viscosidade constante |
| Newtoniano generalizado | power-law e variantes |
| Viscoelástico diferencial | Oldroyd-B, Giesekus, LPTT, GPTT (Mittag-Leffler), FENE-P, e-FENE |
| Viscoelástico integral | K-BKZ, K-BKZ fracionário (amortecimento PSM/UCM) |
| Viscosidade variável / tixotropia | BMP, BMP-solvent, MBM, NM-taup, NM-t |
| Shear banding | VCM (two-species) |
| Elastoviscoplástico | - |
| Suspensões shear-thickening | - |
| Multifásico | VOF com PLIC, ELVIRA, altura-função adaptativa, curvatura por diferenças finitas |
| Eletro-osmótico | Poisson–Nernst–Planck, Poisson–Boltzmann; acoplado a viscoelástico e a multifásico |

### 1.2 Métodos numéricos

- **Temporal:** Euler explícito, RK2, RK3, Euler semi-implícito, Crank–Nicolson semi-implícito, BDF2 semi-implícito
- **Convectivo:** central, upwind 1ª ordem, upwind 2ª ordem; esquemas de alta resolução CUBISTA, QUICK, Modified Coefficient Upwind
- **Espacial:** 2ª ordem (4ª ordem declarada mas **não implementada**)
- **Projeção:** incremental e não-incremental

---

## 2. Inventário quantitativo

### 2.1 Arquivos

| Métrica | Valor |
|---|---|
| Arquivos versionados | 1.507 |
| Working tree | 99 MB |
| **Diretório `.git`** | **180 MB** |
| Fontes C | 167 arquivos |
| Cabeçalhos C | 163 arquivos |
| Fontes C++ | 4 arquivos (`build-fringe.cpp`, `rng.cpp`, `solver-viennacl.cpp`, `term-det.cpp`) |
| Malhas `.amr` | 780 arquivos |
| Configurações YAML | 56 arquivos |
| Scripts shell | 35 arquivos |

### 2.2 Os dez maiores arquivos-fonte

| Linhas | Arquivo |
|---|---|
| 10.187 | `higflow/src/hig-flow-io.c` |
| 3.684 | `higflow/src/hig-flow-step-shear-thickening-suspension.c` |
| 3.380 | `higflow/src/hig-flow-step-viscoelastic-shear-banding.c` |
| 2.657 | `higflow/src/hig-flow-step-viscoelastic-variable-viscosity.c` |
| 2.556 | `higflow/src/hig-flow-step-multiphase.c` |
| 2.349 | `higflow/src/hig-flow-step-viscoelastic-integral.c` |
| 2.163 | `higflow/src/hig-flow-step-viscoelastic.c` |
| 2.070 | `higtree/src/domain.c` |
| 1.995 | `higflow/src/hig-flow-step-electroosmotic.c` |
| 1.957 | `higtree/src/lbal.c` |

### 2.3 Contribuidores

| Commits | Autor |
|---|---|
| 69 | Pedro Coimbra |
| **49** | **juniormar** *(mesma pessoa)* |
| 40 | Daniel Garcia |
| **38** | **Juniormar Organista** *(mesma pessoa)* |
| 23 | kainaas |
| **13** | **juniormarorganista** *(mesma pessoa)* |
| 5 | Kainã |
| 4 | Johnatas |
| 3 | Antonio Castelo Filho |
| **1** | **Juniormar Organista** *(quarto e-mail)* |

**Suas 101 contribuições estão fragmentadas em 4 identidades git.** Um `.mailmap`
consolida isso e faz o GitHub e qualquer ferramenta de estatística atribuírem tudo a
você corretamente. É a correção de maior retorno por esforço de todo este dossiê.

---

## 3. Estado do upstream - trabalho paralelo não merjado

`master` está sincronizada com `upstream/master`, mas existem **quatro branches ativas**
no upstream que nunca foram integradas:

| Branch | Data | Commits à frente | Conteúdo |
|---|---|---|---|
| `Kaina` | **2026-08-05** | 72 | Doxygen, 9 tutoriais com imagens, grafo de dependências, documentação do formato `.amr`, `nn-weights.cpp` (inferência MLP), modelos `.pt`, remoção da dependência LibTorch |
| `PC_Daniel_Mesh` | 2026-06-28 | 75 | Correção de conservação de massa em VOF 3D multifásico |
| `Daniel` | 2025-11-04 | 7 | VOF/PLIC, `hig-flow-step`, `hig-flow-terms` |
| `PC_ImproveDocumentation` | 2025-10-31 | - | `conteiner/Dockerfile`, `conteiner/Dockerfile.petsc`, `install_higflow_ubuntu22.sh`, `install_higflow_arch.sh`, READMEs de `higflow/` e `higtree/src/` |
| `Castelo` | 2023-05-12 | - | Ajustes de Makefile |
| `VersionUpdateOnGithub` | 2021-06-21 | - | Energia cinética |
| `new_order_multiphase_functions` | 2022-03-18 | - | Reordenação multifásico |

**Consequência estratégica:** container, manual de instalação e ML em C++ já têm
trabalho iniciado. A decisão tomada foi **construir em cima dando crédito**, corrigindo
os defeitos reais do material existente. Isso posiciona a contribuição como
"consertei e completei", não "refiz por cima" - o que é decisivo num repositório
acadêmico com autores ativos.

### 3.1 Defeitos no material de container existente (`PC_ImproveDocumentation`)

Auditoria de `conteiner/Dockerfile.petsc`:

| # | Defeito | Severidade |
|---|---|---|
| 1 | `echo '. $HOME/.varsrc"' >> $HOME/.bashrc` - **aspa desbalanceada** gravada no `.bashrc`; todo shell subsequente do container falha ao interpretar | **Crítico** |
| 2 | `--with-debubbing=yes` - typo de `--with-debugging`; o `configure` do PETSc aborta em opção desconhecida | **Crítico** |
| 3 | Cabeçalho diz `# Dockerfile to build Ubuntu22.04x64 + OpenFOAM-9 + Python 3.10` - copiado de outro projeto, não descreve a imagem | Cosmético |
| 4 | `FROM ubuntu_petsc3.14:v01` - tag manual, não é multi-stage; exige o usuário construir e taguear na ordem certa sem que nada documente isso | Alto |
| 5 | `chmod 777 libfyaml-master` | Médio |
| 6 | `COPY . "$HOME/HigFlow/"` sem `.dockerignore` - o tarball de 37 MB e os binários de 40 MB entram na camada permanentemente | Alto |
| 7 | `ENV HOME "/home/hig_user/"` com barra final para caminhos viram `/home/hig_user//HigFlow` | Cosmético |
| 8 | Nenhuma versão fixada (`apt-get install` sem pin, `FROM ubuntu:22.04` sem digest) - build não reprodutível | Médio |

---

## 4. Defeitos confirmados no código-base atual

Cada item abaixo foi verificado diretamente no fonte, com arquivo e linha.

### 4.1 Sistema de build

| # | Local | Defeito |
|---|---|---|
| B1 | `CMakeLists.txt` | **Não existe comando `project()`.** CMake opera em modo degradado |
| B2 | `CMakeLists.txt:11` | `if(${debug})` - sem `-Ddebug=...` na linha de comando isto expande para `if()`, erro de configuração. O mesmo em `if(NOT ${prefix} STREQUAL "")` |
| B3 | `CMakeLists.txt:5-6` | `-DDIM=${dim}` - sem `-Ddim=N` compila com `-DDIM=` (vazio) |
| B4 | `CMakeLists.txt:29-33` | `find_package(BLAS)` e `find_package(LAPACK)` gravam em `MFSIM_DEPENDENCIES` - **variável de outro projeto (MFSim), jamais usada**. BLAS e LAPACK são localizados e silenciosamente descartados |
| B5 | `CMakeLists.txt:13-14` | `-march=native -mtune=native -flto=8 -fno-strict-aliasing` aparecem apenas em `add_link_options`, **não** em `add_compile_options`. As flags que mais importam em tempo de compilação não são aplicadas |
| B6 | `CMakeLists.txt:110,113` | `include(.../CMakeLists.txt)` em vez de `add_subdirectory()` - anti-padrão; escopo de variáveis vaza entre módulos |
| B7 | `higflow/Makefile:40` | `-ltrinilos_zoltan` - typo de `trilinos`. O `example*/Makefile` grafa correto, o da biblioteca não |
| B8 | `higflow/Makefile:45` vs `:92` | Define `ANLIB = gcc-ranlib` mas invoca `$(RANLIB)`. Usa-se o `ranlib` padrão do sistema, não o do GCC - o índice do arquivo `.a` fica inconsistente com objetos LTO |
| B9 | `higtree/Makefile:76,80` | `build-fringe` listado **duas vezes** em `MODULES` |
| B10 | `higtree/Makefile:31-33` | `-I/usr/include/trilinos`, `-I/usr/include/hypre`, `-I/usr/include/petsc` hardcoded - quebra fora de Debian/Ubuntu |
| B11 | CMake vs Makefile | **Compilam conjuntos diferentes de fontes.** Detalhe em §4.2 |

#### 4.2 Divergência entre os dois sistemas de build

| Arquivo | No `Makefile`? | No `CMakeLists.txt`? |
|---|:---:|:---:|
| `hig-flow-timestep.c` | sim | **não** |
| `hig-flow-step-electroosmotic-viscoelastic.c` | sim | **não** |
| `hig-flow-step-multiphase-electroosmotic.c` | sim | **não** |
| `hig-flow-step-multiphase-electroosmotic-viscoelastic.c` | sim | **não** |
| `hig-flow-vof-elvira.c` | **não** | sim |
| `hig-flow-vof-9-cells.c` | **não** | sim |
| `hig-flow-mittag-leffler.c` | **não** | sim |
| `higtree/src/solver.camila.c` | **não** | sim |
| `higtree/src/solver-petsc.camila.c` | **não** | sim |
| `higtree/src/solver-sor.c` | condicional (`SOR=1`) | sempre |
| `higtree/src/solver-viennacl.cpp` | condicional (`VIENNACL=1`) | sempre |

**Nenhum dos dois sistemas compila o projeto inteiro.** Dois usuários com dois métodos
de build obtêm binários com funcionalidades diferentes. Isso é a raiz de boa parte da
dificuldade relatada de "fazer o projeto funcionar".

Nota adicional: `mittag-leffler` é dependência do modelo GPTT, que é oferecido no YAML
de exemplo - quem compila via `Makefile` não tem esse símbolo.

### 4.3 Instalação e ambiente

| # | Local | Defeito |
|---|---|---|
| I1 | `install_higflow_ubuntu22` | Sem shebang, sem `set -euo pipefail`, sem verificação de erro entre etapas |
| I2 | idem | `--with-debubbing=yes` - typo; PETSc aborta |
| I3 | idem | `--PETSC_ARCH=x86_64` - na sintaxe do PETSc é `PETSC_ARCH=` sem hífens |
| I4 | idem | Instala `openmpi-bin` + `libopenmpi-dev` **e** `mpich` **e** ainda pede `--download-openmpi` ao PETSc: três MPIs concorrentes no mesmo sistema |
| I5 | idem | `export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig/libfyaml.pc` - a variável espera **diretório**, recebe caminho de arquivo |
| I6 | idem | O `export` anterior (`.../bibliotecas/hypre/`) é sobrescrito pelo seguinte, e nenhum dos dois sobrevive ao fim do script (não é `source`) |
| I7 | idem | `sudo ./configure` no PETSc - a árvore de build fica propriedade do root |
| I8 | idem | `sudo chmod 777 libfyaml-master` |
| I9 | idem | `pip3 install` global; em Ubuntu ≥23 falha com `externally-managed-environment` |
| I10 | idem | Cria symlink `libHYPRE.so para libHYPRE_krylov.so`, contornando um problema de linkagem em vez de resolvê-lo |
| I11 | `varsrc` | `PETSC_DIR=$(pwd)/bibliotecas/petsc-3.14.0/x86_64` com `PETSC_ARCH=arch-linux-c-debug` - **contradiz** o instalador, que instala em `/opt/petsc-3.14.0-openmnpi-hypre-hdf5` com `PETSC_ARCH=x86_64`. Seguir o README à risca não produz um ambiente funcional |
| I12 | `varsrc` | `$(pwd)` - só funciona se o `source` for feito da raiz do repositório |
| I13 | `README.md` | Bloco corrompido: `"Após re./configure --prefix=... ; sleep 5alizar um dos passos anterior"` - um comando foi colado no meio da palavra "realizar" |
| I14 | `stacks/singularity/scripts/petsc.sh:31` | `cp configure.log /pacotes/.` - diretório `/pacotes` nunca é criado |
| I15 | `higflow_image.def` | `%post` chama `./scripts/pip.sh` com caminho relativo; o diretório de trabalho em `%post` não é a raiz |

### 4.4 Código - correção

| # | Local | Defeito | Impacto |
|---|---|---|---|
| C1 | `higflow/src/hig-flow-kernel.h:18` | `#define DEBUG` **incondicional** num cabeçalho público | Todo TU que inclui o kernel entra no caminho de debug, independentemente de `NDEBUG`. Anula o modo otimizado |
| C2 | `higtree/src/Debug-c.h:9-13` | `static struct timeb __lasttime; static int __debug_flag[100];` **definidos em cabeçalho** | Uma cópia privada por unidade de tradução. `DEBUG_PUSH`/`DEBUG_POP` numa TU não são vistos por outra - a pilha de debug não funciona como projetada |
| C3 | `higtree/src/Debug-c.h:32` | `#define DEBUG_WARNING(x) {fprintf(debugfd, x);}` | **Format string bug.** Qualquer `%` na mensagem lê a pilha |
| C4 | `higtree/src/Debug-c.h:41` | `DEBUG_ASSERT` termina com `*(int *)NULL = 0;` | Comportamento indefinido deliberado em vez de `abort()`. Com otimização o compilador pode eliminar o caminho inteiro |
| C5 | `higtree/src/Debug-c.h:5` | `#include <sys/timeb.h>` / `ftime()` | Removido do POSIX.1-2008, marcado obsoleto no glibc. **Bloqueia portabilidade** |
| C6 | Todos os `Makefile` e `CMakeLists.txt` | `-Ofast` | Implica `-ffast-math` para `-ffinite-math-only`, que autoriza o compilador a **assumir que Inf e NaN nunca ocorrem**. Ver §4.5 |
| C7 | `hig-flow-step-electroosmotic.c:1809` etc. | `real u_min[DIM]={INFINITY}` | Inicializa **apenas o elemento 0**; os demais recebem `0.0`. Atualmente as linhas de uso estão comentadas (defeito dormente), mas o idioma está presente em 6 declarações |
| C8 | `hig-flow-io.c:6114-7696` | `higflow_load_all_controllers_and_parameters_yaml()` - **~1.580 linhas numa única função** | Intestável, irrevisável |
| C9 | Global | 308 chamadas a `exit()` em biblioteca | Uma biblioteca não deve terminar o processo do chamador. Impede tratamento de erro e testes unitários |
| C10 | Global | 306 chamadas a `fopen()` | Auditoria de verificação de retorno pendente |
| C11 | `higflow/include/`, `higtree/include/` | Cópias geradas por `cp src/*.h` **versionadas** | `higtree/include/domain.h` e `solver-petsc.h` **já divergiram** do fonte. `higflow/include/hig-flow-step-multifase.h` é órfão (grafia antiga em português) |

### 4.5 A questão do `-Ofast` - análise numérica

`-Ofast` = `-O3` + `-ffast-math` + `-fallow-store-data-races` + relaxamentos adicionais.
De `-ffast-math` decorre `-ffinite-math-only`, que informa ao compilador que nenhum
valor de ponto flutuante será `Inf` ou `NaN`.

O código **contradiz essa premissa diretamente**. Em
`hig-flow-step-electroosmotic.c:923` e `:932`, e em
`hig-flow-step-multiphase-electroosmotic.c:561,569`:

```c
real max_psi_res_global = INFINITY;
// ... laço de convergência do tipo:  while (max_psi_res_global > tol) { ... }
```

`INFINITY` é usado como sentinela para garantir que o laço execute ao menos uma
iteração. Sob `-ffinite-math-only` o compilador está autorizado a tratar a comparação
como se o operando não pudesse ser infinito - o laço pode ser reordenado ou eliminado.
O comportamento passa a depender de versão do compilador e nível de inline.

Consequências, em ordem de gravidade:

1. **Detecção de divergência quebra.** Numa simulação que diverge, o valor vira `NaN`.
   Sob `-ffinite-math-only`, comparações que detectariam `NaN` são otimizáveis para
   `false`. A simulação continua produzindo lixo silenciosamente em vez de abortar.
2. **Sentinelas `INFINITY` perdem garantia semântica** (caso concreto acima).
3. **Reprodutibilidade destruída.** `-march=native` faz o binário depender da
   microarquitetura da máquina. Dois nós de um mesmo cluster com CPUs diferentes
   produzem resultados diferentes bit a bit. Para um código que gera resultado
   publicável, isso é um problema metodológico, não só de engenharia.
4. **Reassociação de ponto flutuante** (`-funsafe-math-optimizations`) altera a ordem
   de somatórios. Em reduções sobre a malha isso muda o erro de arredondamento
   acumulado de forma não controlada.

O tratamento correto: `-O3` como padrão; `-ffast-math` apenas sob opção explícita e
documentada; `-march=native` nunca como padrão, e sim `-mtune=generic` com
`-march` configurável; e uma opção de build reproduzível com
`-ffp-contract=off`.

### 4.6 Código - desempenho

| # | Achado | Detalhe |
|---|---|---|
| P1 | **292 chamadas a `pow()`**, das quais 83 com expoente inteiro literal | 47× `pow(x,2)`, 35× `pow(x,2.0)`, 1× `pow(x,3.0)`, 3× `pow(x,4.0)`. Em caminho quente, `x*x` é ordens de grandeza mais rápido. Sob `-Ofast` o GCC converte parte, mas isso deixa a otimização refém de uma flag que precisa ser removida (§4.5) |
| P2 | 81 literais de tolerância hardcoded | 42× `1.0e-14`, 18× `1.0e-3`, 6× `1.0e-4`, 4× `1.0e-10`, 3× `1.0e-6`, 3× `1.0e-16` - dispersos, não configuráveis, não documentados |
| P3 | Compilação dupla obrigatória | `-DDIM=2` e `-DDIM=3` geram bibliotecas separadas (`libhig2d.a`, `libhig3d.a`). A dimensão é constante de compilação, não parâmetro de execução |
| P4 | `hig-flow-io.c` com 10.187 linhas | Recompilação integral a cada alteração de I/O |

### 4.7 Higiene do repositório

| Item | Tamanho | Situação |
|---|---|---|
| `bibliotecas/petsc-3.14.0.tar.gz` | **37 MB** | Versionado. `bibliotecas/` está no `.gitignore` - adicionado **depois** do commit, portanto sem efeito |
| `bibliotecas/libfyaml-master.zip` | 451 KB | Idem |
| 11 binários ELF `ns-example`/`ns-complex-3d` | **~40 MB** | Executáveis compilados, versionados |
| `higflow/src.zip` | 779 KB | Snapshot do próprio diretório, versionado |
| `higflow/src/.hig-flow-kernel.h.swp` | 16 KB | Arquivo de swap do vim |
| 7 arquivos `*_old.c` / `*.old.c` | 689 KB | Código morto |
| `higflow/src/src_hugo/` | ~2,4 MB | Cópia paralela de toda a árvore, incluindo `Modifications of files viscoelastic flows with variable viscosity/` - **nome de diretório com espaços**, que quebra Makefiles e scripts ingênuos |
| `higflow/src/contr.flowtype` | 266 KB | Dados no diretório de fontes |
| `Attic/hig-flow-kernel.h` | 24 KB | Cabeçalho antigo |
| No histórico | 45 MB | `example2d_ElectroOsmotic/animation.avi` |
| No histórico | 38 MB + 20 MB | VTKs de saída de VOF |
| No histórico | ~30 MB | Árvore completa de `atf-0.15/` com binários de teste |

`.gitignore` ignora `*.txt`, `*.dat` e `*.vtk` globalmente - o que também mascara
arquivos legítimos e é a razão de 35 `.dat` e 16 `.vtk` estarem versionados
(entraram antes da regra).

### 4.8 Ausências de projeto

| Arquivo | Situação |
|---|---|
| `LICENSE` | **Ausente** - status jurídico das contribuições é indefinido |
| `CONTRIBUTING.md` | Ausente |
| `CODE_OF_CONDUCT.md` | Ausente |
| `CITATION.cff` | Ausente - impede citação automática do software |
| `CHANGELOG.md` | Ausente |
| `.editorconfig` | Ausente |
| `.clang-format` | Ausente |
| `.mailmap` | Ausente |
| `Doxyfile` na raiz | Ausente (existe `higtree/doc/doxygen/HiGTreeDoxy`, não referenciado) |
| CI (qualquer provedor) | **Ausente** |
| Testes automatizados | `higtree/atf-tests/` usa o framework ATF, sem runner integrado nem execução documentada |

---

## 5. O atrito de uso - medição

Este é o problema central relatado: "tem que gerar muita coisa para poder rodar um caso".

### 5.1 O que é necessário hoje para um único caso

| Artefato | Volume |
|---|---|
| `input/<caso>.load.par.contr.yaml` | **261 linhas** - contém os parâmetros de **todos** os modelos, mesmo os não usados |
| `input/<caso>.load.bc.yaml` | **150 linhas** |
| `input/<caso>.load.domain.yaml` | 5 linhas |
| `input/<caso>.load.init.yaml` | 6 linhas |
| **Subtotal YAML** | **422 linhas** |
| Driver `ns-example-2d.c` | **280 a 2.609 linhas**, conforme o caso |
| Cabeçalho `ns-example-2d.h` | ~60 linhas |
| `Makefile` | ~95 linhas |
| `CMakeLists.txt` | ~20 linhas |
| Malhas `amrs/domain/*.amr` | 1 por bloco de domínio |
| Malhas `amrs/bc/*.amr` | 1 por fronteira (4 no caso mais simples) |
| **Compilação** | Obrigatória a cada alteração de fronteira ou condição inicial |

Tamanho dos drivers por exemplo:

| Exemplo | Linhas do driver |
|---|---|
| `example3d_complex` | 2.609 |
| `example2d_BMP` | 929 |
| `example2d_VOF_Gptt` | 649 |
| `example2d_VOF_Oldroyd` | 648 |
| `example2d_VOF` | 616 |
| `example2d_Oldroyd` | 505 |
| `example3d_lid_driven` | 492 |
| `example2d_KBKZ` | 445 |
| `example2d_Newt_contraction` | 329 |
| `example2d_Gptt` | 326 |
| `example2d_Newt` | 280 |

### 5.2 Onde está a duplicação

Os 11 drivers repetem:

- O mesmo `main()`: `higflow_initialize` para `higflow_create` para carregar YAML para registrar
  funções externas para criar domínio para criar solver para laço temporal para destruir.
  Cerca de 120 linhas idênticas em todos.
- As mesmas 10 funções `get_*`: `get_pressure`, `get_velocity`, `get_source_term`,
  `get_facet_source_term`, `get_viscosity`, `get_boundary_pressure`,
  `get_boundary_velocity`, `get_boundary_source_term`,
  `get_boundary_facet_source_term`, `get_boundary_viscosity`.
- Em quase todos, essas funções apenas retornam constantes. O caso Newtoniano tem
  **uma única linha de física real** no driver inteiro:
  `value = 1.5*(1.0 - center[1]*center[1]);` (perfil de Poiseuille na entrada).

**Uma linha de física exige 280 linhas de código e recompilação.** Esse é o
diagnóstico preciso do problema.

### 5.3 Formato `.amr`

Formato posicional sem cabeçalho, 4 linhas:

```
0.0 8.0 -1.0 1.0      # xmin xmax ymin ymax
1                     # número de níveis de refinamento
0.05 0.05 1           # dx dy (ou espaçamento por nível)
1 1 160 40            # índices e número de células por direção
```

Sem validação, sem mensagem de erro compreensível, sem gerador. Um domínio 3D
complexo exige dezenas desses arquivos escritos à mão.

---

## 6. Portabilidade para Windows

### 6.1 O código em si é razoavelmente portável

| Categoria | Ocorrências |
|---|---|
| Cabeçalhos exclusivos de POSIX | 3 (`unistd.h`, `sys/timeb.h`, `sys/resource.h`) |
| `__attribute__` do GCC | 4 |
| Funções aninhadas (extensão GCC) | 0 |
| `system()` / `fork()` / `popen()` | 0 |
| VLAs (arrays de tamanho variável) | 71 - sobretudo em **parâmetros** de função |

### 6.2 O que realmente bloqueia

O obstáculo não é o código: é a pilha de dependências.

| Dependência | Windows nativo |
|---|---|
| PETSc | Possível via MSYS2, laborioso |
| OpenMPI | Não suportado; MS-MPI é a alternativa |
| Zoltan (Trilinos) | Difícil |
| HDF5 | Disponível |
| glib-2.0 | Via MSYS2/vcpkg |
| libfyaml | Sem porte oficial |
| ViennaCL | Header-only, OK |
| **libnuma (`numactl`)** | **Exclusivo de Linux. Sem equivalente.** Exigido pelo `CMakeLists.txt` (`PKG_CHECK_MODULES(NUMACTL REQUIRED numa)`) e usado por `higtree/src/solver-sor.c` |

`libnuma` é `REQUIRED` no CMake, portanto **a configuração falha em Windows antes de
compilar uma única linha** - mesmo que todo o resto estivesse disponível.

### 6.3 Conclusão

O caminho correto para Windows é **WSL2 e containers**, não porte nativo. Um porte
MSYS2 seria meses de trabalho para um resultado frágil e de manutenção cara. As duas
melhorias que valem a pena no código são:

1. Tornar `numactl` opcional no CMake (`solver-sor` já é opcional no Makefile) - o que
   também beneficia macOS e clusters sem libnuma.
2. Substituir `ftime()` por `clock_gettime(CLOCK_MONOTONIC, ...)`, que é POSIX moderno
   e tem equivalente trivial em qualquer plataforma.

### 6.4 Estado da máquina de desenvolvimento

| Ferramenta | Situação |
|---|---|
| Docker | **Não instalado** |
| WSL2 | Kernel presente, **nenhuma distribuição instalada** |
| `make` | **Não disponível** |
| CMake 3.x | Presente (Strawberry Perl) |
| GCC / G++ | Presente (MinGW-w64, via Strawberry Perl) |
| Python 3.14 | Presente |
| Git | Presente |
| ParaView | Não instalado |

Nada do projeto pode ser compilado ou executado localmente hoje. Resolver isto é
pré-requisito de qualquer etapa que envolva execução ou geração de figuras.

---

## 7. Migração C para C++ - avaliação de viabilidade

A conversão é **substancialmente mais tratável do que o tamanho do código sugere**.

| Incompatibilidade CparaC++ | Ocorrências | Esforço |
|---|---|---|
| Palavras reservadas de C++ usadas como identificador | **1** (`struct _local_neighbor *new` em `lbal.c:1189`) | Trivial |
| `malloc`/`calloc`/`realloc` sem cast explícito | 42 | Mecânico |
| VLAs | 71 (majoritariamente parâmetros) | Moderado - `std::vector` ou `span` |
| Literais compostos | 0 | - |
| `restrict` | 0 | - |
| `_Generic`, `_Static_assert`, `_Atomic` | 0 | - |
| Inicializadores designados | 1 | Trivial |

Ou seja: o **Nível 1** (fazer tudo compilar como C++ sem mudar semântica) é uma tarefa
mecânica de dias, não de meses. O trabalho real está nos níveis seguintes:

- **Nível 2 - RAII:** eliminar as 308 chamadas a `exit()`, encapsular `FILE*`,
  ponteiros de MPI, objetos PETSc e alocações em tipos com destrutor. Elimina
  vazamentos e torna o código testável.
- **Nível 3 - tipos e templates:** `DIM` como parâmetro de template em vez de macro
  de compilação, eliminando a necessidade de duas bibliotecas separadas; tipos fortes
  para `Point`, tensores e índices, hoje todos `real[]` ou `int`.

Riscos a controlar em qualquer nível: **regressão numérica silenciosa**. Nenhuma etapa
de C++ deve começar antes de existir a suíte de verificação numérica (Etapa 08 do
roadmap), que é o único mecanismo capaz de provar que a conversão não alterou
resultados.

---

## 8. Síntese - prioridades por impacto

| Prioridade | Item | Justificativa |
|---|---|---|
| 1 | `.mailmap` | Consolida 101 commits fragmentados em 4 identidades. Custo: 5 linhas |
| 2 | Unificar os dois sistemas de build | Raiz da maior parte do "não consigo compilar" |
| 3 | Container funcional + manual | Elimina a barreira de entrada por completo |
| 4 | README com resultados visíveis | Um simulador de reologia sem imagens é invisível |
| 5 | Remover `-Ofast`, tornar `-march` configurável | Correção metodológica, não só de engenharia |
| 6 | Runner genérico com expressões | Transforma 700+ linhas por caso em ~30 |
| 7 | Suíte de verificação numérica + CI | Pré-requisito de qualquer refatoração segura |
| 8 | Limpeza do repositório | 78 MB de artefatos versionados |
| 9 | Correções pontuais (C1-C11) | Baixo custo, alto valor de revisão |
| 10 | Documentação, tradução, Doxygen | Alcance internacional |
| 11 | Migração C++ | Depende de 7 estar pronto |
| 12 | Módulo ML | Exploratório |

---

*Fim do dossiê. O plano de execução está em [`01-ROADMAP.md`](01-ROADMAP.md).*
