# Diário de execução

> Uma entrada por sessão de trabalho. Documento interno, não integra Pull Requests.
>
> Estado atual verificável a qualquer momento com:
> ```bash
> bash tools/contrib/status.sh --fetch
> ```

---

## Estado das etapas

| Etapa | Branch | Estado | PR |
|---|---|---|---|
| E00 - Ambiente de desenvolvimento | - | **concluída** (WSL2 + Ubuntu 22.04 + Docker Engine 29.7.2) | - |
| E01 - Infraestrutura da contribuição | `juniormar/main` | **concluída** (falta `.mailmap`, que vai no PR de E07) | - |
| E02 - README em inglês | `juniormar/pr2-readme` | **enviada** | [#4](https://github.com/antoniocastelofilho/HigFlow/pull/4) |
| E03 - Instalação Linux | - | pendente | - |
| E04 - Containers | `juniormar/pr3-containers` | **enviada** | [#5](https://github.com/antoniocastelofilho/HigFlow/pull/5) |
| E05 - Instalação Windows | `juniormar/pr4-windows-guide` | **enviada** | [#6](https://github.com/antoniocastelofilho/HigFlow/pull/6) |
| E06 - Galeria de resultados | `juniormar/pr5-gallery` | **enviada** | [#7](https://github.com/antoniocastelofilho/HigFlow/pull/7) |
| E07 - Higiene do repositório | `juniormar/pr1-repo-hygiene` | **enviada** | [#3](https://github.com/antoniocastelofilho/HigFlow/pull/3) |
| E08 - Build unificado | - | pendente | - |
| E09 - Correções pontuais | - | pendente | - |
| E10 - Flags numéricas | - | pendente | - |
| E11 - Verificação numérica | - | pendente | - |
| E12 - Integração contínua | - | pendente | - |
| E13 - Runner: scaffolding | - | pendente | - |
| E14 - Runner: driver unificado | - | pendente | - |
| E15 - Runner: expressões | - | pendente | - |
| E16 - Gerador de malha | - | pendente | - |
| E17 - Doxygen | - | pendente | - |
| E18 - Tradução | - | pendente | - |
| E19 - Governança | - | pendente | - |
| E20 - Divulgação | - | pendente | - |
| E21 - Plano C++ | - | pendente | - |
| E22 - C++ nível 1 | - | pendente | - |
| E23 - C++ nível 2 | - | pendente | - |
| E24 - C++ nível 3 | - | pendente | - |
| E25 - Arquitetura ML | - | pendente | - |
| E26 - PoC ML | - | pendente | - |
| E27 - Galeria de pessoas | dobrada na `pr2-readme` | **estrutura enviada** (aguarda consentimentos) | [#4](https://github.com/antoniocastelofilho/HigFlow/pull/4) |

---

## 2026-08-22 - Reorganização do modelo de branches

O modelo anterior era o inverso do pedido: a `juniormar/main` era um tronco de
integração que nunca virava PR, e as branches de etapa saíam da `master`. O pedido é
que a `juniormar/main` seja a branch principal, base das etapas e origem dos PRs.

A tensão a resolver: se as etapas saem da `main`, tudo que está nela entra no diff de
todo PR, e ali estavam as 2.400 linhas de planejamento em português.

**Decidido:** o planejamento vai para `juniormar/notes`, que nunca é base de nada. A
`juniormar/main` fica limpa.

### O que mudou

| Branch | Antes | Agora |
|---|---|---|
| `juniormar/main` | tronco com planejamento, histórico divergente da cadeia | sai da `master`, cresce por merge de cada etapa, sem planejamento |
| `juniormar/notes` | não existia | só `docs/contrib-plan/` e `tools/contrib/` |
| as cinco de PR | cadeia linear da `master` | **inalteradas** |

### As branches de PR não precisaram mudar

Elas já eram cadeia linear saindo da `master`. Se a `main` começa na `master` e cresce
por merge de cada etapa, "a etapa saiu da main" e "a etapa saiu da anterior" descrevem
o mesmo commit. O modelo novo já estava satisfeito pelas branches; o que estava fora do
lugar era a própria `main`.

Consequência prática: **nenhum force-push, e os PRs #3 a #7 intocados.** Confirmado que
os cinco continuam MERGEABLE e que local e origin apontam para o mesmo commit.

### Verificação

- `juniormar/main` reconstruída merjando pr1 a pr5 em ordem, todas limpas
- Conteúdo da `main` idêntico ao da ponta da cadeia: nenhuma divergência
- 1.430 arquivos versionados, zero de `contrib-plan` ou `tools/contrib`
- `juniormar/notes` com 8 arquivos, só o planejamento
- Histórico das notas preservado: 62 commits, e revisões anteriores ainda têm os
  arquivos de projeto

### Dois tropeços no caminho

O primeiro `git rm` em massa quebrou no caminho `src_hugo/Modifications of files
viscoelastic flows with variable viscosity/`, que tem espaços, e a remoção ficou pela
metade. Refeito limpando o índice inteiro e trazendo de volta só os dois diretórios.

Depois disso o checkout para a `main` abortou: os arquivos de projeto tinham virado
untracked na `notes` e seriam sobrescritos. Antes de limpar, conferi que os 494
untracked estavam todos na árvore da `main`, com uma exceção que era um `.pyc` de
cache.

### Onde fica a ferramenta de status

`tools/contrib/status.sh` vive na `notes`. Para rodar sem trocar de branch:

```bash
git show juniormar/notes:tools/contrib/status.sh | bash
```

---

## 2026-08-22 - Revisão de escrita e envio dos pull requests

### Revisão de escrita

Você removeu os travessões no tronco; restavam outros padrões que a skill de escritor
proíbe e que eu não tinha tratado:

| Padrão | Ocorrências | Destino |
|---|---|---|
| Régua de caixa `─` em comentário | 2.692 | removida; o rótulo vira comentário simples |
| Seta `→` como conectivo | 56 | vira palavra |
| Meia-risca em intervalo | 35 | vira hífen |
| Emoji em tabela | 2 | vira "yes" |
| Travessão nas branches de etapa | 27 a 465 por branch | vira hífen |

Preservado por ser notação, não decoração: a meia-risca entre nomes de pessoas
diferentes (Navier-Stokes, Poisson-Nernst-Planck, Crank-Nicolson), os sinais de
multiplicação e menos, e os glifos que desenham árvore de diretório.

O SVG da arquitetura precisou de tratamento à parte: meu corretor pulava `.svg` por
extensão, e os travessões ali são texto renderizado, visível dentro da figura.

### Reestruturação das branches

A topologia anterior era um campo minado. A 05 e a 06 carregavam os arquivos da 04, a
06b carregava tudo, e a 04 e a 07 disputavam o `.gitignore`. Sobreposições medidas:

```
07 x 04    .gitignore
07 x 06b   118 arquivos
02 x 27    README.md, README.pt-BR.md, architecture.svg
04 x 05    11 arquivos
```

Refeitas como cadeia estritamente linear, cada uma sobre a anterior. Isso elimina
conflito por construção, que era a exigência.

```
master
 -> pr1-repo-hygiene      13 commits   118 arquivos
 -> pr2-readme            14 commits     4 arquivos
 -> pr3-containers        13 commits    11 arquivos
 -> pr4-windows-guide      4 commits     2 arquivos
 -> pr5-gallery            5 commits    12 arquivos
```

Dois commits caíram no rebase por terem virado redundantes: o que ancorava os padrões
do `.gitignore`, já feito pela pr1, e um `style` idêntico ao de baixo.

### Verificação antes de enviar

- Ponta da cadeia idêntica ao tronco em todo arquivo que sobe: nenhuma divergência
- Nenhum caractere proibido em nenhuma das cinco branches
- Cada branch aplica limpo sobre a anterior (`git merge-tree`)
- `upstream/master` continua em `f5eb580`, sem drift
- Sintaxe: entrypoint, render.py, compose e SVG
- `.gitignore` mantém semântica em onze caminhos-sonda; nada versionado virou ignorado

### Pull requests

| PR | Título | Estado |
|---|---|---|
| [#3](https://github.com/antoniocastelofilho/HigFlow/pull/3) | Remove build artefacts from tracking | MERGEABLE |
| [#4](https://github.com/antoniocastelofilho/HigFlow/pull/4) | Rewrite README in English | MERGEABLE |
| [#5](https://github.com/antoniocastelofilho/HigFlow/pull/5) | Add reproducible container images | MERGEABLE |
| [#6](https://github.com/antoniocastelofilho/HigFlow/pull/6) | Add a Windows installation guide | MERGEABLE |
| [#7](https://github.com/antoniocastelofilho/HigFlow/pull/7) | Add results gallery, fix example2d_Oldroyd | MERGEABLE |

Limitação conhecida do modelo empilhado entre forks: o GitHub não aceita como base uma
branch que só existe no fork, então os cinco apontam para `master` e o diff de cada um
inclui os anteriores até que sejam merjados. Isso está dito no corpo de cada PR.

Duas correções pegas na revisão final: um parágrafo duplicado no guia de container,
criado pela resolução de conflito, e um comentário no Dockerfile que ainda descrevia o
PETSc configurado com `--with-mpi-dir` depois da troca para os wrappers.

---

## 2026-08-21 - E05 e E06: Windows e galeria

**Branches:** `juniormar/05-install-windows` (de E04), `juniormar/06-gallery` (de E04),
`juniormar/06b-readme-gallery` (do tronco) para todas merjadas em `juniormar/main`

### E05 - guia Windows

358 linhas, escritas a partir da instalação real, não de conselho genérico. Cobre as
quatro surpresas que custaram tempo: a mensagem sobre WSL1 ser irrelevante, o
`VirtualizationFirmwareEnabled: False` ser artefato de leitura, o `DefaultUid = 0`
quando o setup não é concluído, e o `df -h /` reportar o tamanho máximo do disco
virtual em vez do espaço livre real.

### E06 - galeria

8 figuras de 3 casos, todas reproduzíveis. Renderizador próprio (`tools/gallery/render.py`),
sem ParaView: lê o VTK ASCII e desenha as células como estão.

**Números verificados:**

| Caso | Verificação | Resultado |
|---|---|---|
| Poiseuille | erro contra `u = u_max(1-y²)` | L₂ rel. **6,89×10⁻⁴** |
| Poiseuille | gradiente de pressão | −3, batendo com −2μu_max/h² |
| Poiseuille | vazão em 7 estações | constante em 5 partes em 10⁶ |
| Contração | conservação de massa | 0,37% através de razão de área 4:1 |
| Oldroyd-B | N₁ = τxx − τyy | máximo nas paredes, nulo na linha de centro |

### O achado maior da etapa

**Os dez casos versionados declaram `flowphase: singlephase` e `flowtype: newtonian`**,
independentemente do nome e do que o driver chama. Onde o driver chama
`higflow_solver_step()` isso é consistente; onde não, é fatal - o `flowtype` decide
quais propriedades distribuídas o `hig-flow-kernel.c` aloca, então o passo viscoelástico
ou multifásico desreferencia arrays nunca criados.

| Caso | Como versionado | Trocando a linha |
|---|---|---|
| `example2d_Oldroyd` | SIGSEGV no passo 0 | roda completo - **corrigido aqui** |
| `example2d_VOF` | SIGSEGV no passo 0 | ainda falha - **reportado, não corrigido** |

Confirmei a causa empiricamente: com `flowtype: viscoelastic` o caso roda os 101 passos
e o VTK passa a conter `TENSORS τₚ`, o tensor de tensão polimérica que antes nunca era
alocado.

### Correção que fiz de um diagnóstico meu

Cheguei a concluir que o `flowtype: newtonian` tornava o caso Oldroyd newtoniano. Errado:
o **driver** escolhe a física, chamando `higflow_solver_step_viscoelastic()` diretamente.
O `flowtype` não escolhe o solver - ele decide a alocação. Verifiquei antes de escrever
qualquer legenda.

### Nota operacional

O `/tmp` do WSL não sobrevive entre invocações - apagou logs e resultados duas vezes.
Trabalho persistente vai em `/root/hf`.

### Pendências

- Figura multifásica: bloqueada pelo `example2d_VOF`
- Os outros 6 casos com a mesma declaração não foram testados
- `higflow:latest` fica na máquina; reconstruir leva 6min30s

---

## 2026-08-21 - E00 e E04: ambiente e containers

**Etapas:** E00, E04
**Branch:** `juniormar/04-containers` (a partir de `master`) para merjada em `juniormar/main`

### Ambiente (E00)

O WSL já estava completo - WSL 2.6.3.0, kernel 6.6.87.2-1, `WslService` rodando,
hypervisor ativo. Faltava só a distro, que o usuário instalou. Não houve reinício.

Dois pontos do ambiente que valem registro:

- **`DefaultUid = 0`** - o setup inicial do Ubuntu não criou usuário normal, então
  tudo roda como root. Isso quebrou o primeiro teste do container por permissão, e
  também impede o OpenMPI de rodar sem `--allow-run-as-root`. Vale criar um usuário.
- **Docker Engine instalado dentro do WSL**, não Docker Desktop. O `systemd` já
  estava habilitado (`/etc/wsl.conf` com `[boot] systemd=true`), que é o requisito.
  Sem licença, mais leve, e é o caminho documentado no guia.

Validações que economizaram um build longo: `hdf5.pc` **existe** no Ubuntu 22.04
(o Makefile chama `pkg-config ... hdf5`), headers do Zoltan em `/usr/include/trilinos`
como o Makefile do higtree assume, e boost é header-only aqui - `libboost-dev` basta,
não `libboost-all-dev`.

### Container (E04)

11 commits. Imagem multi-stage de 4 estágios, imagem de desenvolvimento, compose com
7 serviços, definição Apptainer, `.dockerignore`, entrypoint e o guia de 473 linhas.

**Medições:**

| | |
|---|---|
| Build a frio | **6 min 30 s** (16 núcleos) - PETSc é 5 min 3 s |
| Contexto enviado | **11,67 MB** de uma árvore de 96 MB |
| Imagem | 299 MB de conteúdo, 1,34 GB em disco |
| Rebuild após mudar fonte | 36 s |

**Casos executados, saída no host:**

| Caso | Ranks | Saída |
|---|---|---|
| `example2d_Newt` | 1 | 101 VTK, 166 MB |
| `example2d_Newt` | 2 | 202 VTK, 166 MB |
| `example3d_lid_driven` | 1 | 880 VTK, 509 MB (interrompido após verificar) |

### Quatro defeitos novos, todos achados porque a imagem não construía ou não rodava

1. **`higflow/Makefile:83` não põe a dimensão no nome do objeto.** Usa
   `hig-flow-%.o`; o `higtree/Makefile:90` usa `%-$(DIM)d.o` corretamente. Logo
   `make DIM=2 && make DIM=3` linka objetos 2D dentro de `libhigflow3d.a`. E
   `make clean` não resolve: apaga `$(HIGFLOW_LIBPATH)/*.a`, levando a outra junto.

2. **Os exemplos não podem ser movidos.** Todos os 11 incluem `../src/hig-flow-*.h`
   e linkam `../src/hig-flow-*.o` diretamente, não `libhigflow<dim>d.a`. Copiar um
   caso para fora de `higflow/` quebra a compilação. Somado ao defeito 1: os 9
   exemplos 2D e os 2 exemplos 3D **não podem estar compiláveis ao mesmo tempo**.

3. **Os exemplos compilam com `gcc`, não `mpicc`.** Os 11 usam `CC = gcc`; as duas
   bibliotecas usam `CC = mpicc`. A única fonte de caminho de MPI do exemplo é o
   `$(PETSC_CC_INCLUDES)`.

   **Isso explica o `--download-openmpi` do instalador**, que eu tinha catalogado
   como parte do problema dos três MPIs. É, em parte, contorno *deste* defeito: o
   MPI próprio do PETSc põe `mpi.h` dentro do prefixo do PETSc, então o
   `PETSC_CC_INCLUDES` cobre por acidente de layout. Apontar para um MPI de sistema
   desfaz o acidente. E no Debian/Ubuntu não existe raiz única de MPI:
   `--with-mpi-dir=/usr` grava include sem `mpi.h`, e o prefixo do OpenMPI faz o
   configure falhar com `Fortran error! mpi_init() could not be located!`.

4. **`.gitignore` bloqueando `docs/install/`** - o defeito previsto na E07,
   atrapalhando de verdade: `git add docs/install/containers.md` não fez nada.

### Correção de processo

Três arquivos (`cases/.gitkeep`, `.gitignore`, `.dockerignore`) entraram no commit
do Dockerfile porque já estavam no índice. Mensagem não batia com conteúdo, contra
a própria regra de "um commit, uma ideia". Como nada havia sido enviado, refiz com
`git reset --mixed` e separei em 4 commits limpos.

### Próxima sessão

**E05** (guia Windows) é a continuação natural - é o registro escrito da experiência
E00 + E04. Ou **E06** (galeria), que agora está desbloqueada: a imagem roda os casos.

---

## 2026-08-21 - E07: higiene do repositório

**Etapa:** E07
**Branch:** `juniormar/07-repo-hygiene` (a partir de `master`) para merjada em `juniormar/main`

### Resultado

**Árvore versionada: 1.507 arquivos / 95,3 MB para 1.406 arquivos / 18,6 MB - redução de 76,7 MB (80%).**

### Commits (12)

| SHA | Assunto |
|---|---|
| `65aea28` | 14 binários ELF destrackeados (37,2 MB) |
| `1d8bf6b` | Tarball do PETSc para `fetch-petsc.sh` com verificação SHA-256 (37,8 MB) |
| `eb23e6a` | `src.zip`, `.swp`, `contr.flowtype` |
| `33466e2` | 7 `*_old.c` + `Attic/` (713 KB) |
| `2673acb` | `include/` gerados (62 arquivos) |
| `e6a654d` | `.gitignore` reescrito |
| `d61de31` | Saídas de simulação versionadas (20 arquivos, 1,7 MB) |
| `4fc6ceb` | `.gitattributes` |
| `09da0ee` | `.mailmap` |
| `a425bb8` | `.editorconfig` |
| `efd69fa` | `.clang-format` (raiz + higtree) |
| `c067502` | `bibliotecas/README.md` |

### Achados novos

1. **`**build**` e `**install**` no `.gitignore` casam com qualquer caminho contendo a
   substring.** Fora das três formas especiais (`**/`, `/**`, `/**/`), asteriscos
   consecutivos colapsam para um só. Consequência verificada:

   ```
   git check-ignore --no-index -v higtree/src/build-fringe.cpp
   .gitignore:73:**build**  higtree/src/build-fringe.cpp
   ```

   `build-fringe.cpp` é compilado pelos dois sistemas de build. Sobrevive apenas por ter
   sido commitado antes da regra existir. A mesma regra engole
   `install_higflow_arch.sh` (da branch `PC_ImproveDocumentation`) e qualquer coisa sob
   `docs/install/` - exatamente onde E03 e E05 vão escrever.

2. **`contr.flowtype` é saída de grep.** 2.721 linhas do tipo
   `hig-flow-bc.c:126:    if (ns->contr...`, resultado de uma busca redirecionada para
   um arquivo com o nome do termo buscado. Cópia idêntica em `src_hugo/old/`.

3. **Havia 14 binários ELF versionados, não 11.** Faltavam três cópias de
   `generate_amr` em `example2d_VOF*/amrs/` - com `generate_amr.c` ao lado.

4. **Saída 3D dentro de exemplo 2D.** `example2d_Newt/output/example-3d.save.pres` e
   `.vel` não podem ter sido produzidos pelo caso que acompanham.

5. **Os dois lados do projeto têm convenções de indentação opostas.** Medido:
   `higflow/src` 98% espaços, `higtree/src` 89% tabs. Por isso dois `.clang-format`, com
   o de `higtree/` herdando o da raiz e sobrescrevendo só a indentação.

### O `.mailmap` funcionou

De 4º lugar com 49 commits para **1º lugar com 128 commits**:

```
antes                                     depois
  69  Pedro Coimbra                        128  Juniormar Organista
  49  juniormar                             69  Pedro Coimbra
  40  Daniel Garcia                         40  Daniel Garcia
  38  Juniormar Organista <usp.br>          28  Kainã
  27  Juniormar Organista <alumni>           4  Johnatas
  23  kainaas                                3  Antonio Castelo Filho
  13  juniormarorganista
   5  Kainã
   1  Juniormar Organista <gmail>
```

### Verificação

- 26 módulos do `higflow/Makefile`, 29 do `higtree/Makefile`, todas as fontes dos dois
  `CMakeLists.txt` e todos os drivers de exemplo continuam versionados
- Nenhum arquivo versionado é casado pelo novo `.gitignore`
- `bash -n` passa nos três scripts tocados
- Cópias removidas da working tree confirmadas byte-idênticas aos blobs antes de apagar

### Decisões deliberadamente não tomadas

- **`src_hugo/` intocado.** 2,4 MB de cópia paralela da árvore, com diretório cujo nome
  contém espaços. É trabalho de outro pesquisador - vira pergunta no PR, não remoção.
- **Histórico não reescrito.** Restam ~160 MB no histórico (`.avi` de 45 MB, VTKs de 38
  e 20 MB, árvore `atf-0.15/`). `git filter-repo` quebraria todos os clones - fica
  registrado para os donos decidirem.
- **`libfyaml-master.zip` mantido.** Snapshot de `master`, não release taggeada; trocar
  por download exige escolher uma versão, o que muda contra o que o projeto compila.

### Rascunhos de PR prontos

`docs/contrib-plan/pr-drafts/E02-readme.md` e `E07-repo-hygiene.md`.

### Próxima sessão

**E00** (WSL2) continua sendo o gargalo: bloqueia E03, E04, E05, E06, E11, E12 e todo o
bloco C++.

---

## 2026-08-21 - E02: README em inglês

**Etapa:** E02
**Branch:** `juniormar/02-readme` (a partir de `master` em `f5eb580`) para merjada em `juniormar/main`

### O que foi feito

- `README.md` reescrito em inglês, 503 linhas, 12 seções
- `README.pt-BR.md` preservado, com o bloco corrompido reparado
- `docs/images/architecture.svg` - diagrama autoral das três camadas de software e dos
  seis estágios do passo de tempo
- Figuras existentes (`all_mesh.png`, `ghosts.png`, `wccmjet.png`) referenciadas **no
  lugar**, sem duplicar binário na árvore

### Commits (11)

| SHA | Mensagem |
|---|---|
| `da50eac` | `docs: preserve the Portuguese README as README.pt-BR.md` |
| `e730508` | `docs: repair corrupted sentence in the Portuguese README` |
| `4a25ae7` | `docs(readme): add English overview and table of contents` |
| `b2d1b33` | `docs(readme): document the constitutive model library` |
| `5776b94` | `docs(readme): document the numerical methods` |
| `69d9244` | `docs(readme): add architecture diagram` |
| `21e4304` | `docs(readme): add architecture section with grid figures` |
| `859e5fb` | `docs(readme): add gallery placeholder and installation instructions` |
| `e5e988e` | `docs(readme): document how to run a case and what a case contains` |
| `1523294` | `docs(readme): add layout, documentation index, citing, contributing and license` |
| `a8a7c0d` | `docs(readme): use contributor names exactly as recorded in git history` |

### Achados novos durante a redação

1. **O esquema convectivo de segunda ordem degrada para primeira ordem na fronteira.**
   `hig-flow-discret.c:131-138` rebaixa `SECOND_ORDER` para `FIRST_ORDER` quando o
   estêncil cruza a fronteira. É escolha deliberada de robustez, mas afeta diretamente
   a ordem de convergência observada - insumo importante para E11.

2. **`ORDER4` é parseado mas nunca usado.** Aparece em `hig-flow-io.c:6229` e `:8015`
   (leitura e escrita do YAML) e na declaração do enum. Nenhuma rotina de
   discretização o consulta. Selecionar `forth_order` não produz esquema de quarta
   ordem - produz silenciosamente segunda ordem.

3. **Nenhuma referência bibliográfica no código.** `higtree/doc/biblio.bib` tem 398
   entradas, nenhuma de reologia. A tabela de modelos do README traz as referências
   canônicas onde a formulação é inequívoca, e declara a lacuna onde não é, em vez de
   atribuir por suposição.

### Correção feita durante a etapa

Na primeira redação da seção Authors atribuí o sobrenome "Sousa" a um contribuidor
cujo histórico registra apenas "Kainã". Corrigido em `a8a7c0d`; a lista agora reproduz
exatamente o que o histórico registra.

### Verificação

- Todos os caminhos de imagem, links de arquivo e âncoras internas resolvem
- Tabelas markdown bem formadas
- SVG validado como XML e sem transbordo de geometria
- Diff do PR contém apenas `README.md`, `README.pt-BR.md` e o SVG - o diretório
  `docs/contrib-plan/` não vazou, como o modelo de branches previa

### Pendências herdadas para etapas seguintes

- Galeria com resultados reproduzíveis para **E06**
- Instalador corrigido, hoje apenas sinalizado no README para **E03**
- Seção de container para **E04**
- `CITATION.cff` e `CONTRIBUTING.md`, hoje declarados ausentes no README para **E19**
- Badges de CI reais para **E12**

### Próxima sessão

**E00** (WSL2) para desbloquear metade do roadmap, ou **E07** (higiene) que é PR
pequeno e independente.

---

## 2026-08-20 - Análise e planejamento

**Etapa:** E01
**Branch:** `juniormar/main` (criada a partir de `master` em `f5eb580`)

### O que foi feito

- Análise completa do repositório: 1.507 arquivos, ~176 mil linhas de C/C++, 99 MB de
  working tree e 180 MB de `.git`
- Levantamento das branches não merjadas do upstream, com auditoria do material de
  container existente
- Catalogação de defeitos com arquivo e linha: 11 de build, 15 de instalação e
  ambiente, 11 de correção de código
- Medição do atrito de uso: 422 linhas de YAML e 280 a 2.609 linhas de driver C por caso
- Avaliação de portabilidade Windows e de viabilidade da migração C++
- Dossiê, roadmap de 27 etapas, modelo de git e ferramenta de status

### Commits

| SHA | Mensagem |
|---|---|
| `b442e4b` | `docs(contrib): add technical analysis dossier for contribution planning` |
| `cb7230b` | `docs(contrib): add staged contribution roadmap` |
| `a1ae76d` | `docs(contrib): define branch, commit and pull request workflow` |
| `c37fc49` | `feat(contrib): add contribution status tool` |
| `f9bd2d1` | `chore(contrib): keep LF endings on shell scripts` |

### Achados que mudaram o plano

1. **O upstream tem trabalho paralelo recente e não merjado.** `Kaina` (2026-08-05,
   72 commits) tem Doxygen, 9 tutoriais e `nn-weights.cpp` com inferência MLP;
   `PC_ImproveDocumentation` já tem Dockerfiles e scripts de instalação;
   `PC_Daniel_Mesh` corrigiu conservação de massa em VOF 3D. Container, manual de
   instalação e ML em C++ já têm trabalho iniciado por outros. Decisão: construir em
   cima, dando crédito.

2. **A máquina local não compila nada do projeto.** Docker ausente, WSL2 sem
   distribuição, `make` ausente. E00 passou a ser pré-requisito de metade do roadmap.

3. **O repositório não tem licença.** Nem `CONTRIBUTING`, `CITATION.cff`, CI ou
   `.clang-format`. Decisão: propor a licença no PR, deixando a escolha aos donos.

4. **`-Ofast` contradiz o próprio código.** A flag implica `-ffinite-math-only`, que
   autoriza o compilador a assumir que `Inf` e `NaN` nunca ocorrem - mas
   `hig-flow-step-electroosmotic.c:923,932` usa `INFINITY` como sentinela em laço de
   convergência, e a detecção de divergência depende de comparar `NaN`. É problema
   metodológico, não apenas de engenharia.

5. **Os dois sistemas de build compilam conjuntos diferentes de arquivos.** Nenhum
   compila o projeto inteiro. Provável raiz de boa parte da dificuldade relatada de
   fazer o projeto funcionar.

6. **101 commits do autor estão fragmentados em 4 identidades git.** Um `.mailmap` de
   5 linhas é a correção de maior retorno por esforço de todo o levantamento.

### Pendências desta etapa

- `.mailmap` - escrito no PR de E07, onde faz sentido junto com `.editorconfig` e
  `.gitattributes`

### Verificação de estado ao fim da sessão

```
branch                juniormar/main
commits ahead master  5
pushed to origin      não (trabalho mantido local, conforme decidido)
master vs upstream    em sincronia
working tree          99 MB
.git                  180 MB
```

### Próxima sessão

Sugestão registrada no roadmap: **E01 concluída para E02 (README v1)** como entrega
visível imediata, com **E00 (WSL2)** em paralelo, já que a compilação do PETSc roda
sozinha por 30 a 90 minutos.

---
