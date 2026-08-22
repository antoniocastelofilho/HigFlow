# HigFlow - Roadmap de Contribuição

> Documento de trabalho interno da contribuição de **Juniormar Organista**.
> Português por ser material de decisão pessoal. **Não integra Pull Requests** -
> tudo que sobe para o upstream é em inglês.
>
> Análise que fundamenta este plano: [`00-ANALISE.md`](00-ANALISE.md).
> Modelo de branches e PRs: [`02-GIT-WORKFLOW.md`](02-GIT-WORKFLOW.md).

---

## Decisões que orientam o plano

| Questão | Decisão |
|---|---|
| Trabalho paralelo no upstream | **Construir em cima, dando crédito.** Trazer o material das branches `PC_ImproveDocumentation` e `Kaina`, corrigir os defeitos reais, creditar os autores |
| Ambiente local | **WSL2 + Ubuntu 22.04 agora**, Docker Desktop numa etapa posterior para validar a imagem |
| Licença | **Propor no PR**, deixando a decisão final aos donos |
| Branches | **Tronco `juniormar/main` + uma sub-branch por etapa**, cada uma um PR independente |
| Simplificação de uso | **Runner genérico com expressões em YAML** - sem recompilar para mudar fronteira |
| Limpeza do repositório | **Remover do HEAD, preservar histórico** - sem reescrita de SHAs |
| Idioma | **Tudo em inglês**, incluindo tradução do material existente em português |
| C++ | **Planejar os 3 níveis agora**, decidir o nível de execução depois |
| CI | **Completo, com testes de regressão numérica** |
| Vitrine do README | Poiseuille 2D, contração 4:1 viscoelástica, VOF multifásico, malha AMR |
| ML | Documento de arquitetura **+ prova de conceito mínima** em C++ puro |

**Regra transversal:** nenhum artefato produzido menciona ferramentas de IA. Commits,
documentação, comentários de código e descrições de PR são redigidos como trabalho
autoral direto.

---

## Mapa das etapas

```
BLOCO 0 - FUNDAÇÃO (local, não vira PR)
  E00  Ambiente de desenvolvimento
  E01  Infraestrutura da contribuição

BLOCO A - APRESENTAÇÃO
  E02  README em inglês (v1)          ──┐
  E03  Manual de instalação - Linux     │
  E04  Containers (Docker + Apptainer)  ├─para E06
  E05  Manual de instalação - Windows   │
  E06  Galeria de resultados + README v2 ┘

BLOCO B - CÓDIGO E BUILD
  E07  Higiene do repositório
  E08  Build unificado
  E09  Correções de defeitos pontuais
  E10  Flags numéricas e reprodutibilidade
  E11  Suíte de verificação numérica     ──┐
  E12  Integração contínua                 │
  E13  Runner: scaffolding                 │
  E14  Runner: driver unificado            ├─para pré-requisito do BLOCO D
  E15  Runner: avaliador de expressões     │
  E16  Gerador de malha .amr               ┘

BLOCO C - DOCUMENTAÇÃO E DIVULGAÇÃO
  E27  Galeria de pessoas (opt-in)
  E17  Doxygen e documentação de API
  E18  Tradução completa PT para EN
  E19  Governança do projeto
  E20  Divulgação

BLOCO D - MIGRAÇÃO C++   (requer E11 concluída)
  E21  Plano de migração
  E22  Nível 1 - compilar como C++
  E23  Nível 2 - RAII
  E24  Nível 3 - templates e tipos fortes

BLOCO E - MACHINE LEARNING (longo prazo)
  E25  Arquitetura de ML no HigFlow
  E26  Prova de conceito - inferência MLP em C++ puro
```

### Grafo de dependências reais

| Etapa | Depende de | Motivo |
|---|---|---|
| E03, E04, E05 | E00 | Precisa de Linux para testar o que é escrito |
| E06 | E04 | Gerar imagens exige um ambiente que compile e execute |
| E11 | E08 | Testes precisam de um build confiável e único |
| E12 | E11, E04 | CI executa a suíte dentro do container |
| E14 | E13 | Driver unificado consome os templates do scaffolding |
| E15 | E14 | Expressões substituem as funções `get_*` do driver |
| E22-E24 | E11 | Só se refatora numérica com rede de proteção |
| E27 | E02 | Edita o README reescrito |
| E26 | E25 | - |

**Todas as demais etapas são independentes** e podem ser executadas em qualquer ordem.

---

# BLOCO 0 - FUNDAÇÃO

Trabalho local de preparação. Não gera Pull Request.

---

## E00 - Ambiente de desenvolvimento

**Objetivo:** ter uma máquina capaz de compilar, executar e depurar o HigFlow.

**Depende de:** nada. **Bloqueia:** E03, E04, E05, E06, E11, E12 e todo o Bloco D.

### Situação atual
Docker ausente, WSL2 sem distribuição, `make` ausente. Nada do projeto compila hoje.

### Tarefas

1. Instalar WSL2 com Ubuntu 22.04
   ```bash
   wsl --install -d Ubuntu-22.04
   ```
2. Configurar usuário, `apt update && apt upgrade`
3. Instalar o conjunto mínimo de build e validar cada dependência isoladamente,
   registrando o que funciona e o que falha - esse registro é insumo direto de E03
4. Decidir e documentar a localização do repositório: manter em `/mnt/c/dev/HigFlow`
   (acesso do Windows, I/O lento) ou clonar em `~/HigFlow` no sistema de arquivos do WSL
   (rápido, mas duas cópias). **Recomendação:** clonar dentro do WSL para trabalho de
   compilação, com o repositório Windows como área de escrita de documentos
5. Instalar ParaView ou VisIt no Windows para pós-processar os VTKs gerados no WSL
6. Docker Desktop fica para o início de E04

### Entregáveis
- Ambiente funcional
- `docs/contrib-plan/notas/E00-ambiente.md` - registro do que foi instalado, versões,
  o que falhou e como foi contornado

### Critério de aceitação
`cd higtree && make DIM=2` produz `libhig2d.a` sem erro. Se falhar, o log do erro é o
primeiro insumo de E08.

### Esforço
Baixo, mas com espera longa: a compilação do PETSc leva de 30 a 90 minutos.

---

## E01 - Infraestrutura da contribuição

**Objetivo:** montar o esqueleto de trabalho e a ferramenta de verificação de estado.

**Depende de:** nada.

### Tarefas

1. Branch `juniormar/main` criada *(feito)*
2. Dossiê de análise *(feito)*
3. Este roadmap *(feito)*
4. `02-GIT-WORKFLOW.md` - modelo de branches, convenção de commits, fluxo de PR
5. `tools/contrib/status.sh` - ferramenta de verificação de estado. Responde, num único
   comando:
   - em que branch estamos e o que há de não commitado
   - quais etapas têm branch criada, quais têm commits, quais foram merjadas no tronco
   - o que existe localmente e **não** foi enviado ao `origin`
   - o que o `origin` tem e o `upstream` não (candidatos a PR)
   - o que o `upstream` ganhou desde o último `fetch` (risco de conflito)
   - estado de cada PR aberto, se o `gh` estiver disponível
6. `docs/contrib-plan/LOG.md` - diário de execução, uma entrada por sessão de trabalho
7. `.mailmap` consolidando as 4 identidades git

### Entregáveis
- `docs/contrib-plan/02-GIT-WORKFLOW.md`
- `docs/contrib-plan/LOG.md`
- `tools/contrib/status.sh`
- `.mailmap` *(vai para o PR de E07)*

### Critério de aceitação
`bash tools/contrib/status.sh` roda em Git Bash no Windows e no WSL, e produz um
panorama legível sem depender de rede.

---

# BLOCO A - APRESENTAÇÃO

O bloco que o usuário pediu para vir primeiro. Objetivo: tornar o projeto atraente,
instalável e compreensível para quem chega de fora.

---

## E02 - README em inglês (versão 1)

**Branch:** `juniormar/02-readme` · **PR:** "Rewrite README in English with project overview and quick start"

**Objetivo:** substituir o README atual - que está em português, contém um bloco de
texto corrompido e não mostra o que o código faz - por um documento que apresente o
projeto.

**Depende de:** nada. Usa figuras já existentes em `higtree/doc/figuras`; as imagens de
resultados reais entram em E06.

### Problemas do README atual
- Português, sem versão em inglês
- Bloco corrompido: `"Após re./configure --prefix=... ; sleep 5alizar um dos passos"`
- Instruções contraditórias com `varsrc` (defeito I11 do dossiê)
- Nenhuma imagem, nenhuma descrição da física simulada
- Nenhum badge, nenhuma citação, nenhuma licença

### Estrutura proposta

```
┌ Logo / banner + badges (build, licença, DOI, versão)
├ Uma frase: o que é o HigFlow
├ Figura de destaque (em E06 vira resultado real; agora, malha AMR)
├ Features - tabela de modelos constitutivos e métodos numéricos
├ Quick start - três caminhos:
│    · Docker (uma linha)
│    · Linux nativo
│    · Windows via WSL2
├ Running your first case - Poiseuille 2D do início ao VTK
├ Gallery - grade de resultados (placeholders até E06)
├ Documentation - links para manuais, API, tutoriais
├ Project structure - o que é higtree, o que é higflow
├ Citing HigFlow
├ Contributing
├ Authors and acknowledgements
└ License
```

### Tarefas
1. Redigir o README em inglês
2. Corrigir toda a informação factualmente errada, cruzando com o que funcionou em E00
3. Tabela de modelos constitutivos com a referência bibliográfica de cada um -
   é o que distingue o HigFlow de um CFD genérico
4. Badges apontando para os workflows que E12 vai criar (placeholders até lá)
5. Preservar o conteúdo em português como `README.pt-BR.md` nesta etapa; a tradução
   integral dos demais documentos é E18

### Entregáveis
`README.md` (inglês), `README.pt-BR.md`, `docs/images/` com as figuras selecionadas

### Critério de aceitação
Alguém que nunca viu o projeto entende, em menos de um minuto de leitura, o que ele
faz, se serve para o problema dele, e qual comando executar primeiro.

### Commits previstos
5 a 8, um por seção maior.

---

## E03 - Manual de instalação - Linux nativo

**Branch:** `juniormar/03-install-linux` · **PR:** "Rewrite Linux installation: idempotent script, correct PETSc flags, verification step"

**Objetivo:** substituir `install_higflow_ubuntu22` por um instalador que funcione.

**Depende de:** E00 (é preciso ter executado a instalação de verdade para escrever sobre ela).

**Trabalho prévio a incorporar:** `upstream/PC_ImproveDocumentation` já converteu o
script para `install_higflow_ubuntu22.sh` e adicionou `install_higflow_arch.sh`.
Partir desses arquivos e creditar Pedro Coimbra.

### Defeitos a corrigir (I1-I12 do dossiê)

| Defeito | Correção |
|---|---|
| Sem shebang, sem `set -euo pipefail` | Adicionar; falhar cedo e com mensagem |
| `--with-debubbing=yes` | `--with-debugging=0` para build de produção |
| `--PETSC_ARCH=x86_64` | `PETSC_ARCH=x86_64` (sintaxe correta do PETSc) |
| Três MPIs concorrentes | Escolher **um**. Recomendação: OpenMPI do sistema, e `--with-mpi-dir=/usr` no PETSc em vez de `--download-openmpi` |
| `PKG_CONFIG_PATH` recebendo arquivo | Apontar para o diretório e **acrescentar**, não sobrescrever |
| `export` que não persiste | Gerar um `env.sh` versionável e instruir o `source` |
| `sudo ./configure` | Construir como usuário, `sudo` só no `make install` |
| `chmod 777` | `chmod 755` |
| `pip3` global | `venv`, ou remover - as dependências de Sphinx não são necessárias para compilar |
| Symlink `libHYPRE.so` | Investigar a causa; resolver via `pkg-config` ou `-L` correto |
| `varsrc` contraditório | Gerar `varsrc` **a partir** do que o instalador de fato instalou |
| `$(pwd)` em `varsrc` | Resolver o caminho do próprio script (`${BASH_SOURCE[0]}`) |

### Tarefas adicionais
1. **Idempotência:** rodar duas vezes não pode quebrar nada, e a segunda execução deve
   pular o que já está pronto
2. **Detecção de distribuição:** Ubuntu 22.04, Ubuntu 24.04, Debian, Arch, Fedora -
   ao menos detectar e avisar quando não suportada, em vez de falhar de forma obscura
3. **Etapa de verificação:** ao final, compilar e rodar o caso Newtoniano, e reportar
   sucesso ou falha explicitamente. Hoje o script termina com "THE END" mesmo se tudo
   falhou
4. **Modo `--dry-run`** listando o que seria instalado
5. Documento `docs/install/linux.md` explicando cada dependência, por que é necessária
   e o que fazer quando falha

### Entregáveis
`install/linux.sh`, `install/arch.sh`, `docs/install/linux.md`, `etc/varsrc.in`

### Critério de aceitação
Numa Ubuntu 22.04 limpa (container descartável), do `git clone` ao VTK do primeiro
caso, sem intervenção manual.

### Risco
`--download-openmpi` versus MPI do sistema é uma decisão com efeito colateral: quem já
tem ambiente montado pode ter o seu quebrado. Tratar como opção, com o padrão sendo o
MPI do sistema.

---

## E04 - Containers

**Branch:** `juniormar/04-containers` · **PR:** "Add reproducible Docker and Apptainer images with multi-stage build"

**Objetivo:** tornar `docker run` o caminho recomendado - nenhuma dependência para
instalar, funciona igual em Windows, Linux e macOS.

**Depende de:** E00. Docker Desktop é instalado no início desta etapa.

**Trabalho prévio a incorporar:** `upstream/PC_ImproveDocumentation` tem
`conteiner/Dockerfile` e `conteiner/Dockerfile.petsc`. Partir deles, corrigir os 8
defeitos catalogados em §3.1 do dossiê, creditar Pedro Coimbra.

### Arquitetura proposta

```
containers/
├── Dockerfile                 multi-stage real, 3 estágios
│     stage 1 "deps"           dependências de sistema (camada estável, cacheável)
│     stage 2 "petsc"          PETSc compilado (camada cara, cacheável)
│     stage 3 "higflow"        o projeto (camada volátil)
├── Dockerfile.dev             + gdb, valgrind, clangd, ParaView headless
├── docker-compose.yml         serviços: build, run, shell, test
├── apptainer.def              conversão do higflow_image.def existente
├── .dockerignore              exclui .git, bibliotecas/, binários, VTKs
└── entrypoint.sh              carrega o ambiente e despacha o comando
```

### Correções obrigatórias sobre o material existente

| Defeito herdado | Correção |
|---|---|
| `echo '. $HOME/.varsrc"' >> .bashrc` | Aspa balanceada. **Este bug quebra todo shell do container** |
| `--with-debubbing=yes` | `--with-debugging=0` |
| Cabeçalho "OpenFOAM-9 + Python 3.10" | Descrever a imagem real |
| `FROM ubuntu_petsc3.14:v01` | Multi-stage real, um único `docker build` |
| Ausência de `.dockerignore` | Criar - evita 78 MB de artefatos entrarem na imagem |
| `chmod 777` | `chmod 755`, usuário não-root |
| Sem versões fixadas | Fixar `FROM ubuntu:22.04@sha256:...` e as versões críticas |

### Melhorias além da correção
1. **Cache de camadas ordenado** por volatilidade: PETSc antes do código-fonte. Alterar
   uma linha de C não deve recompilar o PETSc
2. **Volume para dados:** `docker run -v $PWD/cases:/work` - casos e saídas ficam no
   host, o container é descartável
3. **Suporte MPI:** documentar `--cap-add=SYS_PTRACE` e `--shm-size`, necessários para
   OpenMPI dentro de container
4. **Imagem multi-arquitetura** (`amd64` e `arm64`) - Apple Silicon e clusters ARM
5. **Publicação** no GitHub Container Registry via workflow de E12
6. **Apptainer/Singularity** para clusters HPC (onde Docker não roda por política de
   segurança): converter `stacks/singularity/higflow_image.def`, corrigindo o caminho
   relativo quebrado em `%post` e o `cp` para `/pacotes` inexistente

### Guia detalhado - `docs/install/containers.md`
O usuário pediu explicação detalhada de containers. O documento cobre:
- O que é um container e por que resolve o problema do HigFlow especificamente
- Diferença entre imagem, container e volume, na prática deste projeto
- Docker versus Apptainer: quando usar cada um (Apptainer em cluster, Docker no desktop)
- Instalação do Docker em Windows, Linux e macOS
- Construir a imagem, executar um caso, extrair resultados
- Rodar em paralelo com MPI dentro do container
- Montar o código-fonte para desenvolver dentro do container sem reconstruir a imagem
- Solução de problemas: permissões, `shm`, memória, MPI que não inicia
- Como o mesmo `.def` roda num cluster com Slurm

### Entregáveis
Diretório `containers/`, `docs/install/containers.md`

### Critério de aceitação
```bash
docker build -t higflow containers/
docker run --rm -v "$PWD/cases:/work" higflow run /work/poiseuille.yaml
```
produz VTKs em `./cases`, numa máquina sem nenhuma dependência do HigFlow instalada.

### Esforço
Alto. É a etapa de maior valor por unidade de esforço do bloco.

---

## E05 - Manual de instalação - Windows

**Branch:** `juniormar/05-install-windows` · **PR:** "Add Windows installation guide (WSL2 and Docker Desktop)"

**Objetivo:** documentar o uso em Windows - hoje inexistente.

**Depende de:** E00 (WSL2) e E04 (Docker). É o registro escrito da experiência real
dessas duas etapas.

### Posição técnica a documentar honestamente

O dossiê (§6) estabelece que **porte nativo para Windows não é viável**: `libnuma` é
exclusivo de Linux e está marcado `REQUIRED` no CMake; OpenMPI não tem suporte Windows;
libfyaml não tem porte oficial. O manual afirma isso claramente e apresenta os dois
caminhos que funcionam, em vez de prometer algo que quebraria.

### Conteúdo

1. **Caminho A - Docker Desktop** (recomendado para usar)
   - Instalação, backend WSL2, requisitos de hardware
   - Executar um caso, onde os arquivos aparecem no Windows
   - ParaView no Windows lendo VTKs gerados no container
2. **Caminho B - WSL2 + Ubuntu** (recomendado para desenvolver)
   - `wsl --install -d Ubuntu-22.04`
   - Onde clonar: `/home/user/` versus `/mnt/c/` e o impacto real de I/O
   - Integração com VS Code via extensão WSL
   - Depuração com gdb dentro do WSL
   - Acesso aos arquivos do Linux pelo Explorer via `\\wsl$\`
3. **Caminho C - MSYS2** - documentado como **não suportado**, com a justificativa
   técnica, para que ninguém perca tempo tentando
4. Solução de problemas específicos de Windows: virtualização desabilitada na BIOS,
   consumo de memória do WSL (`.wslconfig`), antivírus interferindo na compilação,
   fim de linha CRLF quebrando scripts shell

### Melhorias de código que esta etapa justifica
Duas correções pequenas que ampliam a portabilidade e entram no PR de E09:
- `numactl` deixa de ser `REQUIRED` no CMake, tornando-se opcional junto com
  `solver-sor` - beneficia também macOS e clusters sem libnuma
- `ftime()` para `clock_gettime(CLOCK_MONOTONIC, ...)`

### Entregáveis
`docs/install/windows.md`, seção Windows no README, `.gitattributes` tratando CRLF

### Critério de aceitação
Um usuário Windows sem experiência em Linux chega ao primeiro VTK seguindo apenas o
documento.

---

## E06 - Galeria de resultados e README v2

**Branch:** `juniormar/06-gallery` · **PR:** "Add results gallery with reproducible cases and figures"

**Objetivo:** o pedido central - README com imagens de problemas resolvidos pelo próprio código.

**Depende de:** E04. Sem ambiente que compile e execute, não há imagens.

### Casos selecionados

| # | Caso | Por quê | Custo |
|---|---|---|---|
| 1 | **Poiseuille 2D** (`example2d_Newt`) | Solução analítica conhecida - a figura mostra numérico versus analítico sobrepostos. Serve de vitrine **e** de teste de validação em E11 | Minutos |
| 2 | **Contração 4:1 viscoelástica** (Oldroyd-B / GPTT) | Benchmark clássico de reologia computacional. Vórtice de canto visível. É o que diferencia o HigFlow de um CFD genérico | Horas |
| 3 | **VOF multifásico** (bolha ou dam break) | Visualmente o mais forte. Interface deformável, campo de fração volumétrica. Candidato a GIF animado no topo do README | Horas |
| 4 | **Malha AMR e decomposição de domínio** | Mostra o diferencial arquitetural do higtree: refinamento hierárquico e partição MPI. Há figuras base em `higtree/doc/figuras` | Baixo |

### Tarefas

1. Executar cada caso no container, registrando parâmetros exatos e tempo de execução
2. Pós-processar em ParaView; salvar o **state file** `.pvsm` junto com a figura -
   é o que torna a figura reproduzível por terceiros
3. Script `tools/gallery/render.py` (pvpython) que regenera todas as figuras a partir
   dos VTKs, sem interação manual
4. Para o caso 1, gráfico de validação: perfil numérico contra
   `u(y) = u_max(1 − y²)`, com a norma do erro anotada na figura
5. GIF ou MP4 curto do caso 3 (limite: 5 MB; se maior, hospedar como release asset e
   linkar)
6. `docs/gallery.md` com, para cada caso: figura, física, parâmetros adimensionais,
   comando exato de reprodução e referência bibliográfica
7. Atualizar o README com a grade de figuras
8. Registrar em `.gitattributes` os binários da galeria e avaliar Git LFS se passar de
   ~10 MB no total

### Entregáveis
`docs/gallery.md`, `docs/images/gallery/`, `tools/gallery/render.py`, states `.pvsm`,
README v2

### Critério de aceitação
Toda figura do README é regenerável por um terceiro a partir do repositório, com o
comando documentado. **Nenhuma imagem sem procedência.**

### Risco
Os casos 2 e 3 podem não convergir com os parâmetros versionados - os defeitos de
build (§4.2 do dossiê) sugerem que nem todos os exemplos estão em estado funcional.
Se um caso falhar, é achado de bug, não fracasso da etapa: vira issue no upstream, o
que também é contribuição.

---

# BLOCO B - CÓDIGO E BUILD

---

## E07 - Higiene do repositório

**Branch:** `juniormar/07-repo-hygiene` · **PR:** "Remove build artifacts and vendored archives from tracking; add editor and attribute configuration"

**Objetivo:** repositório limpo e clonável em segundos.

**Depende de:** nada.

### Remoções (`git rm --cached`, preservando histórico)

| Item | Tamanho | Substituição |
|---|---|---|
| 11 binários `ns-example` / `ns-complex-3d` | ~40 MB | `.gitignore` |
| `bibliotecas/petsc-3.14.0.tar.gz` | 37 MB | Download com verificação de checksum no instalador |
| `bibliotecas/libfyaml-master.zip` | 451 KB | Clone do repositório oficial, tag fixada |
| `higflow/src.zip` | 779 KB | - |
| `higflow/src/.hig-flow-kernel.h.swp` | 16 KB | - |
| 7 arquivos `*_old.c` / `*.old.c` | 689 KB | Preservados no histórico git; `git log --follow` os recupera |
| `higflow/include/`, `higtree/include/` | - | São cópias geradas por `cp src/*.h`. Passam a ser geradas no build |
| `higflow/include/hig-flow-step-multifase.h` | - | Órfão, grafia antiga |
| VTKs e `.dat` de saída versionados | ~6 MB | `.gitignore` |

### Decisão sobre `src_hugo/`
2,4 MB de cópia paralela da árvore, com um diretório cujo nome contém espaços.
**Não remover unilateralmente** - é trabalho de outro pesquisador. Abrir issue
perguntando se pode ir para uma branch de arquivo, e tratar a resposta como decisão dos
donos. Registrar a pergunta no PR.

### Adições

| Arquivo | Função |
|---|---|
| `.gitignore` reescrito | Corrigir `*.txt`/`*.dat`/`*.vtk` globais, que mascaram arquivos legítimos. Regras específicas por diretório |
| `.gitattributes` | Normalização de fim de linha (o repositório já emite avisos CRLF), marcação de binários, `linguist-vendored` para `bibliotecas/` |
| `.mailmap` | Consolida as 4 identidades git em uma |
| `.editorconfig` | Indentação e codificação consistentes |
| `.clang-format` | Baseado no estilo predominante do código atual, não num estilo importado |

### O que **não** será feito
`git filter-repo` para remover os 45 MB de `.avi` e os VTKs de 38 MB do histórico.
Reescreveria todos os SHAs e quebraria os clones de todos os colaboradores - inaceitável
num PR para repositório com autores ativos. Fica registrado como recomendação a ser
executada pelos donos, se decidirem.

### Efeito medido
Working tree: 99 MB para ~21 MB. Clone raso (`--depth 1`): ~78 MB mais leve.

### Critério de aceitação
`git clone --depth 1` e o build continua funcionando: nenhum arquivo removido era
necessário para compilar.

---

## E08 - Build unificado

**Branch:** `juniormar/08-build-system` · **PR:** "Unify build: single canonical CMake configuration"

**Objetivo:** eliminar a divergência entre CMake e Makefile - hoje a raiz da maior
parte das dificuldades de compilação.

**Depende de:** nada (mas E00 ajuda a validar).

### Problema
Os dois sistemas compilam **conjuntos diferentes de arquivos** (tabela em §4.2 do
dossiê). Nenhum compila o projeto inteiro. Dois usuários com dois métodos obtêm
binários com funcionalidades diferentes.

### Decisão de arquitetura
**CMake torna-se canônico.** Os Makefiles permanecem por uma versão, com aviso de
descontinuação, para não quebrar o fluxo de quem já os usa. Motivos: CMake já é usado
no cluster com Lmod, gera `compile_commands.json` (clangd, análise estática), integra
com CTest para E11 e é o que o CI de E12 vai consumir.

### Correções (B1-B11 do dossiê)

| Defeito | Correção |
|---|---|
| Ausência de `project()` | `project(HigFlow VERSION x.y.z LANGUAGES C CXX)` |
| `if(${debug})` quebra sem `-Ddebug` | `option(HIGFLOW_DEBUG "..." OFF)` com padrão |
| `-DDIM=${dim}` vazio sem `-Ddim` | `set(HIGFLOW_DIM 2 CACHE STRING ...)` com validação |
| `MFSIM_DEPENDENCIES` | Renomear; BLAS e LAPACK passam a ser efetivamente linkados |
| `-march=native` só em link options | Mover para compile options (e ver E10 sobre o padrão) |
| `include(.../CMakeLists.txt)` | `add_subdirectory()` |
| `-ltrinilos_zoltan` | `-ltrilinos_zoltan` |
| `$(RANLIB)` versus `ANLIB` | Corrigir o nome da variável |
| `build-fringe` duplicado | Remover a duplicata |
| `-I/usr/include/trilinos` hardcoded | `find_package` / `pkg-config` |
| `numactl REQUIRED` | Opcional, acoplado a `solver-sor` |

### Melhorias estruturais
1. **Uma lista de fontes, uma só vez** - arquivo `sources.cmake` incluído por ambos os
   sistemas enquanto o Makefile existir. Impossibilita nova divergência
2. **Alvos exportados:** `HigFlow::higtree`, `HigFlow::higflow`, com
   `target_include_directories` público - projetos externos consomem com
   `find_package(HigFlow)`
3. `compile_commands.json` habilitado por padrão
4. **Presets** (`CMakePresets.json`): `debug`, `release`, `release-reproducible`,
   `cluster`
5. **Resumo de configuração** ao fim do `cmake`: dimensão, dependências encontradas e
   ausentes, flags efetivas. Hoje o usuário descobre o que faltou por erro de link
6. Resolver o destino de `solver.camila.c` e `solver-petsc.camila.c`: verificar se
   duplicam símbolos de `solver.c`/`solver-petsc.c` e, em caso afirmativo, torná-los
   alternativa explícita ou removê-los

### Critério de aceitação
`cmake -B build && cmake --build build` funciona **sem argumento algum** e produz
binários com o conjunto completo de módulos, idêntico ao que o Makefile produz.

### Risco
Alto risco de conflito com `upstream/Kaina` e `upstream/PC_Daniel_Mesh`, que também
tocam Makefiles. Fazer `git fetch upstream` e conferir imediatamente antes de abrir o PR.

---

## E09 - Correções de defeitos pontuais

**Branch:** `juniormar/09-defect-fixes` · **PR:** "Fix debug macro leakage, format-string bug, and obsolete time API"

**Objetivo:** corrigir os defeitos catalogados C1-C11. Cada um é pequeno, isolado e
individualmente revisável - o formato ideal de PR.

**Depende de:** nada.

### Lista

| # | Correção | Nota |
|---|---|---|
| C1 | Remover `#define DEBUG` de `hig-flow-kernel.h:18` | Passa a ser controlado pelo build. **Atenção:** pode revelar código que dependia silenciosamente do caminho de debug - verificar cada `DEBUG_*` afetado |
| C2 | `Debug-c.h`: mover as `static` do cabeçalho para uma unidade de tradução, expondo por `extern` | A pilha de debug passa a funcionar como projetada |
| C3 | `DEBUG_WARNING(x)` para `fprintf(debugfd, "%s", x)` | Fecha o format string bug |
| C4 | `DEBUG_ASSERT`: `*(int*)NULL = 0` para `abort()` | Remove o comportamento indefinido |
| C5 | `ftime()` para `clock_gettime(CLOCK_MONOTONIC, ...)` | Também remove `<sys/timeb.h>`; ganho de portabilidade |
| C7 | `real x[DIM]={INFINITY}` para laço de inicialização ou `{INFINITY, INFINITY, INFINITY}` | 6 declarações. Defeito dormente mas real |
| C11 | `include/` gerado deixa de ser versionado | Feito em E07; aqui, garantir que o build gere |
| - | `numactl` opcional | Justificado por E05 |
| - | `ns-exemple-3d.c` para `ns-example-3d.c` | Typo em nome de arquivo versionado |

### Fora do escopo desta etapa
- C8 (função de 1.580 linhas) para refatoração natural em E14
- C9 (308 `exit()`) para é trabalho de RAII, vai para E23
- C10 (306 `fopen()`) para auditoria própria; abrir issue separada

### Critério de aceitação
Suíte de E11 passa antes e depois, com resultados numericamente idênticos. Se algum
resultado mudar, é sinal de que o comportamento dependia de um defeito - investigar
antes de prosseguir.

---

## E10 - Flags numéricas e reprodutibilidade

**Branch:** `juniormar/10-numerics-flags` · **PR:** "Replace -Ofast with -O3 and make architecture flags configurable"

**Objetivo:** corrigir o problema metodológico do `-Ofast`. Análise completa em §4.5 do
dossiê.

**Depende de:** E08 (o build precisa estar unificado para a mudança valer em todo lugar).
Idealmente E11 antes, para medir o efeito.

### O problema em uma frase
`-Ofast` implica `-ffinite-math-only`, que autoriza o compilador a assumir que `Inf` e
`NaN` nunca ocorrem - mas o código **usa `INFINITY` como sentinela** em laços de
convergência (`hig-flow-step-electroosmotic.c:923,932`;
`hig-flow-step-multiphase-electroosmotic.c:561,569`), e a detecção de divergência
depende de comparar `NaN`.

### Mudanças

| Antes | Depois |
|---|---|
| `-Ofast` em todos os alvos | `-O3` como padrão |
| `-ffast-math` implícito | Somente sob `-DHIGFLOW_FAST_MATH=ON`, documentado, com aviso de que invalida a detecção de divergência |
| `-march=native` fixo | `HIGFLOW_ARCH` configurável; padrão `-mtune=generic`. `native` continua disponível e continua sendo o recomendado para produção - mas como escolha consciente |
| Sem modo reproduzível | Preset `release-reproducible` com `-O2 -ffp-contract=off`, sem flags específicas de arquitetura |
| `-Werror` com `-Ofast` nos exemplos | `-Werror` apenas em CI |

### Tarefas de medição
1. Rodar o caso Poiseuille com `-Ofast` e com `-O3`: comparar resultado numérico e
   tempo de execução. **Documentar o custo real** - se `-O3` for significativamente mais
   lento, é informação que o maintainer precisa para decidir
2. Verificar se algum resultado publicado do grupo depende de `-Ofast`; em caso
   afirmativo, registrar no PR
3. Substituir `pow(x,2)` por `x*x` nos caminhos quentes (83 chamadas com expoente
   inteiro literal) - sob `-Ofast` o GCC fazia parte disso automaticamente; ao remover
   a flag, a otimização precisa ser explícita
4. Documentar em `docs/numerics/floating-point.md` as garantias de ponto flutuante que
   o projeto oferece

### Critério de aceitação
Resultados idênticos bit a bit entre duas máquinas diferentes usando o preset
`release-reproducible`. Hoje isso é impossível por causa do `-march=native`.

### Risco
Alta sensibilidade política: mexer em flags de otimização de um código em uso para
produzir resultados publicáveis. Apresentar com medição, não com argumento teórico.

---

## E11 - Suíte de verificação numérica

**Branch:** `juniormar/11-verification` · **PR:** "Add numerical verification suite with method of manufactured solutions"

**Objetivo:** rede de proteção que torna possível refatorar sem medo. **Pré-requisito
de todo o Bloco D.**

**Depende de:** E08.

### Por que isto é o item técnico mais valioso do plano
Sem verificação numérica, qualquer refatoração é aposta. Com ela, a migração C++ passa
a ser uma sequência de mudanças cuja neutralidade é demonstrável. É também o padrão-ouro
em CFD acadêmico e um diferencial real para o projeto.

### Camadas

**1. Soluções analíticas**

| Caso | Solução exata | Verifica |
|---|---|---|
| Poiseuille 2D | `u(y) = u_max(1 − y²)` | Difusão, gradiente de pressão, Dirichlet |
| Couette | Perfil linear | Fronteira móvel |
| Stokes primeiro problema | Solução por `erfc` | Termo transiente |
| Poisson com fonte manufaturada | Escolhida | Solver linear isolado |

**2. Método das soluções manufaturadas (MMS)**

Escolhe-se `u`, `v`, `p` suaves; substitui-se nas equações de Navier–Stokes; o resíduo
vira termo-fonte. A solução exata é conhecida por construção, para qualquer modelo.
Refina-se a malha e mede-se a ordem de convergência observada - que deve coincidir com
a ordem formal do esquema.

Isto verifica **o que está implementado**, não o que se acredita ter implementado. Se o
esquema CUBISTA declara 2ª ordem e o MMS mede 1,3, há um defeito - e essa é exatamente
a classe de defeito que passa despercebida por anos num código de pesquisa.

**3. Leis de conservação**
- Massa: `∇·u` em norma discreta deve permanecer em nível de erro de máquina
- Massa em VOF: soma da fração volumétrica constante no tempo
  (a branch `PC_Daniel_Mesh` corrigiu exatamente isso - vale confrontar)
- Energia cinética em escoamento sem forçamento e sem viscosidade

**4. Regressão**
Casos de referência com saída versionada e comparação por tolerância. Detecta mudanças
não intencionais em qualquer refatoração.

**5. Testes unitários**
`higtree/atf-tests/` já existe usando o framework ATF, sem runner integrado. Avaliar:
manter ATF ou migrar para CTest. **Recomendação:** CTest, por integrar ao CMake de E08
e ao CI de E12.

### Entregáveis
`tests/verification/`, `tests/regression/`, `tests/CMakeLists.txt`,
`docs/numerics/verification.md` com as tabelas de ordem de convergência medida

### Critério de aceitação
`ctest` roda a suíte completa. Ordem de convergência medida documentada para cada
esquema temporal e convectivo. **Toda discrepância vira issue no upstream** - o que já
é contribuição científica ao projeto.

### Esforço
O maior do plano. Pode ser fatiado: primeiro Poiseuille (rápido, imediatamente útil em
E06 e E12), MMS depois.

---

## E12 - Integração contínua

**Branch:** `juniormar/12-ci` · **PR:** "Add GitHub Actions: build matrix, verification suite, container publishing"

**Objetivo:** validação automática de cada PR. **É o argumento mais forte que existe
para um maintainer aceitar contribuições** - prova objetiva de que nada quebrou.

**Depende de:** E11 (o que executar) e E04 (onde executar).

### Workflows

| Workflow | Gatilho | Conteúdo |
|---|---|---|
| `build.yml` | push, PR | Matriz: `{DIM: 2,3} × {gcc, clang} × {Debug, Release}`. 8 combinações |
| `verify.yml` | push, PR | Suíte de E11 dentro do container |
| `container.yml` | tag, push em master | Constrói e publica imagem no GHCR, multi-arquitetura |
| `lint.yml` | PR | `shellcheck` nos 35 scripts, `cmake-lint`, `clang-format --dry-run` |
| `docs.yml` | push | Doxygen (E17) publicado no GitHub Pages |
| `codeql.yml` | agendado | Análise estática de segurança |

### Detalhes
- Cache das camadas Docker e do build do PETSc - sem isso cada execução leva mais de uma hora
- Timeout e limite de concorrência (recurso gratuito é finito)
- Badges no README apontando para os workflows reais
- `ccache` no container de CI

### Critério de aceitação
Um PR de teste, contendo um defeito numérico deliberado, é reprovado pelo CI. Prova que
a suíte tem poder de detecção real, e não apenas passa sempre.

---

## E13 - Runner: scaffolding

**Branch:** `juniormar/13-scaffolding` · **PR:** "Add case scaffolding tool"

**Objetivo:** primeira das três etapas que atacam o atrito de uso. Elimina a criação
manual de arquivos, sem tocar na arquitetura.

**Depende de:** nada.

### Situação (§5 do dossiê)
Um caso exige 422 linhas de YAML + 280 a 2.609 linhas de driver C + `Makefile` +
`CMakeLists.txt` + N arquivos `.amr`. No caso Newtoniano há **uma única linha de física
real** em 280 linhas de driver.

### Ferramenta

```bash
higflow new my-case --model oldroyd-b --dim 2 --template channel
```

Gera a estrutura completa a partir de templates, com os valores padrão do modelo
escolhido e **apenas** as seções de YAML relevantes - não as 261 linhas com todos os
modelos.

### Templates iniciais
`channel` (Poiseuille), `cavity` (lid-driven), `contraction` (4:1), `bubble` (VOF),
`custom` (mínimo)

### Implementação
Python 3 (já é dependência do projeto), sem bibliotecas externas além da stdlib.

### Critério de aceitação
`higflow new` para `make` para `make run` produz VTK, sem edição manual de arquivo algum.

---

## E14 - Runner: driver unificado

**Branch:** `juniormar/14-unified-driver` · **PR:** "Add unified solver driver replacing per-case main()"

**Objetivo:** eliminar a duplicação de `main()` nos 11 drivers.

**Depende de:** E13.

### Trabalho
1. Extrair o `main()` comum (≈120 linhas idênticas nos 11 exemplos) para
   `higflow/src/hig-flow-driver.c`
2. Um único executável `higflow-solver` que carrega tudo do YAML
3. Os exemplos passam a fornecer **apenas** as funções de física realmente específicas
4. Sub-comandos: `higflow run`, `higflow validate`, `higflow info`, `higflow mesh`
5. **`higflow validate`** - validação do YAML com mensagem de erro compreensível
   apontando linha e campo. Hoje um YAML malformado produz falha obscura em tempo de
   execução
6. Fatiar `higflow_load_all_controllers_and_parameters_yaml()` (defeito C8: 1.580
   linhas numa função) em uma função por seção do YAML. Requisito prático para que a
   validação seja implementável

### Redução esperada
Driver de 280 linhas para cerca de 30, contendo só física.

### Critério de aceitação
Os 11 exemplos existentes rodam pelo driver unificado com resultados **idênticos** aos
atuais (verificado pela suíte de E11).

---

## E15 - Runner: avaliador de expressões

**Branch:** `juniormar/15-expressions` · **PR:** "Add runtime expression evaluation for boundary and initial conditions"

**Objetivo:** o alvo final da simplificação - **mudar uma condição de fronteira sem
recompilar**.

**Depende de:** E14.

### Antes e depois

Hoje, alterar o perfil de entrada exige editar C e recompilar:
```c
case 0: value = 1.5*(1.0 - center[1]*center[1]); break;
```

Depois, editar o YAML e reexecutar:
```yaml
boundaries:
  - id: 0
    type: dirichlet
    velocity:
      u: "1.5 * (1 - y^2)"
      v: "0"
```

### Implementação
Avaliador de expressões próprio, em C, sem dependência externa:
- Tokenizador e parser de precedência (*shunting-yard*)
- Compilação para bytecode na inicialização; a avaliação por ponto é interpretação de
  bytecode, não reparse
- Variáveis: `x`, `y`, `z`, `t`, `n` (passo), e parâmetros nomeados do caso
- Funções: `sin cos tan exp log sqrt abs pow min max tanh erf erfc atan2 floor ceil`
- Condicional ternário para perfis por partes
- **Requisito de desempenho:** a avaliação ocorre por ponto de fronteira e por passo de
  tempo. Precisa custar dezenas de nanossegundos. Bytecode resolve; parse a cada
  chamada, não

### Por que não usar biblioteca pronta
Adicionar dependência a um projeto cuja instalação já é o maior obstáculo seria
contraproducente. Um avaliador de expressões aritméticas é ~600 linhas de C bem
testado, e o projeto ganha um componente sem custo de dependência.

### Compatibilidade
As funções `get_*` em C continuam funcionando. Expressão em YAML é opcional e tem
precedência quando presente. Nenhum caso existente quebra.

### Critério de aceitação
Poiseuille roda com o perfil definido em YAML, resultado idêntico ao da versão
compilada, com sobrecarga de tempo abaixo de 1%.

---

## E16 - Gerador de malha

**Branch:** `juniormar/16-mesh-generator` · **PR:** "Add mesh generation from declarative description"

**Objetivo:** eliminar a escrita manual de `.amr`.

**Depende de:** E13.

### Problema
Formato posicional de 4 linhas, sem cabeçalho, sem validação, sem mensagem de erro
compreensível. Um domínio 3D complexo exige dezenas de arquivos escritos à mão. Há 780
`.amr` no repositório.

### Ferramenta
```yaml
mesh:
  domain:
    type: box
    bounds: [[0, 8], [-1, 1]]
    cells: [160, 40]
  refinement:
    - region: {type: box, bounds: [[3, 5], [-0.5, 0.5]]}
      level: 2
  boundaries:
    - {id: 0, side: xmin, name: inlet}
    - {id: 1, side: xmax, name: outlet}
    - {id: 2, side: ymin, name: wall}
    - {id: 3, side: ymax, name: wall}
```
`higflow mesh case.yaml` gera todos os `.amr` de domínio e de fronteira.

### Extras
- `higflow mesh --check` valida `.amr` existentes e reporta erro com linha e motivo
- `higflow mesh --preview` exporta VTK da malha para inspeção visual antes de simular
- Aproveitar `higtree/utilities/preProcessing/scripts-python/` e
  `meshConverter/dat_to_amr.py`, que já existem e estão subutilizados
- Documentar o formato `.amr` formalmente - `upstream/Kaina` tem um commit
  "amr documentation" a ser confrontado e creditado

### Critério de aceitação
Os `.amr` dos 11 exemplos são reproduzíveis a partir de descrição declarativa, byte a
byte ou equivalente semanticamente.

---

# BLOCO C - DOCUMENTAÇÃO E DIVULGAÇÃO

---

## E27 - Galeria de pessoas (opt-in)

**Branch:** `juniormar/27-contributors` (a partir de `juniormar/02-readme`) · **PR:** "Add opt-in People gallery"

**Objetivo:** seção final do README com foto e nome de quem contribuiu.

**Depende de:** E02 (edita o README reescrito).

### Por que é uma seção separada de Authors

As duas respondem perguntas diferentes, e misturá-las causaria dano:

| | Authors | People |
|---|---|---|
| Origem | `git shortlog` | escolha da pessoa |
| Completude | total para quem commitou | só quem consentiu |
| Precisa de permissão | não | **sim** |
| Sair da lista | não se aplica | a qualquer momento, sem justificar |

Foto de pessoa identificável é dado pessoal numa página pública. Quem não quiser
aparecer continua creditado em Authors, que vem do histórico e não depende de ninguém.

### O que foi entregue

- Seção `## People` no README, com a `<table>` pronta e **comentada**
- `docs/images/contributors/README.md` - requisitos de imagem, markup, e o texto de
  pedido de consentimento com os quatro pontos que ele precisa deixar concretos
- Requisitos de peso: 400×400 px, `.jpg`, abaixo de 100 KB. O repositório acabou de
  perder 76 MB em E07; uma dúzia de fotos sem otimizar devolveria parte disso

### O que **não** foi feito

Nenhuma foto foi adicionada. Nenhum nome foi inserido na grade. O preenchimento é
manual, feito pelo autor da contribuição depois de obter consentimento de cada pessoa.

### Pendência antes do PR

A seção precisa estar preenchida - ou removida - antes de o PR de E02 ou de E27 subir.
Uma seção "sendo montada" num README de upstream fica pela metade.

---

## E17 - Doxygen e documentação de API

**Branch:** `juniormar/17-api-docs` · **PR:** "Add Doxygen configuration and API documentation build"

**Depende de:** nada. **Trabalho prévio:** `upstream/Kaina` tem commits de Doxygen e um
grafo de dependências - partir dali e creditar.

### Tarefas
1. `Doxyfile` na raiz, unificando `higtree` e `higflow` (hoje há
   `higtree/doc/doxygen/HiGTreeDoxy`, não referenciado por nada)
2. Documentar as APIs públicas: `hig-flow-kernel.h` (1.438 linhas, 37 structs, 53
   typedefs) é a prioridade
3. Grafos de chamada e de inclusão (Graphviz)
4. Publicação em GitHub Pages via workflow de E12
5. Página inicial com a arquitetura: como higtree e higflow se relacionam, o caminho de
   um passo de tempo pelo código

### Critério de aceitação
Documentação navegável publicada, com a estrutura `higflow_solver` e o laço de projeção
compreensíveis sem ler o fonte.

---

## E18 - Tradução completa PT para EN

**Branch:** `juniormar/18-translation` · **PR:** "Translate documentation to English"

**Objetivo:** padronizar o projeto em inglês, conforme decidido.

**Depende de:** E02 e E17 (não traduzir o que ainda vai ser reescrito).

### Escopo - decisão registrada
Foi decidido traduzir **inclusive o material existente**. O risco foi apontado - é o
maior diff do plano e a maior chance de rejeição - e a decisão foi mantida. Mitigação:
fatiar em PRs pequenos e independentes, para que a rejeição de um não derrube os demais.

### Inventário

| Documento | Tamanho | PR |
|---|---|---|
| `higflow/doc/HiG-Flow.tutorial` | 6 KB | 18a |
| `Manual_Euler/manual_euler.tex` | 11 KB + 10 figuras | 18b |
| `Manual_Git/manual_git.tex` | 8 KB + 4 figuras | 18c |
| `higtree/doc/usuario.tex` | 20 KB | 18d |
| `higtree/doc/cientifica.tex` | 32 KB | 18e |
| Comentários em português no código | disperso | 18f |
| Nomes de identificadores em português (`multifase`, `exemple`) | disperso | 18g |

### Diretrizes
- Traduzir **e atualizar**: o tutorial afirma que o viscoelástico está "em
  implementação", o que é falso há anos
- Terminologia técnica consistente com a literatura de reologia computacional
- Manter as versões PT-BR ao lado (`*.pt-BR.*`) - o grupo é brasileiro
- 18f e 18g são os mais arriscados (diff enorme, ruído em `git blame`). Executar por
  último e, se houver resistência, abandonar sem prejuízo do resto

### Critério de aceitação
Nenhum documento voltado ao usuário exige português para ser compreendido.

---

## E19 - Governança do projeto

**Branch:** `juniormar/19-governance` · **PR:** "Add contribution guidelines, citation metadata, and license proposal"

**Depende de:** nada.

### Conteúdo

| Arquivo | Observação |
|---|---|
| `CONTRIBUTING.md` | Como montar o ambiente, convenção de commits, como abrir PR, como rodar a suíte de verificação |
| `CODE_OF_CONDUCT.md` | Contributor Covenant |
| `CITATION.cff` | **Alto valor acadêmico.** GitHub exibe botão "Cite this repository"; ferramentas de gestão bibliográfica leem automaticamente. Requer a lista de autores e a publicação de referência do grupo |
| `SECURITY.md` | Mínimo |
| `.github/ISSUE_TEMPLATE/` | bug report, feature request, **numerical issue** (específico: caso, malha, parâmetros adimensionais, esquema, o que se esperava) |
| `.github/PULL_REQUEST_TEMPLATE.md` | Com checklist de verificação numérica |
| `LICENSE` | **Proposta**, com a decisão explicitamente deferida aos donos |
| `CHANGELOG.md` | Keep a Changelog |

### Sobre a licença
O repositório não tem licença, o que deixa o status jurídico das contribuições
indefinido. A decisão foi **propor no PR** deixando a escolha aos donos.

Opções a apresentar, com o trade-off de cada uma:

| Licença | Implicação |
|---|---|
| **GPL-3.0** | Derivados permanecem abertos. Comum em software científico (FEniCS, deal.II). Impede uso em produto proprietário |
| **BSD-3-Clause** | Permissiva. Adotada por PETSc - dependência direta do projeto. Maximiza adoção industrial |
| **LGPL-3.0** | Meio-termo: uso como biblioteca em software fechado, modificações da biblioteca permanecem abertas |
| **Apache-2.0** | Permissiva com cláusula de patentes |

Recomendação a registrar no PR: **BSD-3-Clause**, por coerência com o PETSc e por
maximizar adoção - mas o texto do PR deixa claro que é sugestão, e que a escolha cabe
ao Prof. Castelo e coautores.

### Critério de aceitação
A aba "Insights para Community Standards" do GitHub fica completa. É um indicador
visível de maturidade do projeto.

---

## E20 - Divulgação

**Branch:** `juniormar/20-outreach` · **PR:** parcialmente fora do repositório

**Objetivo:** dar visibilidade ao projeto e às contribuições.

**Depende de:** E02, E06, E19.

### Ações no repositório
1. **Topics** do GitHub: `cfd`, `computational-fluid-dynamics`, `viscoelastic`,
   `rheology`, `navier-stokes`, `amr`, `mpi`, `petsc`, `finite-difference`,
   `multiphase-flow`, `non-newtonian`
2. Descrição e website do repositório
3. Social preview (imagem de card para redes sociais)
4. `docs/publications.md` - trabalhos que usaram o HigFlow, com DOI. Prova de uso real
   e é o que atrai novos usuários acadêmicos
5. Release com tag semântica e notas

### Ações fora do repositório (requerem aprovação dos donos)
1. **Zenodo** - integração GitHub/Zenodo dá DOI a cada release. Torna o software citável
2. **JOSS** (*Journal of Open Source Software*) - submissão. Revisão por pares, gera
   publicação citável. Requisitos: licença OSI, documentação, testes, guia de
   contribuição - exatamente E19, E11, E17
3. Listagem em diretórios de software científico (CFD Online, awesome-cfd)
4. Apresentação no grupo de pesquisa

### Sobre as recompensas do GitHub
As conquistas relevantes decorrem naturalmente do trabalho, e não de ação específica:

| Conquista | Como se obtém |
|---|---|
| **Pull Shark** | 2 / 16 / 128 PRs merjados. Este plano prevê ~20 PRs - a estratégia de PRs pequenos e independentes maximiza a chance de merge |
| **Quickdraw** | Fechar issue ou PR em menos de 5 minutos |
| **Galaxy Brain** | 2 / 8 / 16 respostas aceitas em Discussions |
| **Starstruck** | 16+ estrelas num repositório próprio |
| **Pair Extraordinaire** | Commits com co-autoria (útil ao incorporar trabalho de Pedro Coimbra e Kainã - o `Co-authored-by` é o crédito formal correto) |

O ponto principal: **PRs pequenos e revisáveis são merjados; PRs gigantes ficam
parados**. A estratégia de branch por etapa serve tanto à qualidade quanto a isso.

---

# BLOCO D - MIGRAÇÃO C para C++

**Pré-requisito absoluto: E11 concluída.** Sem verificação numérica, qualquer conversão
é aposta.

Avaliação de viabilidade em §7 do dossiê: a conversão é substancialmente mais tratável
do que o tamanho do código sugere - **1 colisão de palavra reservada**, 42 `malloc` sem
cast, 71 VLAs.

---

## E21 - Plano de migração

**Branch:** `juniormar/21-cpp-plan` · **PR:** documento apenas

**Objetivo:** documento técnico decidindo o quê, em que ordem e com que garantias.
Conforme decidido, o nível de execução é escolhido depois de ler este documento.

### Conteúdo
1. Inventário completo de incompatibilidades, arquivo por arquivo
2. Ordem de conversão (folhas da árvore de dependências primeiro: `coord.c`, `rect.c`,
   `utils.c`, `allocator.c`)
3. Estratégia de coexistência: `extern "C"` mantém interoperabilidade durante a
   transição, permitindo converter um arquivo por vez com o projeto sempre compilável
4. Protocolo de verificação por arquivo convertido
5. Análise de risco de regressão numérica: onde a semântica C e C++ diferem em ponto
   flutuante (promoção de tipos, ordem de avaliação, `<cmath>` versus `<math.h>` em
   sobrecargas de `float`/`double`)
6. Decisão sobre padrão: C++17 (amplo suporte em compiladores de cluster) versus C++20
   (concepts, `std::span` - mas GCC de cluster costuma ser antigo). **Recomendação:
   C++17**

---

## E22 - Nível 1: compilar como C++

**Branch:** `juniormar/22-cpp-level1` · **PR:** "Make codebase compile as C++ without semantic changes"

**Objetivo:** o projeto compila com `g++`/`clang++` mantendo o design em C.

### Trabalho
| Item | Volume |
|---|---|
| `lbal.c:1189` - identificador `new` | 1 |
| Casts explícitos em `malloc`/`calloc`/`realloc` | 42 |
| VLAs para `std::vector` ou `std::array` | 71 |
| Inicializador designado | 1 |
| `extern "C"` nos cabeçalhos públicos | 163 cabeçalhos |
| Renomear `.c` para `.cpp` | 167 arquivos, gradual |

### Regra inegociável
**Nenhuma mudança semântica.** Cada arquivo convertido produz resultado numérico
idêntico, verificado pela suíte de E11. Se mudar, reverter e investigar.

---

## E23 - Nível 2: RAII

**Branch:** `juniormar/23-cpp-level2` · **PR:** "Introduce RAII for resource management"

### Trabalho
1. **Eliminar os 308 `exit()`** - biblioteca não termina o processo do chamador.
   Substituir por exceções ou `expected`
2. Encapsular em tipos com destrutor: `FILE*`, comunicadores MPI, objetos PETSc
   (`Vec`, `Mat`, `KSP`), alocações de malha
3. Auditar as 306 chamadas a `fopen()` (defeito C10)
4. `std::unique_ptr` com deletor customizado para recursos de biblioteca C

### Ganho
Fim dos vazamentos, código testável unitariamente, erros propagáveis em vez de fatais.

---

## E24 - Nível 3: tipos fortes e templates

**Branch:** `juniormar/24-cpp-level3` · **PR:** "Parameterize dimension as template; introduce strong types"

### Trabalho
1. **`DIM` como parâmetro de template**, não macro de compilação. Elimina a necessidade
   de `libhig2d.a` e `libhig3d.a` separadas e a compilação dupla (defeito P3)
2. Tipos fortes: `Point`, `Tensor`, `CellId`, `FacetId` - hoje todos `real[]` ou `int`,
   sem verificação
3. `std::span` para as travessias de malha
4. Hierarquia de modelos constitutivos com interface comum - hoje cada modelo é um
   arquivo de 2 a 3,7 mil linhas com estrutura repetida

### Cuidado
Este é o nível com maior risco de virar reescrita. Manter escopo restrito e verificação
constante.

---

# BLOCO E - MACHINE LEARNING

Projeto de longo prazo. Exploratório.

---

## E25 - Arquitetura de ML no HigFlow

**Branch:** `juniormar/25-ml-design` · **PR:** documento apenas

**Contexto que muda o ponto de partida:** `upstream/Kaina` já tem `nn-weights.cpp`
(inferência MLP), modelos `.pt` versionados, e o commit mais recente (2026-08-05)
**removeu a dependência de LibTorch**. Entender por que é o primeiro passo.

### Conteúdo do documento
1. Auditoria do que existe em `upstream/Kaina`: o que `nn-weights.cpp` faz, como os
   modelos são consumidos, por que LibTorch foi removido (provável: peso da dependência
   e dificuldade de build em cluster - exatamente o problema que este plano combate)
2. **Onde ML faz sentido no HigFlow**, com fundamentação:

| Aplicação | Fundamentação | Maturidade |
|---|---|---|
| Fechamento de modelo constitutivo | Aprender a relação tensão-deformação a partir de dados experimentais, para fluidos sem modelo analítico adequado | Pesquisa ativa na literatura |
| Escolha adaptativa de passo de tempo | Prever o maior `dt` estável, substituindo heurística de CFL conservadora | Aplicável |
| Precondicionador aprendido | Acelerar a convergência do solver de pressão - o gargalo dominante | Pesquisa recente |
| Critério de refinamento AMR | Prever onde refinar, em vez de gradiente heurístico | Bom encaixe com higtree |
| Modelo substituto (*surrogate*) | Prever resultado sem simular, para varredura de parâmetros | Alto valor prático |

3. **Opções em C++ puro**, sem dependência pesada: inferência própria (o caminho da
   Kaina), ONNX Runtime, `frugally-deep`, `tiny-dnn`
4. Requisito arquitetural: inferência **por célula por passo de tempo** exige latência
   de microssegundos. Isso descarta qualquer solução que serialize ou chame Python
5. O problema da verificação: um modelo constitutivo aprendido não tem solução
   analítica. Como verificar? Ligação direta com E11

---

## E26 - Prova de conceito

**Branch:** `juniormar/26-ml-poc` · **PR:** "Add standalone MLP inference module"

**Depende de:** E25.

### Escopo deliberadamente mínimo
1. Inferência MLP em C++ puro, sem dependências: camadas densas, ReLU/tanh/sigmoid,
   leitura de pesos de arquivo
2. Formato de pesos simples e documentado, com script Python de exportação a partir de
   PyTorch
3. Testes unitários comparando saída contra a referência do PyTorch, com tolerância
4. Benchmark de latência por inferência
5. **Não integra ao solver.** Módulo isolado e testado, pronto para integração futura

### Critério de aceitação
Inferência bate com o PyTorch dentro de `1e-6`, com latência abaixo de 1 µs para uma
rede pequena.

---

## Resumo executivo

| Bloco | Etapas | PRs | Esforço relativo |
|---|---|---|---|
| 0 - Fundação | E00-E01 | 0 | Baixo |
| A - Apresentação | E02-E06 | 5 | Médio-alto |
| B - Código e build | E07-E16 | 10 | Alto |
| C - Documentação | E17-E20, E27 | 5+ | Médio |
| D - C++ | E21-E24 | 4 | Alto |
| E - ML | E25-E26 | 2 | Médio |
| **Total** | **28** | **~26** | - |

### Sugestão de ordem para as primeiras sessões

1. **E01** - fechar a infraestrutura (workflow, ferramenta de status, `.mailmap`)
2. **E02** - README v1: entrega visível imediata, sem depender de ambiente
3. **E00** - ambiente WSL2 em paralelo (a compilação do PETSc roda sozinha)
4. **E07** - higiene: PR pequeno, fácil de aceitar, efeito grande
5. **E04** - containers: a etapa de maior valor por esforço
6. **E06** - galeria: fecha o bloco de apresentação com o README completo

A partir daí, qualquer ordem. Cada etapa é autocontida.

---

*Estado de execução em [`LOG.md`](LOG.md). Verificação de estado:
`bash tools/contrib/status.sh`.*
