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
| E00 — Ambiente de desenvolvimento | — | **pendente** | — |
| E01 — Infraestrutura da contribuição | `juniormar/main` | **concluída** (falta `.mailmap`, que vai no PR de E07) | — |
| E02 — README em inglês | — | pendente | — |
| E03 — Instalação Linux | — | pendente | — |
| E04 — Containers | — | pendente | — |
| E05 — Instalação Windows | — | pendente | — |
| E06 — Galeria de resultados | — | pendente | — |
| E07 — Higiene do repositório | — | pendente | — |
| E08 — Build unificado | — | pendente | — |
| E09 — Correções pontuais | — | pendente | — |
| E10 — Flags numéricas | — | pendente | — |
| E11 — Verificação numérica | — | pendente | — |
| E12 — Integração contínua | — | pendente | — |
| E13 — Runner: scaffolding | — | pendente | — |
| E14 — Runner: driver unificado | — | pendente | — |
| E15 — Runner: expressões | — | pendente | — |
| E16 — Gerador de malha | — | pendente | — |
| E17 — Doxygen | — | pendente | — |
| E18 — Tradução | — | pendente | — |
| E19 — Governança | — | pendente | — |
| E20 — Divulgação | — | pendente | — |
| E21 — Plano C++ | — | pendente | — |
| E22 — C++ nível 1 | — | pendente | — |
| E23 — C++ nível 2 | — | pendente | — |
| E24 — C++ nível 3 | — | pendente | — |
| E25 — Arquitetura ML | — | pendente | — |
| E26 — PoC ML | — | pendente | — |

---

## 2026-08-20 — Análise e planejamento

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
   autoriza o compilador a assumir que `Inf` e `NaN` nunca ocorrem — mas
   `hig-flow-step-electroosmotic.c:923,932` usa `INFINITY` como sentinela em laço de
   convergência, e a detecção de divergência depende de comparar `NaN`. É problema
   metodológico, não apenas de engenharia.

5. **Os dois sistemas de build compilam conjuntos diferentes de arquivos.** Nenhum
   compila o projeto inteiro. Provável raiz de boa parte da dificuldade relatada de
   fazer o projeto funcionar.

6. **101 commits do autor estão fragmentados em 4 identidades git.** Um `.mailmap` de
   5 linhas é a correção de maior retorno por esforço de todo o levantamento.

### Pendências desta etapa

- `.mailmap` — escrito no PR de E07, onde faz sentido junto com `.editorconfig` e
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

Sugestão registrada no roadmap: **E01 concluída → E02 (README v1)** como entrega
visível imediata, com **E00 (WSL2)** em paralelo, já que a compilação do PETSc roda
sozinha por 30 a 90 minutos.

---
