# Modelo de Branches, Commits e Pull Requests

> Documento de trabalho interno. Não integra Pull Requests.

---

## 1. Topologia dos remotos

```
upstream   https://github.com/antoniocastelofilho/HigFlow.git   (leitura; repositório oficial)
origin     https://github.com/juniormarorganista/HigFlow.git    (seu fork; escrita)
```

Fluxo dos Pull Requests: `origin/juniormar/NN-slug` para `upstream/master`.

---

## 2. Topologia das branches

```
upstream/master ──●-----------------------------------------------►
                  │
                  ├── juniormar/02-readme ──●──●──●──┐   PR #1
                  │                                  │
                  ├── juniormar/04-containers ─●──●──┤   PR #2
                  │                                  │
                  ├── juniormar/07-repo-hygiene ─●---┤   PR #3
                  │                                  │
                  └── juniormar/main ◄---------------┘
                        (integração local; nunca vira PR)
```

### As três categorias de branch

| Branch | Origem | Destino | Papel |
|---|---|---|---|
| `master` | `upstream/master` | - | Espelho do upstream. **Nunca receber commits.** |
| `juniormar/main` | `master` | nenhum | Tronco de integração local. Reúne o planejamento e o merge de todas as etapas. Permite ver o estado combinado do trabalho |
| `juniormar/NN-slug` | **`master`** | PR para `upstream/master` | Uma etapa, um PR. Contém **apenas** os commits daquela etapa |

### A regra que garante PRs limpos

As branches de etapa saem de **`master`**, não de `juniormar/main`.

Motivo: `juniormar/main` contém o diretório `docs/contrib-plan/`, que é material de
trabalho pessoal e não deve aparecer em nenhum PR. Ramificando de `master`, cada branch
de etapa contém exclusivamente o que pertence àquela etapa - o diff do PR é limpo por
construção, sem necessidade de `cherry-pick` ou rebase interativo na hora de enviar.

`juniormar/main` recebe o merge de cada etapa concluída, apenas para consulta local do
estado agregado.

### Etapas com dependência real

Quando uma etapa depende tecnicamente de outra (ver grafo em `01-ROADMAP.md` §Mapa),
a branch sai da etapa da qual depende, e não de `master`:

```bash
# E06 (galeria) depende de E04 (containers)
git checkout juniormar/04-containers
git checkout -b juniormar/06-gallery
```

No GitHub, o PR de E06 é aberto **contra a branch de E04** enquanto E04 estiver
pendente. Depois que E04 for merjada no upstream, o PR de E06 é reapontado para
`master` - o GitHub faz isso automaticamente ao detectar o merge da base.

---

## 3. Ciclo de trabalho de uma etapa

```bash
# 1. Sincronizar com o upstream antes de qualquer coisa
git fetch upstream
git checkout master
git merge --ff-only upstream/master

# 2. Criar a branch da etapa a partir de master
git checkout -b juniormar/07-repo-hygiene

# 3. Trabalhar, com commits frequentes e pequenos
#    (ver §4)

# 4. Conferir o estado antes de considerar concluída
bash tools/contrib/status.sh

# 5. Integrar no tronco local para ver o estado agregado
git checkout juniormar/main
git merge --no-ff juniormar/07-repo-hygiene

# 6. Quando decidir enviar (e só então):
git push -u origin juniormar/07-repo-hygiene
gh pr create --repo antoniocastelofilho/HigFlow \
             --base master \
             --head juniormarorganista:juniormar/07-repo-hygiene
```

**Nada é enviado ao `origin` sem decisão explícita.** O trabalho permanece local até
que a etapa esteja pronta e revisada.

---

## 4. Convenção de commits

Formato Conventional Commits, em inglês.

```
<tipo>(<escopo>): <assunto no imperativo, minúscula, sem ponto final>

<corpo: o porquê, não o quê - o diff já mostra o quê>

<rodapé: Refs, Co-authored-by>
```

### Tipos

| Tipo | Uso |
|---|---|
| `feat` | Nova funcionalidade |
| `fix` | Correção de defeito |
| `perf` | Melhoria de desempenho |
| `refactor` | Mudança sem alteração de comportamento |
| `build` | Sistema de build, dependências |
| `ci` | Integração contínua |
| `docs` | Documentação |
| `test` | Testes |
| `chore` | Manutenção, higiene do repositório |
| `style` | Formatação apenas |

### Escopos do projeto

`higtree`, `higflow`, `build`, `cmake`, `install`, `container`, `docs`, `examples`,
`ci`, `runner`, `mesh`, `numerics`, `vof`, `viscoelastic`, `eo` (eletro-osmótico),
`contrib`

### Exemplos

```
fix(higtree): remove unconditional DEBUG define from public header

hig-flow-kernel.h defined DEBUG unconditionally, forcing every
translation unit that includes it into the debug path regardless of
the build type. This defeated NDEBUG in release builds and pulled
per-TU static state from Debug-c.h into every object file.

DEBUG is now controlled by the build system.
```

```
build(cmake): link BLAS and LAPACK that were previously discarded

find_package(BLAS) and find_package(LAPACK) appended to
MFSIM_DEPENDENCIES, a variable belonging to a different project and
never referenced here. Both libraries were located and silently
dropped.
```

### Regras

1. **Um commit, uma ideia.** Se a mensagem precisa de "e", provavelmente são dois commits
2. **Muitos commits pequenos**, conforme solicitado. Um PR de 15 commits legíveis é
   revisado; um de 1 commit gigante, não
3. **Nunca mencionar ferramentas de IA** em nenhuma parte do commit
4. **Sem `Co-authored-by` automático.** Usar apenas para crédito real de pessoas -
   por exemplo, ao incorporar trabalho de `PC_ImproveDocumentation` ou `Kaina`
5. O corpo explica **por que**, incluindo o defeito concreto que motivou a mudança.
   Num repositório acadêmico, o revisor precisa entender a motivação sem ter feito a
   análise

### Crédito a trabalho incorporado

Ao trazer material das branches do upstream:

```
Co-authored-by: Pedro Coimbra <pedro.coimbra@wikki.com.br>
```

Além disso, o corpo do commit e a descrição do PR declaram explicitamente a origem.
Isso é correção de atribuição e é o que torna a contribuição bem recebida.

---

## 5. Estrutura de um Pull Request

```markdown
## Summary
<uma ou duas frases: o que muda e por quê>

## Motivation
<o defeito ou lacuna concreta, com arquivo:linha quando aplicável>

## Changes
- <lista objetiva>

## Verification
<como foi testado: comandos, plataforma, resultado>

## Based on previous work
<quando aplicável: qual branch, qual autor, o que foi aproveitado e o que foi corrigido>

## Notes for reviewers
<pontos que exigem decisão do maintainer, e não do contribuidor>
```

### Diretrizes

- **PRs pequenos.** Alvo: menos de 400 linhas de diff efetivo. Acima disso, fatiar
- **Um assunto por PR.** Correção de build e tradução de documentação não convivem
- **Nada de mudança oportunista.** Reformatar código não relacionado destrói o
  `git blame` e é a causa mais comum de rejeição
- **Decisões que cabem aos donos** (licença, remoção de `src_hugo/`, escolha de esquema
  padrão) são apresentadas como pergunta na seção *Notes for reviewers*, nunca decididas
  no diff

---

## 6. Verificação de estado

`bash tools/contrib/status.sh` responde, num único comando:

1. Branch atual, alterações não commitadas
2. Quais etapas do roadmap têm branch, quais têm commits, quais foram integradas ao tronco
3. O que existe localmente e **não** foi enviado ao `origin`
4. O que o `origin` tem e o `upstream` não - candidatos a PR
5. O que o `upstream` ganhou desde o último `fetch` - risco de conflito
6. Estado dos PRs abertos, se o `gh` estiver autenticado

Executar sempre **antes** de iniciar uma etapa e **antes** de abrir um PR.

---

## 7. Sincronização com o upstream

O upstream tem branches ativas (`Kaina` com atividade em agosto de 2026). Antes de
qualquer etapa que toque em build, `domain.c`, `higtree.h` ou Makefiles:

```bash
git fetch upstream --prune
bash tools/contrib/status.sh
```

Se o upstream avançou em arquivos que a etapa toca, decidir conscientemente entre:
- Rebase da branch de etapa sobre o novo `master`
- Ajustar o escopo da etapa
- Adiar a etapa

Etapas de maior risco de conflito, por ordem: **E08** (build), **E09** (correções em
higtree), **E18** (tradução, diff amplo).

---

## 8. Higiene do histórico

| Situação | Ação |
|---|---|
| Commit com erro de digitação na mensagem, ainda não enviado | `git commit --amend` |
| Vários commits de tentativa e erro numa etapa | `git rebase -i` para consolidar **antes** do push |
| Já enviado ao `origin` e com PR aberto | **Não reescrever.** Novo commit de correção |
| Branch de etapa merjada no upstream | Apagar local e remota |

O `git rebase -i` não é executável neste ambiente (é interativo). Consolidação de
commits, quando necessária, é feita manualmente com `git reset --soft` seguido de novo
commit.

---

## 9. O que nunca subir

| Item | Motivo |
|---|---|
| `docs/contrib-plan/` | Material de trabalho pessoal, em português |
| `tools/contrib/status.sh` | Ferramenta pessoal de acompanhamento |
| Binários compilados | Ver E07 |
| Saídas de simulação (VTK, DATA) | Exceto figuras curadas da galeria |
| Caminhos absolutos da máquina local | Quebram para qualquer outra pessoa |
| Qualquer menção a ferramentas de IA | Decisão registrada |

O diretório `docs/contrib-plan/` só existe em `juniormar/main`. Como as branches de
etapa saem de `master`, ele nunca aparece em um PR - a proteção é estrutural, não
depende de lembrar.
