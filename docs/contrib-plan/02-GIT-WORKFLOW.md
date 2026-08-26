# Modelo de Branches, Commits e Pull Requests

> Documento de trabalho interno. Vive na branch `juniormar/notes` e nunca integra
> Pull Requests.

---

## 1. Topologia

```
upstream/master  ────────────────────────────────────────────────────────►
      │
      └── juniormar/main ──●──────●──────●──────●──────●───────────────►
                            \    / \    / \    / \    / \    /
                             \  /   \  /   \  /   \  /   \  /
                              pr1    pr2    pr3    pr4    pr5   ...

juniormar/notes  ────●───●───●───●─────────────────────────────────►
                     (só o material de planejamento)
```

| Branch | Papel |
|---|---|
| `master` | Espelho do `upstream/master`. Nunca recebe commit |
| `juniormar/main` | **A branch principal do trabalho.** Sai da `master` e cresce por merge de cada etapa concluída. É dela que saem as branches de etapa e para ela que elas voltam |
| `juniormar/<etapa>` | Uma etapa, um Pull Request. Sai da `juniormar/main` e volta por merge quando pronta |
| `juniormar/notes` | Planejamento, análise, log e ferramentas de acompanhamento. Só isso. Nunca vira PR |

## 2. Por que o planejamento fica separado

A `juniormar/main` é a base das branches de etapa, e as branches de etapa viram Pull
Request. Tudo que estiver na `main` aparece no diff de todo PR que sair dela.

O material de planejamento são cerca de 2.400 linhas de notas internas em português.
Num PR para um repositório acadêmico isso dilui a revisão e atrapalha a aprovação, que
é o objetivo. Por isso vive na `juniormar/notes`, que nunca é base de nada.

Para consultar o plano sem trocar de branch:

```bash
git show juniormar/notes:docs/contrib-plan/LOG.md
git show juniormar/notes:docs/contrib-plan/01-ROADMAP.md
```

Para editar, troque de branch, edite, commite e volte.

## 3. Ciclo de uma etapa

```bash
# 1. sincronizar
git fetch upstream
git checkout master && git merge --ff-only upstream/master

# 2. atualizar a main com o que o upstream aceitou, se houve merge lá
git checkout juniormar/main && git merge master

# 3. abrir a etapa a partir da main
git checkout -b juniormar/e11-verification

# 4. trabalhar, com commits pequenos

# 5. quando pronta, integrar na main
git checkout juniormar/main
git merge --no-ff juniormar/e11-verification -m "merge: E11 numerical verification"

# 6. quando decidir enviar
git push -u origin juniormar/e11-verification
gh pr create --repo antoniocastelofilho/HigFlow --base master \
             --head juniormarorganista:juniormar/e11-verification
```

**Nada vai para o `origin` sem decisão explícita.**

## 4. A consequência de a main acumular tudo

Como a `juniormar/main` contém todas as etapas ainda não aceitas pelo upstream, uma
etapa nova que sair dela carrega as anteriores. O diff do PR mostra tudo que ainda não
foi merjado lá.

Isso não é um defeito do modelo, é o que significa trabalhar num fork à frente do
upstream. Duas coisas o tornam administrável:

- Cada PR diz no corpo quais commits são dele e quais vieram dos anteriores
- Cada merge no upstream encolhe automaticamente os PRs seguintes

A alternativa seria abrir cada etapa a partir da `master` pura, o que dá diffs
mínimos mas volta a permitir conflito entre etapas que tocam o mesmo arquivo. Foi
justamente isso que precisou ser desfeito em agosto de 2026.

## 5. Convenção de commits

Conventional Commits, em inglês.

```
<tipo>(<escopo>): <assunto no imperativo, minúscula, sem ponto final>

<corpo: o porquê, não o quê>

<rodapé: Refs, Co-authored-by>
```

Tipos: `feat`, `fix`, `perf`, `refactor`, `build`, `ci`, `docs`, `test`, `chore`,
`style`.

Escopos: `higtree`, `higflow`, `build`, `cmake`, `install`, `container`, `docs`,
`examples`, `ci`, `runner`, `mesh`, `numerics`, `vof`, `viscoelastic`, `eo`, `gallery`,
`notes`.

Regras:

1. Um commit, uma ideia. Se a mensagem precisa de "e", provavelmente são dois
2. Muitos commits pequenos: um PR de 15 commits legíveis é revisado, um de 1 gigante não
3. Nunca mencionar ferramentas de IA em nenhuma parte
4. `Co-authored-by` só para crédito real de pessoas
5. O corpo explica **por que**, incluindo o defeito concreto que motivou a mudança

## 6. Escrita

Sem travessão longo, sem seta como conectivo, sem régua decorativa em comentário, sem
emoji, sem meia-risca em intervalo. Preservado por ser notação: meia-risca entre nomes
de pessoas diferentes (Navier-Stokes), sinais de multiplicação e menos, glifos que
desenham árvore de diretório.

## 7. O que nunca sobe num PR

| Item | Motivo |
|---|---|
| `docs/contrib-plan/`, `tools/contrib/` | Material de trabalho pessoal, em português. Vive na `notes` |
| Binários compilados | Ver E07 |
| Saídas de simulação | Exceto figuras curadas da galeria |
| Caminhos absolutos da máquina local | Quebram para qualquer outra pessoa |
| Menção a ferramentas de IA | Decisão registrada |

## 8. Estado atual

| Branch | Papel | Enviada |
|---|---|---|
| `juniormar/main` | integração, base das etapas | sim |
| `juniormar/notes` | planejamento | sim |
| `juniormar/pr1-repo-hygiene` | E07 | PR #3 |
| `juniormar/pr2-readme` | E02 e E27 | PR #4 |
| `juniormar/pr3-containers` | E04 | PR #5 |
| `juniormar/pr4-windows-guide` | E05 | PR #6 |
| `juniormar/pr5-gallery` | E06 | PR #7 |

As cinco branches de PR foram criadas antes desta reorganização, saindo da `master` em
cadeia linear. Como a `juniormar/main` foi reconstruída merjando exatamente essa
cadeia, o resultado é idêntico ao que o modelo novo produziria, e nenhuma delas
precisou ser reescrita. Os PRs não foram tocados.
