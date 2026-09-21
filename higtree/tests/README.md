# A suíte da HiGTree

Esta suíte afirma **valor**, não forma. Ela existe porque a suíte do HiGFlow
(`ci/run_suite.py`) **constrói** a HiGTree e nunca a **testa**: a biblioteca é
tratada como pré-requisito dos exemplos, e um defeito dentro dela só aparece
como número errado algumas camadas adiante — quando aparece.

Não é hipótese. Todos os defeitos corrigidos no `domain.c` entre 15 e 18 de
setembro de 2026 moravam na HiGTree e foram encontrados por uma simulação 3D
abortando, não por um teste.

A distinção entre **valor** e **forma** é o que justifica a suíte existir ao lado
da ATF de 2020. Aquela já chamava `sd_get_stencil` em pontos fora do domínio, com
contornos montados, e verificava:

```c
ATF_CHECK(stn_get_numelems(stn) >= 3);
ATF_CHECK(FLT_NE(stn_get_rhs(stn), 0.0));
```

isto é, "tem pelo menos três elementos" e "o lado direito não é zero". Os nove
defeitos de fechamento daquela semana passariam pelas duas sem exceção: todos
produziam estêncil com a forma certa e os números errados. **Uma asserção só vale
se o defeito a violaria.**

## Como rodar

```bash
python3 ci/run_higtree_tests.py
```

Precisa de `PETSC_DIR` no ambiente — rode `source ./varsrc` antes. O driver
reconstrói a biblioteca por dimensão (`make clean` entre 2D e 3D), porque os
objetos **não carregam a dimensão no nome**: depois de um build `DIM=3`, o
`libhig2d.a` simplesmente não existe, e o `make` responde `No rule to make
target`, que se lê como fonte faltando.

Para incluir os testes do t8code, `--t8code`. Sem a opção eles são pulados e a
suíte roda exatamente como antes — o t8code não é dependência da HiGTree.

Cada binário imprime uma linha por caso e um resumo:

```
caso <nome> PASS
caso <nome> FAIL <o que falhou, com os números>
resumo <n> casos, <k> falharam
```

Nesta árvore, hoje: **170 de 170 casos** e **18 de 18 cláusulas** do contrato.
Com `--t8code` e a biblioteca construída, sobe para 316 casos.

O driver **exige** a linha `resumo`. Um binário que morre no meio imprime casos
que passaram e sai — sem o resumo, isso se leria como sucesso parcial. É a
mesma família de armadilha que o `--verificar` do instalador tinha.

## O que cada teste afirma

Os cabeçalhos de cada arquivo explicam o *porquê* de cada um, muitas vezes
nomeando o defeito que o motivou. O que segue é só o mapa.

**A malha responde certo sobre si mesma**
| arquivo | afirma |
|---|---|
| `test-cell-queries.c` | centro, tamanho, identificador e iterador — as quatro consultas com mais de 850 usos na HiGFlow, até então cobertas só por acidente |
| `test-point-location.c` | a célula devolvida contém o ponto, e o desempate é **definido**, não emergente |
| `test-point-class.c` | dentro / sobre o contorno / fora — o despacho que escolhe o ramo do fechamento |
| `test-boundary-faces.c` | onde o domínio termina, dito por dois mecanismos opostos que têm de concordar |

**O estêncil produz o valor certo**
| arquivo | afirma |
|---|---|
| `test-stencil-value.c` | reprodução polinomial exata — o oráculo que não precisa de referência gravada |
| `test-stencil-support.c` | a malha fornece o suporte e o ajuste é o mesmo, testados **separados** |
| `test-stencil-selection.c` | **qual** parede fecha o estêncil, medido pelo valor, sem espiar o interior |
| `test-boundary-path.c` | o ramo `ON_BOUNDARY`, que nenhum dos outros alcança — foi ali que três dos sete sítios de um defeito sobreviveram, verdes por ausência de teste |
| `test-level-jump.c` | estêncil atravessando salto de nível maior que 2:1, em malha não graduada |

**A partição não muda a resposta**
| arquivo | afirma |
|---|---|
| `test-partition-independence.c` | a mesma geometria, descrita de duas formas, dá os mesmos valores |
| `test-fringe-support.c` | o que a franja entrega ao estêncil — e o que ela **não** entrega ao iterador |
| `test-fringe-sync.c` | a troca de franja entrega o **valor** certo, não só a estrutura |
| `test-fringe-parallel.c` | a franja sob particionamento real |
| `test-facet-domain-serial.c` | a montagem serial de um domínio de facetas ainda funciona (ver abaixo) |

**O instantâneo descreve a malha que existe**
| arquivo | afirma |
|---|---|
| `test-mesh-snapshot.c` | o instantâneo e os dois backends que o preenchem |
| `test-domain-snapshot.c` | quando ele nasce, e o que acontece se a malha mudar depois |
| `test-snapshot-oracle.c` | que o oráculo diferencial **acusa** — um detector que nunca acusou não é detector |
| `test-mesh-snapshot-parallel.c` | o mesmo sob particionamento real |

**O t8code cumpre o contrato** (só com `--t8code`)
| arquivo | afirma |
|---|---|
| `test-partition-t8code.c` | o t8code sob particionamento real |
| `test-production-t8code.c` | o t8code **produzindo** a malha que o domínio usa |
| `test-rank-production-t8code.c` | produção por rank, em caixas completas |
| `test-particao-t8code.c` | a partição final sendo a do t8code, com a franja como oráculo |

O contrato que essas cláusulas verificam está em `higtree/src/hig-mesh-contract.h`
— 14 cláusulas `C` e 4 `P`. O driver lê as cláusulas de lá e exige que cada uma
tenha caso passando: uma cláusula que perdeu o teste sai como `SEM TESTE QUE
RODE`, uma cujo teste falhou sai como `REPROVADA`.

## O que o verde NÃO prova

Esta é a parte que o número no fim da execução não diz.

**Nenhuma física roda aqui.** A suíte cobre malha, localização, estêncil e
partição. Não há Navier-Stokes, não há passo no tempo, não há sistema linear.
Verde aqui significa que a camada de malha responde certo; não diz nada sobre a
discretização que se apoia nela. Quem cobre aquilo é o `ci/run_suite.py`, e
aquilo tem os próprios pontos cegos.

**O paralelismo vai até np = 3.** Sete dos testes rodam em np = 1, 2 e 3;
dezesseis chamam MPI. Nada exercita contagem alta de ranks, e o `MAXPARTITIONS`
é 256. Um defeito que só aparece com muitos vizinhos não seria visto.

**Nada aqui chega perto do teto de 40 árvores por domínio.** O `sfd->sfbi[]` é
arranjo fixo de `MAXHIGTREESPERDOMAIN` e é indexado sem conferência (ver o
cabeçalho do `domain.h`). As malhas destes testes são pequenas; o caso do
repositório que mais se aproxima é o `example3d_complex`, com 33 blocos, e ele
não roda aqui.

**Os oráculos são analíticos, e isso corta nos dois sentidos.** Reprodução
polinomial e conservação em corte transversal não dependem de referência
gravada — é a força desta suíte. Mas um defeito que **preserva** a reprodução
polinomial é invisível para ela. Foi por isso que o `test-stencil-selection`
precisou ser escrito com paredes de valores diferentes: a seleção errada de
condição de contorno ainda reproduzia polinômios.

**O VTK não é testado aqui, e não deve ser.** Os escritores VTK são
dependentes de partição por construção. Ver a seção de VTK no `AGENTS.md`.

**A suíte ATF de 2020 não roda.** São 71 casos em `../atf-tests` que precisam
de `atf-c`, não instalado em lugar nenhum — cinco anos sem executar. Ficam como
material de referência, não como suíte.

## Uma correção que vale registrar

Até 2026-09-19 o `sfd->sfbi[]` só era preenchido pelo caminho **particionado**,
então um `sim_facet_domain` montado sem MPI estourava na primeira consulta de
faceta — metade da interface de malha não podia ser exercitada sem a camada de
partição. Hoje são duas metades: `sfd_compute_sfbi` (no `domain.c`) preenche os
blocos das árvores locais e não precisa de MPI; `psfd_compute_sfbi` chama essa e
depois troca os blocos de franja com os vizinhos.

É a cláusula P3 do contrato, e o `test-facet-domain-serial.c` existe para que a
regressão não volte em silêncio. Se você encontrar comentário afirmando que não
há caminho serial, é anterior a essa correção.
