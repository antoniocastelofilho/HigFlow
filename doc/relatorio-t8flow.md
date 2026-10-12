# t8Flow — relatório de avanços

Estado em **11/10/2026**. Branch `t8flow`, worktree
`/home/castelo/HiGFlow-System/t8flow-wt`. Dois commits: `ea83f27` (esqueleto) e
`11239e0` (solver e discretizador 2D).

    suite do lado 2d    15 casos, 0 falharam
    contrato            9/9 clausulas (N1..N6, B1, F1, F2)
    codigo de saida     0

O projeto completo está em [`projeto-t8flow.md`](projeto-t8flow.md); este
relatório é o resumo do que avançou e do que está medido.

---

## 1. As sete decisões, na ordem em que foram tomadas

Todas do usuário, todas em 11/10/2026. A ordem importa porque cada uma restringe
as seguintes.

| | decisão | consequência imediata |
|---|---|---|
| D1 | herdar o núcleo independente de malha da HiGTree — solver e MLS | 4822 das 21581 linhas de `higtree/src` entram sem porte |
| D2 | partição e domínio do t8code | saem `lbal.c`, `build-fringe.cpp` e metade do `pdomain.c` |
| D3 | a suíte nasce com o código | contrato cobrado pelo driver desde o primeiro commit |
| D4 | um domínio só | some o domínio uniforme de 375 mil células |
| D5 | C++ na norma do t8code (`gnu++20`) | uma norma só; nada de `extern "C"` |
| D6 | vive em worktree, branch `t8flow` | isolado das sessões na árvore principal |
| D7 | código separado 2D/3D, começando pelo 2D | `src/comum`, `src/2d`, `src/3d`; sem `#if DIM` no código novo |

E a ordem de trabalho, também do usuário: **solver e discretizador primeiro**. Foi
a escolha certa — são o que a D1 já garantiu existir, adaptá-los é fachada e não
formulação, e o resultado é testável sem malha nenhuma.

---

## 2. O que está medido

Nada nesta seção é estimativa.

### 2.1 O núcleo herdado é de fato independente de malha

Cada unidade de tradução compilada com os cabeçalhos de malha substituídos por
arquivos vazios — e o teste teve de ser refeito, porque a primeira versão
compilava no próprio diretório do fonte e `#include "x.h"` procura ali antes de
qualquer `-I`: ela declarou `domain.c` independente de malha.

| | resultado |
|---|---|
| referências a `hig_cell`/`higcit`/`sim_domain`/`psd_` em todos os 7 backends do solver | **zero** |
| `#include "higtree.h"` em `solver.h` | vestigial, nada usado |
| dependências de `wls.h` | só `coord.h` e `types.h` |
| núcleo independente | **4822 de 21581 linhas** (~22%) |

### 2.2 A norma do t8code não custou porte

O t8code **impõe** a norma ao consumidor:
`target_compile_features(T8 PUBLIC c_std_11 cxx_std_20)`, que reaparece como
`INTERFACE_COMPILE_FEATURES` na interface instalada.

| | gnu++17 | gnu++20 |
|---|---|---|
| 16 unidades do núcleo herdado | ok | **ok, sem uma alteração** |
| as mesmas, sem `-fpermissive` | — | 15 limpas; só `solver-sor.c` falha (11 conversões `void*`) |

Daí a recomendação adotada: compilar sem `-fpermissive`, e resolver o SOR com
casts ou deixá-lo de fora (é backend opcional; o PETSc é o de produção).

### 2.3 Por que um domínio só

O domínio uniforme do VOF é a malha cheia por construção —
`higflow_create_amr_info_mult` copia os limites inteiros e usa o `delta` do
último nível.

| | células |
|---|---|
| `dpfracvol`, domínio uniforme | 375 000 (50×50×150) |
| `dpu2`, árvore adaptada | ~81 000 |

4,6× com **um** nível de refino, e a razão cresce 8× por nível em 3D. O
escalonamento se aplicaria à parte que não domina o custo.

### 2.4 A franja não pode ser o ghost do t8code

| | profundidade | vizinhança |
|---|---|---|
| ghost do t8code | 1 camada | só `T8_GHOST_FACES` ("currently only T8_GHOST_FACES is supported") |
| franja do HiGFlow | 5 células | `pg_set_fringe_size(pg, 5)` |

Não é limitação do t8code: ghost de uma camada é o que volumes finitos clássico
pede, e o suporte largo é consequência do MLS que o t8Flow herda por decisão. A
franja vai por caixas, com profundidade pedida.

---

## 3. O código que existe

    t8flow/Makefile                       alvos `2d` e `3d`, gnu++20, sem -fpermissive
    t8flow/src/comum/arranque.{hxx,cxx}   a ordem de arranque, num objeto
    t8flow/src/comum/t8flow-contract.h    9 clausulas
    t8flow/src/2d/mls2d.{hxx,cxx}         o discretizador
    t8flow/src/2d/solver2d.{hxx,cxx}      o solver
    t8flow/tests/testing.h                arcabouco, zero dependencia externa
    t8flow/tests/2d/test-mls2d.cxx        5 casos
    t8flow/tests/2d/test-solver2d.cxx     6 casos
    t8flow/tests/2d/test-floresta2d.cxx   4 casos
    ci/run_t8flow_tests.py                driver, com quatro guardas

**O discretizador** acrescenta ao herdado tipos 2D concretos (`Ponto2` é
`std::array<double,2>`, com a travessia sustentada por `static_assert`), tempo de
vida, e a condição de uso afirmada: abaixo do mínimo de pontos para a ordem,
`pesos()` recusa em vez de aproximar.

**O solver** concentra quatro armadilhas do herdado: tempo de vida; a ordem de
montagem (só existe `monta()`, e não há como fazer a ordem errada que custou um
SIGSEGV); o nome do motor (a aritmética do deslocamento entre as duas listas do
`solver.h` em um lugar, cobrada por teste); e toda linha preenchida — linha
faltando não é sistema singular, é abort dentro da pré-alocação do PETSc.

**O contrato** não é prosa: o driver lê as cláusulas do arquivo e exige, para
cada uma, um caso que saiu `PASS`. Cláusula sem cobertura, com caso ausente ou
reprovado conta como falha. Uma cláusula entra quando tem caso que a afirme,
nunca antes — assim o vermelho continua significando "quebrou".

---

## 4. Os controles

É a seção que dá sentido ao verde. Cada guarda foi removida para ver o que
acontece.

| sabotagem | resultado |
|---|---|
| ordem do ajuste 2 → 1 | N1 **reprova** nas duas metades |
| centroides do suporte deslocados meia célula | F2 **continua verde** — limitação conhecida |
| cláusula apontando caso inexistente | driver: `CASO AUSENTE`, contrato 4/5, saída 1 |
| guarda de suporte mínimo removida | N3 **reprova**: "suporte de 5 pontos (mínimo 6) passou sem reclamar" |
| guarda de linha faltando removida | o processo **morre** no abort do herdado, sem linha de resumo |

A segunda linha é a mais informativa: ela mede que a reprodução polinomial
verifica a **discretização** e não a malha — por isso "o suporte cai na grade e é
local" precisa ser cláusula própria (F3), não subproduto.

A última justifica as fachadas existirem: elas convertem um abort três camadas
abaixo em exceção que nomeia a linha que falta. E foi a guarda do driver que
exige a linha de resumo que pegou o processo morto, em vez de reportar "nenhuma
falha".

**Dois oráculos que passavam por acidente, pegos pelo próprio contraste.** Com o
ponto de consulta no centro de simetria da nuvem, o ajuste de ordem 2 reproduziu
um campo **cúbico** com erro 1,1e-16 — os termos ímpares cancelam contra pesos
simétricos. O caso reprovou dizendo "o oráculo não discrimina". E uma guarda do
driver nasceu inoperante: o `make clean` exigia `DIM`, falhava com aviso, e a
limpeza por dimensão não acontecia.

---

## 5. Três achados do primeiro contato com o núcleo herdado

Nenhum previsto no projeto.

1. **`Point` é array cru.** `typedef coordtype Point[DIM]` (`coord.h:31`) não é
   copiável nem construível in-place, então `std::vector<Point>` não compila — e o
   erro sai dentro do `<bits/stl_construct.h>`, sem dizer isso. A ponte é
   `std::array<double,2>`, de layout garantido: conversão de tipo, não de dados.
   Foi a **única** fricção de tipo.
2. **O núcleo herdado tem estado global de arranque.**
   `_try_initialize_petsc` (`utils.c:34`) desreferencia um `static char **_argv`
   que só `higtree_initialize` preenche, e a flag `is_petsc_initialized` é
   `static` — inicializar o PETSc por fora não a desarma. Custou um SIGSEGV sem
   nenhuma linha apontando para o nosso código. Virou `t8flow::Arranque`.
3. **`solver.c` acopla os motores em tempo de link.** Referencia
   `_slv_create_hypre` e `_slv_create_debug_write` sem guarda de macro, enquanto
   SOR e ViennaCL ficam atrás de `#ifdef`. Não dá para levar só o PETSc.

---

## 6. O que o verde ainda não prova

Escrito desde o primeiro dia, em `t8flow/README.md`, para a suíte não ser lida
como garantia:

- nenhuma física roda; np vai só até 1; não há partição, franja, refino local,
  domínio nem malha multiblocos;
- o motor do solver **instalado** não é consultável — não existe
  `slv_get_engine`. A cláusula N5 cobra a tabela de nomes, não a instalação;
- um defeito que **preserve** a reprodução polinomial é invisível para esse
  oráculo;
- o lado 3D não existe.

---

## 7. Próximos passos, em ordem de retorno

1. **A malha multiblocos 2D** — uma árvore por bloco, com nível por bloco. É o
   multiblocos de verdade (resolução por região), e o `cmesh` multiárvore do
   t8code o dá nativamente. Cláusulas a acrescentar: contagem de árvores,
   folhas por bloco com níveis diferentes, e centroide dentro da caixa do seu
   bloco.
2. **A numeração global** contígua por rank, alcançada por **busca por ponto** —
   comparar com o percurso que a gerou é tautologia, e está medido que a suíte
   inteira fica verde com `n_visitadas++` no lugar de `mp_lookup`.
3. **Medir escalabilidade em 3D** antes de comprometer o desenho. O 2D rendeu
   1,62× em np=2 e 1,68× em np=4 por tamanho do caso (9377 células), não por
   defeito; o domínio 3D tem 375 mil. Sem essa medida, o escalonamento é
   promessa.
4. **A franja por caixas**, na profundidade de produção (5) e em np=3 — que é o
   menor np em que pareamento transposto aparece: com uma caixa por rank, 1×1
   esconde o erro.
5. **A cláusula F3**, "o suporte cai na grade e é local", que é o que a
   reprodução polinomial não vê.

Em aberto, e são decisões do usuário: o `cel_min` da banda fina, e se o
fechamento de contorno em DIM-1 também sobe para o reuso.
