# t8Flow — relatório de avanços

Estado em **11/10/2026**. Branch `t8flow`, worktree
`/home/castelo/HiGFlow-System/t8flow-wt`. Dois commits: `ea83f27` (esqueleto) e
`11239e0` (solver e discretizador 2D).

    suite do lado 2d    52 casos, 0 falharam
    contrato            24/24 clausulas (N1..N6, B1, F1, F2, M1..M5, R1..R5, G1..G5)
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
    t8flow/src/comum/t8flow-contract.h   24 clausulas
    t8flow/src/2d/mls2d.{hxx,cxx}         o discretizador
    t8flow/src/2d/solver2d.{hxx,cxx}      o solver
    t8flow/src/2d/malha2d.{hxx,cxx}       malha, refino local, numeracao global
    t8flow/tests/testing.h                arcabouco, zero dependencia externa
    t8flow/tests/2d/test-mls2d.cxx        5 casos
    t8flow/tests/2d/test-solver2d.cxx     6 casos
    t8flow/tests/2d/test-floresta2d.cxx   4 casos
    t8flow/tests/2d/test-malha2d.cxx     14 casos
    t8flow/tests/2d/test-refino2d.cxx    14 casos
    t8flow/tests/2d/test-numeracao2d.cxx  9 casos
    t8flow/bench/escalabilidade.cxx       ferramenta, fora da suite
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

**A malha multiblocos** é uma árvore por bloco, com a caixa e o nível de cada
bloco — resolução por região. As interfaces coincidentes são **unidas**: sem o
join, cada bloco é uma ilha e toda face na borda de um bloco aparece como "face
sem vizinho", que é o critério de contorno — dois blocos lado a lado do mesmo
domínio dariam contorno físico no meio do escoamento, em silêncio. Sobreposição
parcial de faces é recusada, porque o join do t8code não a representa.

**O refino local** refina até o nível pedido toda folha cujo centro caia numa
caixa, e **reequilibra**: depois dele a razão entre níveis de células vizinhas é no
máximo 2:1. Sem o `t8_forest_set_balance` a malha commita, as contagens fecham, e
o salto atravessa o estêncil do discretizador sem aviso — por isso o caso de
refino pede **dois** níveis de uma vez, que é o que faz o equilíbrio importar.

**A numeração global** é o que o solver pede da malha — um id por incógnita,
contíguo por processo —, e a curva de preenchimento do t8code já a entrega: as
folhas de um rank ocupam um intervalo na ordem da curva. O trabalho aqui foi
**verificá-la sem tautologia**: no HiGFlow, comparar a numeração com o percurso
que a gerou deixa a suíte inteira verde mesmo trocando o mapeador por um contador.
As três vias independentes usadas: aritmética de Morton (calculada no teste, sem
olhar a malha), bijeção sobre o intervalo local, e "escreve por percurso, lê por
geometria" com `id_do_ponto`.

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
| `t8_cmesh_set_join` removido | M3 **reprova no conteúdo** (0 com vizinho, 4 sem) — e a contagem de junções continua certa |
| vértices do QUAD v1 ↔ v2 | M2 reprova **só** na ordem de percurso; o centroide-na-caixa continua verde |
| nível uniforme global em vez de por árvore | M1 **reprova** em dois blocos: 64 folhas onde se esperavam 16 e 4 |
| `t8_forest_set_balance` removido | R3 **reprova** com salto de 2 níveis — e só ela |
| adapt não recursivo | R1 **reprova**: nível máximo 3 em vez de 4 |
| `id_global` sem o offset da árvore | G2 **reprova** (4 ids repetidos) e o processo morre na guarda do solver |
| `id_global` deslocado de 1 | G2 **reprova**: 1 id fora de `[0, 20)` |
| `id_do_ponto` com intervalo fechado | G4 **reprova**: o desempate passa a depender da ordem de visita |

A segunda linha é a mais informativa: ela mede que a reprodução polinomial
verifica a **discretização** e não a malha — por isso "o suporte cai na grade e é
local" precisa ser cláusula própria (F3), não subproduto.

A última justifica as fachadas existirem: elas convertem um abort três camadas
abaixo em exceção que nomeia a linha que falta. E foi a guarda do driver que
exige a linha de resumo que pegou o processo morto, em vez de reportar "nenhuma
falha".

**Dois testes meus estavam errados antes de o código estar, na rodada do
refino.** "Toda folha de nível máximo tem o centro dentro da caixa" é **falso** — o
critério julga a célula *antes* de dividir e os filhos herdam a decisão do pai,
então a região refinada é a caixa dilatada por até uma célula do nível final. E um
`t8_forest_unref` a mais deixou a malha **vazia**: o `commit` consome a floresta de
origem, liberar de novo não dá erro, as contagens vão a zero. Quem pegou foi a
conservação de área (0 contra 1) — a contagem de folhas não pegaria.

**Um oráculo meu estava frouxo, e só a sabotagem mostrou.** O comentário do
`test-malha2d` afirmava que "centroide dentro da caixa" pegava a troca da ordem
dos vértices do QUAD. Não pega — a troca espelha o mapeamento *dentro* da mesma
caixa, e todo centroide continua dentro dela. Fechei a lacuna com um caso que
ancora a orientação (a curva varre x primeiro, logo a folha 1 está deslocada em X
da folha 0), e ele reprova nos três blocos com a assinatura certa.

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

## 7. A escalabilidade, medida

A versão anterior deste relatório dizia "sem essa medida, o escalonamento é
promessa". Medido com `t8flow/bench/escalabilidade.cxx` — t8code puro,
`COMM_WORLD`, melhor de cinco execuções, restrito aos cores livres da máquina.

**3D**, HEX de nível 7 (2,1 M células; 3,0 M após o refino):

| etapa | np=2 | np=4 | desequilíbrio |
|---|---|---|---|
| construção | 1,89× (95%) | 1,95× (49%) | 1,00× |
| **percurso** (trabalho local) | **1,93× (97%)** | **3,70× (93%)** | 1,00× |
| refino + equilíbrio | 1,65× (83%) | 2,27× (57%) | **1,89× / 2,77×** |
| repartição | 1,64× (82%) | 2,17× (54%) | volta a 1,00× |
| percurso na malha adaptada | 1,90× (95%) | 3,66× (92%) | 1,00× |

**2D**, QUAD de nível 11 (4,2 M células; 6,2 M após o refino):

| etapa | np=2 | np=4 | desequilíbrio |
|---|---|---|---|
| **percurso** | **2,04× (102%)** | **3,82× (96%)** | 1,00× |
| refino + equilíbrio | 1,71× (86%) | 2,53× (63%) | **1,96× / 2,92×** |
| repartição | 1,33× (66%) | 2,02× (51%) | volta a 1,00× |

### 7.1 As três conclusões

**O trabalho local escala quase idealmente — 93% a 97%.** É o teto do que a
divisão pode dar, e ele está lá.

**A dimensão não era o fator; o tamanho era.** O 2D com 4,2 M células escala igual
ao 3D com 2,1 M. Logo o 1,62× do HiGFlow em np=2 não é limitação de 2D: eram 9377
células, ~2300 por processo, e a comunicação comendo o ganho. A explicação que
estava registrada como suposição agora está medida.

**O refino local desequilibra, e isso exige repartição** — o achado acionável. A
caixa de refino cai num pedaço da curva, e a partição anterior não sabe disso: em
np=4 o rank mais carregado ficou com **2,8× mais células**, nos dois DIM. Como o
relógio é o do rank mais lento, isso sozinho explica a queda (93% no percurso
contra 57% no refino). Reparticionar devolve o equilíbrio **exato** — 756 427
contra 756 428 folhas — e custa **3% a 6%** do tempo do refino.

### 7.2 O que a medida não diz

- **Não há física, solver, franja nem estêncil.** É um **teto**, não uma previsão:
  o que entra depois é justamente o que comunica mais.
- **np vai só até 4.** A máquina tem **8 cores físicos** (16 threads — `nproc`
  engana), e duas corridas de outra sessão estavam presas às CPUs 1 e 2.
- **Nível 6 e abaixo não é interpretável.** Com 262 k células os tempos ficam em
  10–140 ms e o cache domina — a medida deu ganho *superlinear* de 10× em np=2.
  Vale para as linhas acima de 100% nas tabelas.

### 7.3 Duas lições de medição

**Sem excluir as CPUs das corridas alheias, a medida deu 18× em np=8** —
impossível para trabalho local puro, e foi só por ser absurdo que duvidei. O
script passou a imprimir o estado da máquina junto com os números.

**Com três repetições um outlier passou:** o percurso adaptado em np=4 apareceu
como 1,05× numa rodada e 3,66× com cinco repetições. Eu teria reportado uma perda
que não existe.

---

## 8. Próximos passos, em ordem de retorno

1. **A franja por caixas**, na profundidade de produção (5) e em np=3 — que é o
   menor np em que pareamento transposto aparece: com uma caixa por rank, 1×1
   esconde o erro.
2. **A cláusula F3**, "o suporte cai na grade e é local", que é o que a
   reprodução polinomial não vê.

3. **A repartição depois de cada adaptação**, que a seção 7 transformou de
   suposição em recomendação com número.
4. **A guarda contra árvore de uma célula** — que não tem interior, toda faceta
   dela é de fronteira, e a franja declarada é de cinco. Entra junto com a franja.

Em aberto, e são decisões do usuário: o `cel_min` da banda fina, e se o
fechamento de contorno em DIM-1 também sobe para o reuso.
