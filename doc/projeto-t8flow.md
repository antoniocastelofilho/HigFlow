# t8Flow — o escalonamento do t8code, o multiblocos do HiGFlow

Batizado pelo usuário em **11/10/2026**. Sistema novo, construído sobre o
t8code, com duas exigências postas desde o início:

- as características de **escalonamento** do t8code — a floresta sendo a
  estrutura, a partição saindo da curva de preenchimento, nenhum processo
  segurando a malha inteira;
- **multiblocos** como o HiGFlow — e essa é a peça que casa melhor, porque o
  `cmesh` do t8code é nativamente multiárvore.

O levantamento de viabilidade está em [`projeto-sistema-t8code.md`](projeto-sistema-t8code.md)
(commit `7b7c769`). Este documento trata do desenho.

---

## 1. Decisão: o t8Flow herda o núcleo independente de malha da HiGTree

Pedido do usuário, em letras: *"um sistema que aproveita o que é independente de
malha do HiGTree, tal como Solver e Aproximações por MLS"*.

É a decisão certa e barata, e o motivo é que esses dois módulos **já não sabem o
que é uma malha** — não é refatoração a fazer, é propriedade a constatar.

### 1.1 O solver

`higtree/src/solver.h` é uma vtable sobre o sistema linear. A abstração fala em
índices:

```c
solver* slv_create(SolverEngine engine, int first_gid, size_t local_size);
void (*set_Ai)(struct solver *cs, int i, int numjs, const int *j, const real *v);
```

O próprio cabeçalho declara o contrato: *"everything above this file assembles
rows and asks for a solution without knowing which library solves it."* E a
recíproca também vale — nada abaixo dele sabe de onde vêm as linhas. Medido:
`solver.c`, `solver-petsc.c`, `solver-sor.c`, `solver-centralized.c` e
`solver-debug-write.c` têm **zero** referências a `hig_cell`, `higcit`,
`sim_domain` ou `psd_`. O `#include "higtree.h"` em `solver.h` é vestigial.

Backends disponíveis: PETSc, hypre, SOR, ViennaCL, centralized, hollow,
debug-write, escolhidos em tempo de execução por `HIGTREE_SOLVER_ENGINE`.
Advertência que o cabeçalho registra e que o t8Flow herda junto: as duas listas
(`engine_names[]` e o enum `SolverEngine`) estão **deslocadas de uma posição**, e
nada confere isso em tempo de compilação — errar roda silenciosamente um motor
diferente do nome pedido.

### 1.2 As aproximações por mínimos quadrados móveis

`higtree/src/wls.h` depende apenas de `coord.h` e `types.h`. A API inteira é
geométrica:

```c
real wls_set_points_and_calc(wls_interpolator *wls, int numpts, Point pts[], Point x, real w[]);
```

Pontos de suporte e um ponto alvo entram; pesos saem. Isso é exatamente a
**divisão malha / discretização** que o contrato de malha de 18 cláusulas já
adota: quem reproduz polinômio é o ajuste, igual para qualquer malha; o que a
malha entrega é o **suporte** — quais células ficam perto.

A ponte malha→WLS existe hoje em `domain.c` (14 chamadas) e em
`hig-flow-io.c` (20, para interpolar na saída). É essa ponte que o t8Flow
reescreve; o ajuste, não.

### 1.3 O que mais vem de graça

Inventário **medido**, não estimado — cada unidade de tradução compilada com os
cabeçalhos de malha substituídos por arquivos vazios:

| módulo | papel |
|---|---|
| `wls.{c,h}` | mínimos quadrados móveis, ordem e tipo de peso configuráveis |
| `solver.{c,h}` + `solver-petsc.c`, `solver-sor.c`, `solver-centralized.c`, `solver-debug-write.c` | sistema linear atrás de vtable |
| `coord.{c,h}`, `rect.{c,h}`, `dim.h`, `types.h` | geometria de ponto e de caixa |
| `point-cloud.{c,h}`, `point-mapper.{c,h}` | nuvem de pontos e mapa ponto→índice |
| `utils.{c,h}`, `allocator.{c,h}`, `uniqueid.{c,h}`, `rng.cpp`, `Debug-c.h` | infraestrutura |
| `isend-pool.{c,h}`, `term-det.cpp` | envio assíncrono e detecção de terminação em MPI |

São **4822 das 21581 linhas** de `higtree/src` — cerca de 22%, e é a fração que
não tem nada a reescrever.

Do lado acoplado, pelo mesmo teste: `higtree.c`, `higtree-iterator.c`,
`domain.c`, `pdomain.c`, `lbal.c`, `build-fringe.cpp`, `mapper.c`,
`mapper-syncer.c`, `global-stencil.c`, `higtree-io.c`, `higtree-parallel.c`,
`higtree-serialize.c`, `hig-mesh-snapshot.c`. Esses são a estrutura — é deles que
o t8Flow se separa.

> **O teste que mediu isso tinha um modo de falha, e ele disparou na primeira
> tentativa.** Compilar o fonte *no próprio diretório* com `-I<sombra>` não
> substitui nada: `#include "x.h"` procura primeiro no diretório do arquivo
> fonte, antes de qualquer `-I`. A primeira rodada declarou `domain.c` e
> `higtree.c` independentes de malha — resultado absurdo que denunciou o
> defeito. Corrigido copiando cada fonte para fora antes de compilar, o veredito
> passou a discriminar. Vale a regra geral deste projeto: **um oráculo que não
> reprova nada ainda não foi verificado.**

---

## 2. Decisão: o particionador e o domínio são do t8code

Pedido do usuário: *"utilizar o particionador e domínio do t8code"*. É a
**inversão de posse** — a floresta deixa de ser produtora de malha e passa a ser
a estrutura. É também a única das duas alternativas que entrega o escalonamento,
e já está prototipada e verificada neste repositório.

### 2.1 O que o t8Flow deixa de herdar

Com a partição e o domínio vindo do t8code, saem de cena os três maiores
arquivos da camada de estrutura:

- `lbal.c` (57 KB) — o balanceador. Sem `lb_calc_partition`: cada rank fica com
  o que a curva de preenchimento lhe deu.
- `build-fringe.cpp` — a construção da franja por pareamento negociado.
- a metade de `pdomain.c` (54 KB) que calcula partição e troca blocos de franja.

O que **não** sai é o contrato que esses arquivos cumprem. Ele continua valendo,
e está lido do código no cabeçalho de `higtree/tests/t8code/t8-particao-grafo.h`:

> quem envia — para cada `to_send[i]`, enumera os elementos de `local_tree` no
> intervalo `[lo_idx, hi_idx)` e manda, concatenados na ordem de `i`;
> quem recebe — para cada `to_recv_trees[k]`, enumera a árvore **inteira** e
> consome o buffer na ordem de `k`.

A regra que faz isso sair igual nos dois lados sem negociação: **a caixa dona do
dado no laço de fora, a que precisa dele no de dentro**, com um `allgather` das
caixas e nada acertado por mensagem. Errar isso não quebra nada visível — entrega
à célula de franja o valor de outra posição. Medido com os laços transpostos:
np=2 **passa** (uma caixa por rank, 1×1 esconde o erro) e np=3 **reprova**, com
228 de 304 células de franja contendo valor de outra posição.

### 2.2 O que já existe, pronto

`t8_monta_dominio_particionado(lo, hi, nb, sd, pg)` em
`higtree/tests/t8code/t8-particao-grafo.cxx`, exposta no exemplo como
`HIGFLOW_MALHA=t8code-particao`, verificada em np = 1, 2 e 3 nos dois DIM pelo
`test-particao-t8code`. Ela falha **abortando**, de propósito: domínio com
vizinhança errada dá resultado plausível.

O que a tornou possível: a região de um rank na curva de preenchimento é quase um
retângulo, e o grafo representa o que se envia como faixa retangular. Medido no
brick 160×40 — caixas por rank: np=2 dá 1, np=3 dá 2 a 3, np=4 e np=8 dão 1.

### 2.3 O achado que esta decisão traz: a franja não pode ser o ghost do t8code

Medido nos dois lados, e é o ponto de projeto mais importante desta seção:

| | profundidade | vizinhança |
|---|---|---|
| ghost do t8code | **1 camada** | só `T8_GHOST_FACES` — o cabeçalho diz *"currently only T8_GHOST_FACES is supported"* |
| franja do HiGFlow | **5 células** (`pg_set_fringe_size(pg, 5)` em `malha-adaptativa.c:684` e `hig-flow-kernel.c`) | o padrão da HiGTree é 1; o 5 é escolha do HiGFlow |

A camada nativa do t8code **não cobre** o suporte que o MLS pede. Isso não é
limitação do t8code: ghost de uma camada é o que um método de volumes finitos
clássico precisa, e o suporte largo é consequência da aproximação por mínimos
quadrados — a mesma que o t8Flow herda por decisão.

Portanto a franja do t8Flow é construída **por caixas**, com profundidade
pedida, e não pelo ghost nativo. É exatamente o que o protótipo já faz: ele lê
`pg_get_fringe_size(pg)` e monta as faixas. O ghost nativo continua útil para
outra coisa — distinguir face de contorno de face de franja, que sem ele seria
indistinguível (`t8-partition.cxx:42`).

**Consequência para a suíte:** o caso que afirma a franja tem de rodar com a
profundidade de produção, não com 1. Uma franja de profundidade 1 verificada e
uma de 5 em produção é a forma clássica de passar verde e divergir na simulação.

### 2.4 Armadilhas desta frente, cada uma medida

- **`set_for_coarsening=1` não é recursivo.** Ele impede que uma família seja
  dividida entre ranks, mas preserva a família do nível **mais fino**; uma célula
  base refinada duas vezes contém família de famílias. Medido em 3D com np=2:
  oito células base saíram divididas, a materialização criou as que faltavam, e o
  domínio ficou com 541 folhas contra as 533 da floresta. A decomposição desce
  até o nível em que a **posse fecha**.
- **`set_adapt`, `set_ghost` e `set_partition` no mesmo commit do t8code**
  derrubam com np>1 e **funcionam** em np=1. Refino num commit, partição com
  ghost em outro.
- **O `cmesh` em `COMM_WORLD` sai distribuído**: o rank 0, produzindo em série,
  via só a parte dele — 47 folhas em np=2 contra 79. Para produção serial, cmesh
  **e** floresta em `COMM_SELF`.
- **A varredura por ponto alcança também as árvores de franja**, então dois ranks
  reivindicam o mesmo ponto. O filtro é `lid < n_local`.
- **O índice que `t8_forest_leaf_face_neighbors` devolve é da floresta local, não
  da árvore.** Tratá-lo como índice de árvore entrega um elemento existente, só
  que **outro** — erro silencioso. Use `t8_forest_get_leaf_element`.
- **`t8_forest_element_points_inside` devolve TRUE para dois vizinhos** quando o
  ponto está na face comum. A primitiva não desempata; a convenção é imposta
  sobre ela, e o caso tem de inverter o percurso para não passar por sorte.

---

## 3. O que o t8Flow tem de fornecer ao núcleo herdado

Esta é a consequência de projeto mais útil do inventário, porque é curta. O
núcleo independente precisa de **duas** coisas da malha, e só:

1. **Uma numeração global das incógnitas, contígua por processo** — é o que
   `slv_create(engine, first_gid, local_size)` pede, e nada além disso.
2. **O suporte de cada ponto** — a lista de centroides vizinhos que vai ao
   `wls_set_points_and_calc`.

As duas o t8code dá, e a primeira ele dá *melhor* do que o caminho atual: a
curva de preenchimento já numera as folhas de um rank contiguamente, que é
precisamente a forma que o solver quer. Hoje essa numeração vem do mapeador
(`mp_assign_from_celliterator`) sobre o iterador de células.

Fato já medido e que vale para o desenho: **a identidade entre o id do mapeador e
o contador do percurso é verdadeira por construção, e nenhum teste da suíte a
distingue** — trocar `mp_lookup` por `n_visitadas++` deixa a suíte inteira verde.
Consequência direta para o t8Flow: o oráculo que autorizar a numeração nova
**não** pode compará-la com o percurso que a gerou; tem de alcançar a célula por
outra via (busca por ponto).

---

## 4. A suíte nasce com o sistema

Pedido do usuário: *"gerar também uma suíte de testes assim como foi feito na
HiGTree inicialmente"*. Na HiGTree a suíte foi escrita **depois**, e o custo
disso está medido: todos os defeitos de `domain.c` corrigidos entre 15 e
18/09/2026 viviam na biblioteca e foram achados por uma simulação 3D abortando,
não por teste. Um deles tinha **cinco anos** e caiu na primeira execução da
suíte nova. No t8Flow a suíte entra junto com o código.

### 4.1 O que replicar do desenho da HiGTree

Cinco peças, cada uma vinda de uma armadilha real:

1. **Afirmar VALOR, não forma.** A suíte ATF de 2020 (71 casos) verificava
   `numelems >= 3` e `rhs != 0`; os nove defeitos de fechamento daquela semana
   passariam por todas elas. Nunca aceitar "retornou algo plausível".
2. **Contrato lido do próprio cabeçalho.** `higtree/src/hig-mesh-contract.h` tem
   18 cláusulas, e `ci/run_higtree_tests.py` **lê as cláusulas do arquivo** e
   exige caso passando para cada uma: cláusula sem teste sai como *falha*, não
   como omissão. É o que impede a suíte de parecer completa por silêncio.
3. **As três guardas do driver.** (a) exigir a variável de ambiente que o build
   precisa e *dizer* qual falta, em vez de reportar 0 de N; (b) `make clean` por
   dimensão — os objetos não carregam o `DIM` no nome, e um build 2D sobre
   árvore 3D liga objetos 3D e falha como segfault, não como erro de build;
   (c) exigir a **linha de resumo** do binário: sem ela, "nenhuma falha" e
   "nenhum caso rodou" são indistinguíveis — que é exatamente como se
   apresentam binário obsoleto e segfault no meio.
4. **O alvo de teste depende da biblioteca.** Com `%-2d: %.c testing.h` o make
   considera o binário atualizado quando só a biblioteca mudou; o teste roda
   obsoleto e dá os números da versão anterior. Aconteceu no dia em que a suíte
   foi escrita e quase enterrou uma correção como "não funciona".
5. **Validar nos dois sentidos.** Verde que nunca foi vermelho não prova nada, e
   o padrão de falha tem de bater com o mecanismo. Quando a seleção de parede
   entrou, o controle foi: com o critério antigo restaurado, `test-stencil-selection`
   falha e `test-partition-independence` **continua passando** — a segunda linha
   vale tanto quanto a primeira, porque se as duas falhassem a primeira estaria
   medindo seleção por acidente.

Arcabouço: um `testing.h` sem dependência externa (a ATF morreu por depender do
`atf-c`, que não está instalado em lugar nenhum), binário que imprime
`caso <nome> PASS|FAIL <números>` e uma linha `resumo <n> casos, <k> falharam`.
Teste que chama `MPI_Init` vai **sempre** por `mpirun`, inclusive em np=1 —
rodado direto ele trava sem imprimir e aparece como timeout, não como falha.

### 4.2 Os oráculos, e os que não discriminam

Este é o capital mais valioso a transportar, porque cada linha abaixo foi
descoberta descobrindo que um teste verde não provava nada:

- reprodução polinomial **não** detecta franja ausente (campo linear sai exato
  com suporte de um lado só, por extrapolação);
- reprodução polinomial **não** detecta parede de fechamento errada com campo
  linear (Lagrange em qualquer eixo acerta um linear) — daí o segundo campo, de
  grau *acima* da ordem;
- reprodução polinomial **não** detecta suporte deslocado: deslocando todos os
  centroides meia célula ela continua exata, porque o campo é avaliado nos
  próprios pontos. Ela verifica a **discretização**, não a malha. Quem verifica
  a malha é *"o suporte cai na grade e é local"*;
- contagem de células e volume **não** detectam pareamento de franja trocado —
  quem detecta é o **conteúdo** da franja depois do sync, com coeficientes
  distintos por direção para que o valor identifique a posição. Medido: a versão
  transposta passa em np=2 (uma caixa por rank, 1×1 esconde) e reprova em np=3
  com 228 de 304 células de franja contendo valor de outra posição;
- comparar a numeração contra o percurso que a gerou é **tautologia**: trocar
  `mp_lookup` por `n_visitadas++` deixa a suíte inteira verde, nos dois DIM e nos
  três np. A numeração do t8Flow tem de ser alcançada por busca por ponto;
- desempate por ordem de percurso passa **por sorte** quando a curva de
  preenchimento visita a célula de menor coordenada primeiro. O caso tem de
  inverter o percurso e exigir a mesma resposta.

E o truque que dispensa gancho de depuração, que vale para o t8Flow inteiro:
**valores constantes e diferentes por parede** (ou coeficientes distintos por
direção) fazem o valor medido identificar sozinho a escolha interna, e o teste
fica sobre comportamento em vez de implementação.

### 4.3 A suíte inicial que proponho

Opinião, a revisar com o usuário. A ordem é de dependência, e cada linha traz o
oráculo — um teste sem oráculo que discrimine não entra:

| teste | afirma | oráculo |
|---|---|---|
| `test-nucleo-herdado` | solver e WLS funcionam fora da HiGTree | sistema linear de solução analítica; reprodução polinomial com campo de grau acima da ordem |
| `test-motor-do-solver` | o motor que roda é o que o nome pede | as duas listas de `solver.h` estão deslocadas de uma posição e nada confere em build — o caso confere em execução |
| `test-numeracao-global` | ids contíguos por rank, sem furo nem repetição | alcançar a célula por **busca por ponto**, nunca pelo percurso que numerou |
| `test-suporte-local` | o suporte cai na grade e é local | geometria, não valor interpolado — é o que pega suporte deslocado |
| `test-multiblocos` | interface entre blocos do mesmo domínio **não** é contorno | face sem vizinho contra caixa: a decisão D3 já tomada na HiGTree |
| `test-franja-conteudo` | o pareamento envia/recebe na mesma ordem | conteúdo após sync, coeficientes distintos por direção, np = 1, 2 **e 3**, na profundidade de produção (5) e não em 1 |
| `test-refino-2a1` | razão no máximo 2:1, nenhuma árvore de uma célula | geométrico; a árvore de uma célula não tem interior e a franja é de cinco |
| `test-posse-distribuida` | nenhum rank segura a malha inteira | contagem local contra global — é a propriedade que motivou o t8code |
| `test-particao-sem-lbal` | a partição final é a da curva de preenchimento | o domínio se monta sem `lb_calc_partition`, e o conteúdo da franja confere — o caso tem de rodar em np=3, porque np=2 esconde o pareamento transposto |

Com um contrato próprio (`t8flow-contract.h` ou nome equivalente) no mesmo
padrão: cláusulas numeradas no cabeçalho, driver lendo dali, cláusula sem teste
contando como falha.

### 4.4 O que o verde não vai provar

A HiGTree mantém essa seção escrita em `higtree/tests/README.md`, e ela é a parte
mais útil do documento. No t8Flow, de saída: nenhuma física roda na suíte, np vai
até 3 e não até dezenas, nada chega perto do teto de árvores por domínio, e um
defeito que **preserve** a reprodução polinomial é invisível para esse oráculo.
Escrever isso desde o primeiro dia evita que a suíte seja lida como garantia.

---

## 5. Decisão: um domínio só

**Decidido pelo usuário em 11/10/2026.** Dois domínios paralelos — um adaptado e
um uniforme no nível mais fino — não escalam, e há três razões independentes,
cada uma suficiente por si. Elas estão em 5.1 a 5.3; o que a decisão compra e o
que ela obriga está em 5.5 e 5.6.

### 5.1 O domínio uniforme é a malha cheia, por construção

`higflow_create_amr_info_mult` não aproxima nada: copia os limites do domínio
inteiro, põe `numlevels = 1`, `numpatches = 1` e usa o `delta` do **último**
nível —

```c
mi_mult->levels[0].delta[dim] = mi->levels[mi->numlevels - 1].delta[dim];
int nc = (int)round((mi_mult->h[dim] - mi_mult->l[dim])/mi_mult->levels[0].delta[dim]);
```

Logo o segundo domínio tem exatamente o número de células que a malha teria se
**tudo** fosse refinado ao máximo. Ele não encolhe com o AMR: medido no caso 3D,
`dpfracvol` em 375 000 células (50×50×150) contra ~81 000 da árvore adaptada —
4,6× com **um** nível, e em 3D a razão cresce 8× por nível adicional.

A consequência é a que mata a ideia: o escalonamento do t8code se aplicaria à
parte que **não** domina o custo. O gargalo migra para o domínio uniforme, que é
O(volume) por definição, e o AMR deixa de comprar o que foi comprado para
comprar.

### 5.2 A partição conjunta de hoje depende justamente do `lbal`, que sai

O que faz os dois domínios casarem hoje é um recurso do balanceador: **grupos**.

```c
load_balancer *lb = lb_create(MPI_COMM_WORLD, 2);
for(...) lb_add_input_tree(lb, higio_read_from_amr_info(mi[i]),      true, 0);  // adaptado
for(...) lb_add_input_tree(lb, higio_read_from_amr_info(mi_mult[i]), true, 1);  // uniforme
lb_calc_partition(lb, pg);
```

Uma única `lb_calc_partition` reparte os dois grupos, e os dois domínios recebem
o **mesmo** `partition_graph` — é daí que sai a coerência entre as fronteiras de
rank das duas malhas.

Com a decisão 2 (partição do t8code) o `lbal` sai de cena, e esse recurso não tem
equivalente: a curva de preenchimento de **uma** floresta não reparte duas
resoluções diferentes de forma consistente, e duas florestas dão duas curvas com
fronteiras de rank distintas — a célula de fração que um rank precisa passa a
estar em outro rank. Seria preciso inventar um terceiro mecanismo de troca, entre
duas partições independentes. É exatamente o tipo de camada onde as falhas
silenciosas deste projeto moram.

### 5.3 O caminho multifásico não consulta a fonte de malha hoje

Está declarado no próprio código (`hig-flow-kernel.c`, por volta de 1795):
`higflow_partition_domain_multiphase` lê os `.amr` **sempre**, e quem pedir
`HIGFLOW_MALHA=t8code` num caso multifásico tem a escolha ignorada — o código
apenas **denuncia** a discrepância em stdout. Ou seja, a alternativa de dois
domínios não está só mais lenta: ela está, hoje, fora do alcance do t8code por
construção.

### 5.4 A ressalva, e ela é de formulação

"Um domínio só" não é pura engenharia. Se a formulação do VOF precisa de
vizinhança de tamanho uniforme para reconstruir a interface, a resposta escalável
**não** é voltar ao domínio uniforme global — é **manter a interface inteira no
nível mais fino com adaptação dinâmica**, que dá a mesma vizinhança uniforme onde
ela importa a custo O(área) em vez de O(volume). É o que o t8code faz
nativamente (adapt + 2:1 + re-partição), e já existe no projeto em outra frente
(ver `doc/projeto-amr-dinamico.md`).

Essa escolha — vizinhança uniforme por domínio separado, ou por nível uniforme na
interface — é de formulação, e é do usuário.

Com a decisão tomada por um domínio só, o caminho é o segundo: **a vizinhança
uniforme é obtida mantendo a interface inteira no nível mais fino**, não por um
domínio separado. É o que a seção 5.6 transforma em invariante.

### 5.5 O que a decisão compra

1. **A classe de defeito desaparece estruturalmente, não por cuidado.** Não há
   dois índices para confundir. A dimensão do que sai: nove propriedades criadas
   sobre `psdmult` (`dpfracvol`, `dpfracvolaux`, `dpvisc`, `dpdens`,
   `dpcurvature`, `dpdistance`, `dpIF[]`, `dpnormal[]`), 204 sítios em 22
   arquivos, e — o pior — **dois estênceis**: `hig-flow-step-multiphase.c` usa
   `ed.mult.stn` 114 vezes contra 1 de `ed.stn`, e o caminho 3D fazia o inverso.
   Foi essa divergência que sobreviveu anos.
2. **A decisão 2 fecha sozinha.** Uma floresta, uma partição, uma numeração, uma
   camada de franja. Nada a casar entre duas curvas de preenchimento, e nenhum
   terceiro mecanismo de troca a inventar.
3. **O custo passa a ser O(células adaptadas).** Deixa de existir o vetor de 375
   mil posições que não encolhe com o refino.

### 5.6 O que a decisão obriga: a banda fina vira invariante de malha

Em troca, a malha passa a ter uma obrigação que antes era do segundo domínio:
**toda célula a menos de `cel_min` células da interface tem de estar no nível
mais fino.** No HiGFlow isso é convenção de exemplo; no t8Flow é cláusula de
contrato e caso de suíte.

A boa notícia é que a maquinaria existe e é reaproveitável, em
`higflow/examples-common/malha-adaptativa.{c,h}`:

- `malha_adapt_limiares(h_base, niveis, cel_min, thr)` — os limiares de distância
  **derivados** da regra das células mínimas, não digitados;
- `malha_adapt_mede_banda(ns, niveis, h_base)` — o **oráculo**: mede na malha de
  verdade a menor distância de uma semente de interface a uma célula que não está
  no nível mais fino, já descontada a incerteza da semente, e tem de dar
  `>= cel_min`;
- `malha_adapt_reune_sementes` — `Allgather` + `Allgatherv` das sementes, sem o
  que o critério sai refinado só em torno do pedaço de interface de quem o
  escreveu.

**E o próprio código registra a armadilha que a decisão elimina de brinde:** a
medida da banda feita no domínio errado é um oráculo que nunca pode falhar. Medir
no `sdmult` — uniforme por desenho — significa que não existe célula grossa
nenhuma, e a busca devolvia −1 **sempre**; o comentário em
`malha-adaptativa.c:589` diz que a distinção *"custou uma medida vazia"*. Com um
domínio só, não há onde cometer esse erro.

**Para a suíte:** um caso que afirme a regra das células mínimas, com controle no
sentido inverso — pedir `cel_min` acima do que a malha entrega tem de
**reprovar**. Oráculo que só foi visto verde não vale.

**O que a decisão não resolve, e segue para o AMR dinâmico:** a razão 2:1 na
borda da banda, e o custo de re-adaptar quando a interface se move. Previsto em
`doc/projeto-amr-dinamico.md`.

---

## 6. Decisão: C++ na norma do t8code

Decidido pelo usuário: o t8Flow é escrito em C++, compilado com o mesmo
compilador e a mesma norma do t8code. Isso não é preferência de estilo — é o que
remove a fronteira que hoje atravessa a integração.

### 6.1 Qual norma, e por que ela não é escolha nossa

O t8code **impõe** a norma ao consumidor, e isso está no alvo e não na
documentação:

```cmake
# bibliotecas/t8code/src/CMakeLists.txt:29
target_compile_features( T8 PUBLIC c_std_11 cxx_std_20 )
```

Sendo `PUBLIC`, a exigência reaparece na interface instalada —
`INTERFACE_COMPILE_FEATURES "c_std_11;cxx_std_20"` em
`install/cmake/T8CODE-targets.cmake:61` — e todo consumidor que linka pelo CMake
a herda. Com `CMAKE_CXX_EXTENSIONS` no padrão, o GCC recebe **`-std=gnu++20`**, e
é esse o alvo: `gnu++20`, não `c++20` estrito, porque o código herdado usa
extensões GNU.

**A fronteira que desaparece:** hoje a HiGTree constrói em `-std=gnu++17`, as
unidades do adaptador t8code compilam separadas em C++20, e só assinaturas
`extern "C"` atravessam entre as duas — um dos três custos de integração
registrados. No t8Flow há **uma** norma, e nada a atravessar.

### 6.2 Medido: o núcleo herdado já compila nessa norma

Não é promessa. As **16 unidades** do núcleo independente de malha compiladas com
`mpic++ -std=gnu++20`, sem uma única alteração no código:

| | gnu++17 | gnu++20 |
|---|---|---|
| `wls.c`, `solver.c`, `solver-petsc.c`, `solver-centralized.c`, `solver-sor.c`, `solver-debug-write.c` | ok | **ok** |
| `coord.c`, `rect.c`, `point-cloud.c`, `point-mapper.c`, `utils.c` | ok | **ok** |
| `allocator.c`, `uniqueid.c`, `isend-pool.c`, `rng.cpp`, `term-det.cpp` | ok | **ok** |

A decisão 1 (herdar solver e MLS) e a decisão 6 (norma do t8code) são portanto
compatíveis sem trabalho de porte. Isso era o risco óbvio desta combinação, e ele
não se materializou.

### 6.3 Recomendação: não levar o `-fpermissive`

A HiGTree compila C como C++ e precisa de `-fpermissive` para as conversões
implícitas `void*` → `T*`. Medido quanto custa abandonar isso no t8Flow:

    sem -fpermissive, em gnu++20:
      solver-sor.c     11 erros   (todos void* -> T*)
      os outros 15     limpos

Ou seja, **15 das 16 unidades já são C++ válido sem concessão nenhuma**. Recomendo
compilar o t8Flow sem `-fpermissive` desde o primeiro dia e resolver o `solver-sor.c`
com 11 casts explícitos — ou deixá-lo de fora, já que é backend opcional
(`USE_SOR`) e o PETSc é o motor de produção. Herdar o `-fpermissive` significaria
desligar permanentemente uma checagem do compilador para um arquivo que nem é
usado em produção.

*Observação de honestidade:* compilar o núcleo herdado em C++20 **não** o torna
C++ — ele continua C idiomático, com ponteiros crus e `malloc`. A camada nova do
t8Flow é que pode ser C++ de verdade (RAII nos recursos do t8code, que é
justamente onde vazamento de floresta e de cmesh dói). O núcleo entra como está,
e não há `extern "C"` a escrever: o mesmo compilador, a mesma norma, dos dois
lados.

---

## 7. Decisão: código separado para 2D e 3D, começando pelo 2D

Decidido pelo usuário, e a razão está no t8code: ele tem malhas 2D e 3D
**nativas**. `T8_ECLASS_QUAD` e `T8_ECLASS_HEX` são classes distintas, com 4 e 8
vértices, 4 e 6 faces, e esquemas próprios. O código novo do t8Flow é escrito
duas vezes, explicitamente, em vez de uma vez com `#if DIM == 3`.

    t8flow/src/comum/     sem dimensão: arranque, contrato
    t8flow/src/2d/        o lado 2D
    t8flow/src/3d/        o lado 3D, quando existir
    t8flow/tests/2d/

E o `make` passou a ter um alvo por lado — `make 2d`, e um `make 3d` que **diz**
que ainda não existe em vez de falhar obscuro.

**O que esta decisão NÃO muda:** o núcleo herdado continua parametrizado por
macro, porque `Point` **é** `real[DIM]` (`coord.h:31`) — a interface dele é
assim. As unidades herdadas são compiladas com `-DDIM=2` para o lado 2D, em
objetos separados. Isso não é custo novo: a HiGTree sempre construiu duas
bibliotecas. O que a separação compra é no código novo: nenhum `#if` no caminho
quente, tipos concretos (`Ponto2` é `std::array<double,2>`, copiável, que entra
em container da STL — coisa que `Point` não faz), e o compilador vendo
constantes.

### 7.1 A ordem: solver e discretizador primeiro

Também do usuário: *"adaptar primeiramente o Solver e o Discretizador com
Mínimos quadrados"*. É a ordem certa, porque esses dois são o que a decisão 1 já
garantiu existir — adaptá-los é fachada, não formulação, e o resultado é
testável sem malha nenhuma.

**`src/2d/mls2d.hxx` — o discretizador.** O herdado está certo e é independente
de malha; a fachada acrescenta o que falta nele: tipos 2D concretos, tempo de
vida (`wls_create`/`wls_destroy` → construtor/destrutor), e **a condição de uso
afirmada** — abaixo do mínimo de pontos para a ordem, `pesos()` **recusa** em vez
de aproximar. Numeros ruins em silêncio são o modo de falha que este projeto mais
paga.

**`src/2d/solver2d.hxx` — o solver.** A fachada concentra quatro armadilhas do
herdado: tempo de vida; **a ordem de montagem** (há só `monta()`, que faz
rhs→output→matriz, e não há como fazer a ordem errada que custou um SIGSEGV); **o
nome do motor** (as duas listas do `solver.h` estão deslocadas de uma posição e
nada confere isso em tempo de compilação — a aritmética fica em um lugar, cobrada
por teste); e **toda linha preenchida** (o caminho de pré-alocação do PETSc afirma
`col_ids[i] && vals[i]` para toda linha local, então linha faltando não é sistema
singular, é abort três camadas abaixo — `monta()` recusa antes, dizendo qual linha
falta).

**Lacuna do herdado que nenhum teste fecha daqui:** não existe `slv_get_engine`,
então qual motor foi **instalado** não é consultável, e com `Motor::Qualquer` a
escolha vem de variável de ambiente ou heurística. A fachada devolve o motor
**pedido**, e diz isso; quem precisa de garantia pede o motor explicitamente.

**O que a fachada não esconde, de propósito:** a numeração. `primeiro_id` e
`n_local` são a única coisa que o solver precisa da malha, e é exatamente o que a
curva de preenchimento dá naturalmente (seção 3).

---

## 8. O que continua a decidir

1. **Onde o t8Flow vive** — repositório próprio, ou dentro do HiGFlow-System
   reaproveitando a suíte.
2. **O `cel_min` da banda fina** — quantas células finas de cada lado da
   interface. Há valor derivado em `malha_adapt_limiares`, mas o número é de
   formulação.
3. **Até onde sobe o reuso.** Decidido para solver e MLS. Em aberto: fechamento
   de contorno em DIM-1 (projeção e interpolação, que também são de
   discretização e também independem da estrutura) e a camada de E/S.

---

## 9. O obstáculo real, já identificado

Com refino, a faixa que o remetente envia e a árvore de franja que o receptor
enumera precisam ter a **mesma subdivisão interna** — é trocar estrutura, não
índices. O `higtree/tests/test-particao-t8code.c` é onde isso vai aparecer, e o
oráculo dele é o certo: confere o **conteúdo** da franja depois do `dp_sync`, com
coeficientes distintos por direção, porque pareamento trocado monta o domínio,
roda o sync e entrega à célula de franja o valor de outra posição — sem erro
nenhum.

---

## 10. O estado do código, em 11/10/2026

Vive no worktree `t8flow`, branch `t8flow`
(`/home/castelo/HiGFlow-System/t8flow-wt`), com instruções em
`t8flow/README.md`. O primeiro commit (`ea83f27`) trouxe o esqueleto; o estado
atual tem as duas fachadas da seção 7.1.

    t8flow/Makefile                     alvos `2d` e `3d`, gnu++20, sem -fpermissive
    t8flow/src/comum/arranque.{hxx,cxx} a ordem de arranque, em um objeto
    t8flow/src/comum/t8flow-contract.h  9 clausulas: N1..N6 B1 F1 F2
    t8flow/src/2d/mls2d.{hxx,cxx}       o discretizador
    t8flow/src/2d/solver2d.{hxx,cxx}    o solver
    t8flow/tests/testing.h              arcabouco, zero dependencia externa
    t8flow/tests/2d/test-mls2d.cxx          5 casos
    t8flow/tests/2d/test-solver2d.cxx       6 casos
    t8flow/tests/2d/test-floresta2d.cxx     4 casos
    ci/run_t8flow_tests.py              driver, com as quatro guardas

    lado 2d: 15 casos, 0 falharam; contrato 9/9

### 10.1 Três achados do primeiro contato com o núcleo herdado

1. **`Point` é array cru.** `typedef coordtype Point[DIM]` (`coord.h:31`) não é
   copiável nem construível in-place, então `std::vector<Point>` **não compila** —
   e o erro sai dentro do `<bits/stl_construct.h>`, sem dizer isto. A ponte é
   `std::array<double,2>`, que a norma garante ter o mesmo layout: a travessia é
   conversão de **tipo**, não de dados, e um `static_assert` a sustenta. Foi a
   **única** fricção de tipo encontrada, e é por isso que `Ponto2` existe.
2. **O núcleo herdado tem estado global de arranque.**
   `_try_initialize_petsc` (`utils.c:34`) desreferencia um `static char **_argv`
   que só `higtree_initialize` preenche, e a flag `is_petsc_initialized` é
   `static` — inicializar o PETSc por fora não a desarma. Custou um SIGSEGV sem
   nenhuma linha apontando para o nosso código. Virou `t8flow::Arranque`, um
   objeto: a ordem fica no construtor, uma vez, e não em cada `main`.
3. **`solver.c` acopla os motores em tempo de link.** Referencia
   `_slv_create_hypre` e `_slv_create_debug_write` **sem** guarda de macro,
   enquanto o SOR e o ViennaCL ficam atrás de `#ifdef`. Não dá para levar só o
   PETSc: o link exige hypre, debug-write e centralized. Uma guarda por motor, no
   padrão que o SOR já usa, é limpeza barata e está na lista.

### 10.2 Os controles, porque verde que nunca foi vermelho não prova nada

    ordem do ajuste 2 -> 1            N1 REPROVA nas duas metades
    centroides deslocados meia h      F2 CONTINUA VERDE  (limitacao conhecida)
    clausula -> caso inexistente      driver: CASO AUSENTE, contrato 4/5, saida 1
    guarda de suporte minimo removida N3 REPROVA: "suporte de 5 pontos
                                        (minimo 6) passou sem reclamar"
    guarda de linha faltando removida o processo MORRE no abort do herdado, sem
                                        linha de resumo -- que e' o que a guarda 3
                                        do driver pega

As duas últimas linhas medem o valor das fachadas: elas convertem um abort três
camadas abaixo em exceção que nomeia a linha que falta. A segunda linha mede o
contrário — que a reprodução polinomial verifica a discretização e **não** a
malha, e é por isso que "o suporte cai na grade e é local" precisa ser cláusula
própria (F3), não subproduto.

Dois oráculos que passavam por acidente foram pegos pelo próprio contraste: com o
ponto de consulta **no centro de simetria** da nuvem, o ajuste de ordem 2
reproduziu um campo cúbico com erro 1,1e-16 (os termos ímpares cancelam contra
pesos simétricos). Os pontos de consulta são assimétricos de propósito. E uma
guarda do driver nasceu inoperante e caiu na primeira execução: o `make clean`
exigia `DIM`, falhava com aviso, e a limpeza por dimensão não acontecia.

---

## 11. Antes de comprometer o desenho: medir

Escalabilidade no regime 3D. O 2D rendeu 1,62× em np=2 e 1,68× em np=4, e isso é
tamanho do caso — 9377 células por quatro dão ~2300 por processo, e a comunicação
come o ganho. O domínio 3D da fração tem 375 mil células, quarenta vezes mais
trabalho por processo. **Sem essa medida, o escalonamento do t8Flow é promessa.**

Modos de medição vazia já catalogados neste projeto e que se aplicam: build sem
`-DHIGFLOW_COM_T8CODE` mede malha uniforme e parece medir t8code; processos
`ns-example` órfãos já falsearam uma medida por fator de cinquenta.
