# Escalonamento no cluster Euler — projeto e linha de base

Documento de **projeto**, não de resultado: registra as decisões tomadas, a
linha de base medida nesta máquina, e o que falta saber do ambiente. Os números
daqui são de uma estação de 8 núcleos — servem de ponto de partida e de
previsão a ser confrontada, não de resultado de cluster.

---

## 1. O alvo

**Artigo com o Poisson**, escalonamento **forte e fraco**, e o caso da cavidade
(ou equivalente) como aplicação de Navier–Stokes. Decidido pelo usuário em
19/10/2026, junto com **começar o NS 3D** — que foi feito e está verificado por
solução manufaturada (segunda ordem; ver `relatorio-ns-fv-2d.md` e o contrato,
cláusulas U1–U15).

O artigo já existe em `doc/artigo-t8flow/`, com seções `08-paralelo` e
`09-tresdimensoes`. **A convenção dele é que nenhum número é digitado no
texto**: as tabelas saem de `dados/*.dat`, gerados por `gera-dados.sh`. Estender
para o cluster é gerar `.dat` novos, não reescrever prosa.

---

## 2. O Euler: dois tipos de nó

| | núcleos | acelerador |
|---|---|---|
| nó A | **256** | — |
| nó B | **64** | **H100** |

Isso muda o planejamento em dois pontos, e os dois são consequência das medidas
da seção 3.

**256 núcleos num nó exige caso muito maior do que os de hoje.** O caso complexo
tem 162 mil células em **33 subdomínios**, e o teto natural de processos é 33 —
ainda ganhando em 8 nesta máquina, mas a 33% de eficiência. E o teto de
processos úteis é ainda mais baixo que 33 — é **6,75×**, posto pelo maior bloco
e não pela contagem deles (§4). Uma curva de escalonamento forte até 256 exige
malha gerada para isso, não a existente.

**O nó H100 não tem caminho de software hoje**, e é importante não contar com
ele por engano:

- o **PETSc desta árvore foi construído SEM CUDA** (verificado em
  `petscvariables`). Os solves lineares — que dominam o tempo — rodam em CPU;
- o t8Flow tem back-end de **OpenCL**, não CUDA, e ele cobre o **lote de
  estênceis do MLS**, não o solver linear. O ganho possível ali é sobre uma
  fração do tempo, não sobre o laço;
- o HiGFlow não tem caminho de GPU.

Usar o nó B hoje significa usar 64 núcleos de CPU e deixar a H100 ociosa. Fazer
a H100 render exige, no mínimo, PETSc com CUDA — o que é trabalho de ambiente,
não de código do solver. A discussão de prioridade está em
`HigFlow/doc/relatorio-projeto-cpu-gpu/`.

---

## 3. A linha de base

`example3d_complex`: 162 mil células, 33 subdomínios, monofásico newtoniano,
CUBISTA, projeção. Estação com **8 núcleos físicos** (i9-9880H).

**Dois comprimentos por `np`** (2 e 12 passos), e a diferença dá o laço sem o
arranque. Misturar os dois mediria a soma e atribuiria o resultado ao solver.
Sem gravar VTK nem salvamento durante a medida — escrita por rank é E/S que não
escala. Corridas **sequenciais**: duas disputando núcleo dão um tempo que não é
de nenhuma das duas.

**Seis repetições por ponto**, em duas varreduras independentes de três — a
segunda precedida de uma corrida de aquecimento descartada. Os valores são
**medianas**.

### A linha de base que vale, com `-pc_type hypre`

Três repetições, aquecimento descartado. `hypre` é o que o próprio caso
recomenda em `mesh-channel3D-script.sh`; a primeira linha de base usou
`bjacobi`, herdado do Makefile, e está logo abaixo com o que ela ainda ensina.

| `np` | laço (10 passos) | [mín, máx] | dispersão | ganho | eficiência | arranque |
|---|---|---|---|---|---|---|
| 1 | 207,1 s | 206 – 209 | 1,4% | 1,00× | 100% | 313,4 s |
| 2 | 139,9 s | 133 – 140 | 4,7% | 1,48× | 74% | 215,5 s (1,45×) |
| 4 | 97,6 s | 93 – 99 | 6,9% | 2,12× | 53% | 141,1 s (2,22×) |
| 8 | 67,1 s | 64 – 68 | 6,5% | **3,08×** | 39% | 93,5 s (3,35×) |

**A dispersão caiu de 21–35% para 1,4–6,9%.** Com `bjacobi` eu havia concluído
que os dois primeiros pontos eram ruidosos demais para publicar e que a causa
era desconhecida. Continua desconhecida em detalhe, mas está **ligada àquela
configuração**: as corridas de `hypre` são três vezes mais curtas e muito mais
reprodutíveis. Não é propriedade da estação, como eu supus.

### A primeira linha de base, com `bjacobi` — superada, e o que ela ensina

Seis repetições em duas varreduras independentes de três. **Os tempos absolutos
não valem como desempenho do caso** (§5); o que segue vale pelo método.

| `np` | laço (10 passos) | [mín, máx] | dispersão | ganho | eficiência |
|---|---|---|---|---|---|
| 1 | 388,0 s | 357 – 440 | 21% | 1,00× | 100% |
| 2 | 263,6 s | 225 – 318 | 35% | 1,47× | 74% |
| 4 | 187,1 s | 184 – 221 | 20% | 2,07× | 52% |
| 8 | 146,7 s | 144 – 150 | 4% | 2,64× | 33% |

As **eficiências batem** entre as duas configurações — 100/74/53/39 com `hypre`
contra 100/66/49/31 com `bjacobi`, a mesma forma de queda. É coerente com o teto
de banda de memória (§5) sendo a causa comum: ele não depende do
precondicionador.

### O que estes números dizem

**A curva não satura em 8.** De 4 para 8 o laço cai 40,4 s, contra espalhamento
de 6 s no `np = 8`. O ganho é pequeno em eficiência, mas é real e está fora do
ruído.

**O arranque escala**, e essa previsão minha estava errada: eu esperava que não,
porque o `relatorio-sistema` aponta que um processo materializa a malha global
antes de repartir. Medido, 259,2 → 81,1 s, **3,20×**. Com 33 arquivos, o laço de
leitura reparte arquivos entre processos. O impedimento daquele relatório vale
para caso de **árvore única** — a cavidade `example3d_lid_driven` cai nele.

### Duas conclusões que a repetição derrubou

Com **uma corrida por ponto** eu havia publicado nesta mesma seção: *97% de
eficiência em `np = 2`* e *"o joelho está entre 4 e 8, ganho nulo"*, esta
segunda com explicação causal pronta — 33 subdomínios sobre 8 processos, a
partição deixando de dividir trabalho novo.

**As duas eram ruído.** O `np = 2` dá 74%, e o ganho de 4 para 8 é sete vezes o
espalhamento. A explicação causal não era falsa por ser mal raciocinada: ela era
boa demais para um dado que não existia, e foi isso que a fez convincente.

### A dispersão, e a hipótese que não se sustentou

Os dois primeiros pontos dispersam 21% e 35%; os dois últimos, 20% e 4%.

Minha hipótese era **térmica** — a primeira medida da varredura 3 pegou a
máquina fria (295,5 s contra 339 e 361 depois), e como o protocolo **subtrai**
um par de corridas, basta que as duas metades estejam em estados diferentes para
o laço sair inflado por temperatura e não por cálculo.

Acrescentei o aquecimento descartado e **passei a registrar a temperatura do
pacote em toda corrida**, para não trocar uma afirmação por outra. O resultado:

- **O aquecimento não eliminou a dispersão** — no `np = 2` ela subiu, de 31%
  para 35%;
- **a temperatura não ordena os tempos**: no `np = 2` a corrida que partiu mais
  fria (67 °C) foi a mais lenta, e a mais quente (76 °C) a mais rápida;
- **o aquecimento consertou uma coisa**: o `np = 4` tinha uma corrida longa
  anômala de 304 s contra 352 na varredura 3; na 4 ficou 355/352/347.

Então a causa do resto **continua desconhecida**, e dizer isso é o estado
honesto. O que dá confiança no resultado não é a explicação, é a
**reprodutibilidade entre varreduras independentes**: as medianas de `np = 4` e
`np = 8` concordam dentro de 1% entre a varredura 3 e a 4, e as de `np = 1` e
`np = 2` dentro de 5% e 3%. A forma da curva é estável mesmo com os pontos
individuais ruidosos.

Um detalhe que o registro separado expõe e a média esconderia: no `np = 2` a
corrida **curta** é firme (190/200/199 s) e a **longa** é que balança
(426/518/467 s) — e o mesmo padrão aparece na varredura 3. Variabilidade que
cresce com a duração da corrida não é arranque; é algo que se acumula.

### O que estes números NÃO dizem

- **O `np = 8` é pessimista por motivo de máquina**: oito ranks em oito núcleos
  físicos, sem folga para o sistema. Parte dos 33% é contenção, e eu não separei
  as duas coisas.
- **Nada sobre escalonamento fraco** — todos os pontos são da mesma malha.
- **Nada sobre vários nós.** Tudo aqui é memória compartilhada; a rede do
  cluster não foi exercitada.
- **A dispersão de 20–35% nos dois primeiros pontos é alta demais para
  publicar.** Num nó dedicado de cluster deve cair; se não cair, o número de
  repetições tem de subir, e o protocolo já aceita isso por parâmetro.

---

## 4. A granularidade da partição, e as malhas geradas

### O teto não era 33 — é 6,75×

A seção 3 dizia que o teto de processos úteis deste caso era **33**, o número de
blocos. Está errado, e o número que vale é bem menor.

A higtree reparte **árvores inteiras**: `lb_add_input_tree()` recebe um
`hig_cell *`, e Zoltan distribui essas árvores. Então o **maior bloco é um piso
para o rank mais carregado**, e aqui ele tem 24 000 das 162 000 células —
**14,8% da malha**. Nenhum número de núcleos passa de **6,75×**, e o teto já é
atingido em `np = 8`:

| `np` | maior rank | fatia ideal | teto |
|---|---|---|---|
| 4 | 41 000 | 40 500 | 4,0× |
| 8 | 24 000 | 20 250 | **6,8×** |
| 64 | 24 000 | 2 531 | 6,8× |
| 256 | 24 000 | 632 | 6,8× |

Isso parecia recontextualizar a medida — em `np = 8` o melhor possível era 6,8×
e saiu 2,64×. **Parecia**: medido depois, levantar o teto não mudou o tempo, e a
granularidade não explica nada aqui. Ver o fim desta seção.

E corrige o que "gerar malhas maiores" significa: **refinar sozinho não levanta
o teto**. Multiplica todos os blocos igualmente, e o maior continua com 14,8%.

### O gerador: dois botões que não se confundem

`ci/tools/gera-malha-complexa.py`:

- **`--m`** — células por 0,1 de comprimento. Todas as coordenadas dos 33 blocos
  são múltiplos de 0,1, então qualquer `m` inteiro dá contagem inteira em todo
  bloco. `m = 10` reproduz a malha original, e isso é verificado: caixa e
  contagem idênticas nos 33 blocos e nos 169 retalhos, e o yaml de contorno sai
  byte a byte igual fora dos `path:`.
- **`--blocos` / `--para-np`** — corta o maior bloco ao meio, em índices
  **inteiros** de célula, então a união é preservada exatamente e não há
  arredondamento a conferir depois.

**As condições de contorno não são cortadas, e não precisam ser.** A busca em
`domain.c` é geométrica — varre `dirichlet_bcs` procurando células numa caixa em
torno do ponto — sem associação bloco↔contorno. Um retalho pode cobrir dois
blocos. Por isso os 169 retalhos mantêm número, ordem e id, e o `switch (id)`
compilado em `ns-complex-3d.c` continua válido sem reescrita.

### Cortar blocos é neutro para a física — medido, não suposto

Mesma malha (`m = 10`), `np = 1`, um passo, 33 contra 66 e contra 264 blocos:

| | maior diferença | % da escala do campo |
|---|---|---|
| `p` | 2,4e-06 | 1,6e-08% |
| `u`, pior componente | 1,4e-09 | 4,4e-08% |

Todos os pontos casaram um a um. A diferença está no nível de `ksp_atol = 1e-10`,
não no da discretização.

**Esse número só significa algo porque não veio do VTK.** O escritor interpola a
velocidade nos cantos de cada célula com `compute_facet_value_at_point`, que
chama `sfd_get_stencil` — a mesma máquina que, num ponto sem célula no domínio,
fecha por contorno distante. Neste domínio com vãos os cantos caem fora
justamente junto a eles, e a **escrita** fica contaminada onde a solução é
difícil: 0,96% pela lente do VTK contra 0,084% no campo instrumentado, num caso
já medido. Daí `HIGFLOW_DESPEJA` (valor bruto de `dpp`/`dpu`, coordenada como
chave) e `ci/tools/compara-despejo.py`.

### O teto de 40 árvores, que era real e deixou de ser

Com 66 blocos em `np = 1` o caso **abortou**:
`Number of HigTrees (66) is greater than MAXHIGTREESPERDOMAIN = 40`.

`sim_facet_domain.sfbi` era um vetor **fixo de 40** indexado pelo índice da
árvore **sem conferência de limite**; o guarda em `hig-flow-ic.c` existia para
isso e fez o certo — abortou alto em vez de deixar escrever fora da struct.

O comentário de `domain.h` trazia a medição que justificava o 40 (`np = 1` dava
33 árvores, "sete abaixo do teto") e concluía *"Treat 40 as a real ceiling"*. A
medição estava certa **para as malhas que existiam**; o que ela não cobria era
malha com mais blocos.

Agora `sfbi` cresce no único ponto que escreve (`sfd_set_sfbi`), como `higtrees`
já crescia em `sd_add_higtree`. `MAXHIGTREESPERDOMAIN` virou capacidade inicial.

**`MAXBCSPERDOMAIN = 200` continua sendo teto de verdade** — `dirichlet_bcs[]` é
alocado com ele e só um `assert` protege, e este caso usa 169. Como os retalhos
não são cortados nem multiplicados, nem o corte nem o refino encostam nisso. Mas
é estreito, e está anotado no cabeçalho.

### Um erro meu que o teste pegou, e o que ele ensina sobre `assert`

A primeira versão pôs `assert(i < sfd->sfbi_cap)` em `sfd_get_sfbi`.
`test-facet-domain-serial` morreu com SIGSEGV nos dois DIM.

Duas coisas de uma vez. **A higtree compila com `-DNDEBUG`** (`higtree/Makefile:18`),
então o `assert` não existe em produção e o que sobrou foi desreferência de
ponteiro nulo — eu havia olhado a linha de compilação do *exemplo*, que não tem
`-DNDEBUG`, e generalizado para a biblioteca. E **"fora de faixa devolve NULL" é
contrato**, não conveniência: o vetor fixo era zerado na criação, e o teste
verifica isso como "o estado de partida".

A proteção de verdade não é `assert` nenhum: é o vetor crescer no único ponto
que escreve, o que torna o índice fora de faixa **impossível por construção**.

Sinal que quase passei batido: a contagem foi de 162 para 170 casos e o contrato
de malha de 17/18 para 18/18 cláusulas. Os binários mortos levavam consigo os
casos que ainda não tinham rodado, e a cláusula P3 que o relatório acusou como
"sem teste que rode" era o meu defeito numa segunda contagem — eu li como lacuna
antiga.

Depois da correção: **170/170 e 18/18**, suíte **33/33**, e `example3d_complex`
com `--include-slow` **ok em `np = 1, 2, 3`**.

### E levantar o teto NÃO acelerou nada — a medida que derruba a minha premissa

Mesma malha (`m = 10`, 162 000 células), `np = 8`, três repetições com
aquecimento descartado e configurações intercaladas. Medianas do laço de 10
passos:

| blocos | laço | [mín, máx] | dispersão | teto | contra 33 |
|---|---|---|---|---|---|
| 33 | 153,1 s | 141,5 – 157,9 | 10,7% | 6,8× | — |
| 66 | 153,1 s | 150,2 – 159,6 | 6,1% | 46,3× | 0,0% |
| 264 | 157,4 s | 155,7 – 161,6 | 3,7% | 180,0× | +2,8% |

**O teto subiu 26 vezes e o tempo não se mexeu.** As três medianas estão dentro
da dispersão umas das outras.

E o dado que refutava a minha premissa estava na tabela que eu mesmo imprimi:
em `np = 8` o teto era **6,75×** e a medida era **2,64×**. O desequilíbrio de
partição vale no máximo 18% (24 000 contra 20 250 da fatia ideal); a lacuna até
2,64× é de **2,5 vezes**. Granularidade nunca foi o termo dominante aqui, e isso
dava para ler antes de cortar um bloco sequer.

**O que isso faz com o trabalho das malhas:** não o invalida, mas muda o que ele
é. 6,75× é teto **absoluto** — em 256 núcleos ele ata, e nenhuma curva até lá
existe sem malha mais fina. O gerador é **pré-condição para o regime de muitos
processos**, não ganho hoje. Vendê-lo como ganho seria exatamente o erro da
seção 3: explicação boa para um efeito que não foi medido.

**E fica uma pergunta aberta, que é a que importa:** se não é granularidade, o
que prende `np = 8` em 2,64×? Candidatos que o código sugere e nenhum medido —
a varredura de `dirichlet_bcs` é O(169) por ponto de estêncil, a troca de franja
cresce com o número de blocos, e oito ranks em oito núcleos físicos não deixam
folga para o sistema. **Isso precisa de perfilamento, não de outra malha**, e
vale resolver antes do Euler: o que prende em 8 provavelmente prende em 256.

Nota sobre um número que eu citei antes: em `np = 1` e corrida única, o passo
custou 41,2 s com 33 blocos e 54,2 s com 264, e eu tratei isso como "cortar
custa 32%". **Não foi confirmado** — em `np = 8` o custo de 264 blocos é 2,8%,
dentro do ruído. Os dois regimes são diferentes (em `np = 1` as 264 árvores caem
todas num rank), então não se contradizem; mas o 32% segue sendo corrida única e
não deve ser citado como medida.

Repete-se também o padrão de a **primeira medida sair rápida** (141,5 s contra
153,1 e 157,9), agora com aquecimento descartado antes. É a terceira varredura
em que aparece, continua sem explicação, e é mais uma razão para nunca publicar
ponto único.

### As famílias geradas

`conteiner/gera-familia-de-malhas.sh`. As malhas **não vão versionadas**: 4 318
arquivos de ~150 bytes que o roteiro reproduz exatos em segundos.

**Fraco** — `np = k³` e `m = 5k` dá **20 250 células por rank exatamente** em
todos os pontos, porque células `= 162·m³` e o cubo cancela. Com `np` em
potências de 2 o `m` inteiro não fecha e a carga por rank variaria até 12% entre
pontos; como a eficiência fraca é uma razão de tempos, isso vira 12% de erro
direto nela.

| `np` | 1 | 8 | 27 | 64 | 125 | 216 |
|---|---|---|---|---|---|---|
| `m` | 5 | 10 | 15 | 20 | 25 | 30 |
| células | 20 250 | 162 000 | 546 750 | 1 296 000 | 2 531 250 | 4 374 000 |
| blocos | 33 | 66 | 132 | 264 | 528 | 1056 |

O ponto `np = 8` é `m = 10`: a malha que já tem linha de base medida.

**Forte** — uma malha, `m = 32`, 5 308 416 células em 1056 blocos: 20 736 por
rank em 256, de propósito igual à carga do fraco, para os dois testes falarem da
mesma quantidade de trabalho por processo.

---

## 5. Perfilando o `np = 8`: o que prende, e não era nada do que eu procurava

A seção 4 terminou com a pergunta certa — se não é granularidade, o que prende
`np = 8` em 2,64×? Perfilado com `-log_view` do PETSc, **tirando o laço por
subtração de dois comprimentos dentro da mesma janela de medição**.

Essa última parte não é detalhe. Minha primeira leitura foi de uma corrida de
**3 passos** e disse que o solver linear era 26% do tempo. Era o arranque — que
é grande e escala bem — diluindo o solver. Com a subtração, `KSPSolve` é **56%
do laço em `np = 1` e 70% em `np = 8`**. O mesmo erro que o protocolo de dois
comprimentos existe para evitar, cometido por ler uma corrida única sem subtrair.

| `np` | laço | ganho | `KSPSolve` | ganho | fora do PETSc | ganho | iterações |
|---|---|---|---|---|---|---|---|
| 1 | 461,6 s | 1,00× | 259,8 s | 1,00× | 201,8 s | 1,00× | 13 238 |
| 2 | 384,8 s | 1,20× | 216,1 s | 1,20× | 168,7 s | 1,20× | +21,2% |
| 4 | 222,9 s | 2,07× | 164,8 s | 1,58× | 58,1 s | 3,47× | +21,4% |
| 8 | 186,8 s | 2,47× | 131,1 s | 1,98× | 55,7 s | 3,62× | +18,1% |

Corrida única por ponto, e o `np = 2` mostra o custo disso: dá 1,20× aqui contra
1,47× na linha de base de seis repetições. Nenhuma conclusão abaixo depende dele.

### Descartados por medida, não por argumento

- **Desequilíbrio de carga**: as razões máx/mín do PETSc entre ranks são 1,0–1,1
  em tempo, flops e memória. Confirma por outro caminho o que a varredura de
  granularidade disse.
- **Comunicação**: `VecScatterBegin` é 1,7% do laço em `np = 8` (3,2 s de 186,8).
- **O código fora do PETSc** — estêncil, contorno, CUBISTA, montagem: escala
  **3,62×**, melhor que o conjunto. A varredura O(169) de `dirichlet_bcs` que eu
  apontei como suspeita está dentro desses 3,62×, e não é ela que prende.

### O teto da própria estação, que nenhum perfil do solver mostra

Sem isto não há como dizer se 2,47× é o solver que não escala ou a máquina que
não entrega — duas conclusões opostas, e o número do solver sozinho não as
distingue. `conteiner/teto-da-maquina.c`, dois regimes:

| `np` | tríade (banda de memória) | núcleo (ULA, cabe no cache) |
|---|---|---|
| 2 | 1,21× | 1,63× |
| 4 | 1,26× | 3,63× |
| 8 | **1,31×** | **6,82×** |

Esta estação entrega **1,31×** de banda agregada com oito núcleos: satura em
~32 GB/s, e **um núcleo sozinho já tira 24,4 GB/s**. Trabalho limitado por ULA
escala 6,82×.

Matvec esparso e substituição triangular são limitados por banda. O ganho por
iteração do solver (2,34× em `np = 8`) fica **entre os dois tetos** — não é
desperdício de paralelismo no código. O `MatSolve` chega a 3,06× justamente
porque em `np = 8` a matriz local é oito vezes menor e passa a caber no cache.

**Consequência para o Euler:** este controle tem de rodar nos dois tipos de nó
**antes** de qualquer medida de solver. Um nó de 256 núcleos tem muito mais
banda absoluta, mas a razão banda/núcleo pode ser pior que a deste notebook, e é
ela que manda aqui.

### E o achado que não é sobre paralelismo

O precondicionador é o único fator do lado do código neste gargalo. Testado:

| pc | `np` | laço | `KSPSolve` | iterações |
|---|---|---|---|---|
| `bjacobi` | 1 | 461,6 s | 259,8 s | 13 238 |
| `bjacobi` | 8 | 186,8 s | 131,1 s | 15 634 |
| `asm` | 8 | 195,6 s | 155,8 s | 13 766 |
| **`hypre`** | 1 | **210,2 s** | **52,6 s** | **310** |
| **`hypre`** | 8 | **64,0 s** | **24,1 s** | **300** |

**43 vezes menos iterações.** O laço em `np = 8` cai de 186,8 s para 64,0 s, e o
ganho de paralelismo sobe de 2,47× para **3,28×**.

Verificado que é a **mesma solução**, e não convergência mais frouxa: campo por
coordenada contra o despejo do `bjacobi`, pior diferença **4,4e-08% da escala** —
o mesmo patamar da comparação de blocos, nível de `ksp_atol = 1e-10`.

E o `hypre` **não paga o penal de iterações**: 310 em `np = 1` contra 300 em
`np = 8`. Aquele +20% era específico do `bjacobi`, cuja qualidade cai com o
número de blocos; multigrid algébrico não depende disso.

Correção de algo que eu disse antes de ter a curva: afirmei que o penal de
iterações "piora com `np`" e que era o maior risco para 256 núcleos. **Não
piora** — é salto único ao sair de um bloco para muitos, e plano de 2 a 8
(+21%, +21%, +18%). E com `hypre` não existe.

### A linha de base foi medida com o precondicionador que o caso não recomenda

`bjacobi` veio do alvo `run` do Makefile do exemplo. Mas há **três receitas em
desacordo no mesmo repositório**:

| onde | `pc_type` |
|---|---|
| `example3d_complex/Makefile:104` | `bjacobi` |
| `ci/run_suite.py:79` | `bjacobi` |
| `amrs/Channel_3d_script/mesh-channel3D-script.sh:387` | **`hypre`** |

Peguei a que estava mais à mão. **Então a tabela da seção 3 — 388, 264, 187,
147 s — descreve uma configuração que o próprio caso não recomenda**, e os
tempos absolutos de lá não devem ser citados como desempenho do caso. A forma da
curva se sustenta (o `hypre` escala melhor, não pior); os números não.

### O que fica para resolver

1. ~~Refazer a linha de base com `hypre`~~ — **feita** (§3): `np = 8` em 67,1 s
   contra 146,7 s, ganho 3,08× contra 2,64×, e dispersão de 1,4–6,9% contra
   21–35%.
2. **Unificar as três receitas.** Duas apontam para `bjacobi` e uma para
   `hypre`, e a que vale é a que der o caso resolvido mais rápido com a mesma
   solução. Isso muda a suíte, então é mudança a decidir, não a aplicar de
   passagem.
3. **O risco em 256 muda de lugar, não desaparece.** `bjacobi` degrada com o
   número de blocos mas não comunica globalmente; AMG é o contrário — a
   construção da hierarquia e o nível mais grosso têm comunicação global, que é
   o ponto fraco clássico em muitos processos. O que prende aqui não é o que vai
   prender lá.
4. **Com `hypre`, o alvo de otimização se inverte**: em `np = 8` o laço é 64,0 s
   com `KSPSolve` em 24,1 s (38%), então **62% passa a ser o código fora do
   PETSc**, que escala 3,95×. Com o solver consertado, é ali que sobra tempo.

---

## 6. O que falta saber do ambiente

Cinco itens, e cada um muda o roteiro:

1. **Nós acessíveis de cada tipo**, e a rede entre eles;
2. **Escalonador e teto de tempo por fila.** Se o teto for baixo, o ponto de
   salvamento deixa de ser conforto e vira requisito — o do t8Flow está feito e
   verificado (reprodução bit a bit, ver `relatorio-bifasico-ft-2d.md` §7); o
   HiGFlow tem o seu próprio caminho de retomada;
3. **PETSc e t8code: há módulo, ou construo?** Qual MPI. Se eu construir o
   PETSc, é a oportunidade de ativar CUDA e tornar o nó B útil;
4. **Sistema de arquivos.** Um VTK por processo em centenas de processos é um
   padrão que alguns sistemas paralelos punem severamente;
5. **Até quantos núcleos a curva deve ir** — é isso que dimensiona a malha a
   gerar.

---

## 7. Sequência proposta

1. ~~Gerar malhas maiores~~ — **feito**, e o teto não era 33 e sim 6,75×, posto
   pelo maior bloco (§4). As duas famílias saem de
   `conteiner/gera-familia-de-malhas.sh`.
2. ~~Repetir a linha de base com três corridas por ponto~~ — **feito**, seis
   repetições em duas varreduras (§3). Derrubou duas conclusões.
3. **Poisson 3D do t8Flow** como veículo do artigo: já é paralelo e verificado
   em `np = 1, 2, 4`, e o escalonamento dele não depende do NS.
4. **Antes de qualquer tempo, verificar que o resultado paralelo é idêntico ao
   serial em caso NÃO exato.** No caso exato tudo é zero de máquina e nada
   discrimina — e foi assim que três defeitos de franja sobreviveram verdes
   (ver `t8flow-franja-paralela`).
5. ~~Escalonamento fraco exige a família de malhas~~ — **feita** (§4), com
   células por processo exatamente constante.
6. ~~Perfilar o `np = 8`~~ — **feito** (§5). Não era `dirichlet_bcs`, nem
   desequilíbrio, nem comunicação: é banda de memória sobre álgebra esparsa,
   mais um precondicionador que o caso não recomenda. `hypre` corta o laço de
   `np = 8` a um terço.
7. **Rodar `conteiner/teto-da-maquina.c` nos dois tipos de nó do Euler antes de
   qualquer medida de solver.** Sem a razão banda/núcleo do nó, nenhum número de
   escalonamento lá se interpreta — foi esse controle que mostrou que 2,47× nesta
   estação é majoritariamente a máquina.
8. **Decidir qual das três receitas de `pc_type` vale** (§5). Muda a suíte.

---

## 8. Ressalva de método

Comparar VTK entre `np` diferentes **não é teste válido** em domínio
multibloco: o desempate do WLS depende da ordem de visita das árvores, que muda
com a partição. Para comparar decomposições, instrumentar os campos por faceta
com a coordenada como chave — isso deu 1,0e−10 contra 0,56% do VTK.

---

## 9. O que está preparado, e o que falta confirmar

**O contêiner já existe** (`conteiner/`): Apptainer com Ubuntu 26.04, OpenMPI 5,
PETSc 3.25.4, HYPRE, t8code 4.0.0 e HiGTree/HiGFlow, com `preparar-fonte.sh`
(empacota por `git archive`, só o versionado) e `verificar-euler.sh` (portões em
ordem crescente de custo). A construção do contêiner está marcada para a semana
de 26/10/2026.

**Acrescentado em 20/10/2026:**

| arquivo | o que é |
|---|---|
| `conteiner/escalonamento-euler.sh` | o protocolo: forte e fraco, repetições com mediana e dispersão, dois comprimentos por ponto, sequencial |
| `conteiner/euler-escalonamento.slurm` | submissão, com os dois tipos de nó marcados |

O protocolo começa por um **portão 0**: o paralelo tem de dar o mesmo que o
serial **num caso não exato**, antes de cronometrar. No caso exato tudo é zero
de máquina e nada discrimina — foi assim que três defeitos de franja
sobreviveram a uma suíte verde. E se o portão não conseguir ler os erros, ele
**falha** em vez de passar: um portão que não consegue olhar não é um portão que
aprovou.

### As incógnitas, agora numeradas junto com as do contêiner

`[E1]`–`[E4]` já estavam na definição do contêiner (ISA dos nós, PMIx/Slurm,
interconexão, montagens). Acrescentei:

- **`[E5]` nomes e limites das partições** — `sinfo -s` no nó de submissão.
  Palpite errado numa diretiva do Slurm **não falha**: a fila aceita, o job roda
  num nó que não era o pretendido, e o número sai sem dono.

E o `[E2]` ganhou uma guarda explícita no script de submissão, porque o modo de
falha dele é traiçoeiro: **se o PMIx não casar, cada rank acha que é o rank 0** —
a corrida "funciona", e mede N cópias do problema serial. O resultado tem a
forma de uma curva de escalonamento, com ganho 1,0. A guarda imprime o tamanho
do comunicador antes de qualquer medida.

### O que ainda não está preparado

**As malhas grandes.** O teto de 33 subdomínios do caso complexo satura em 8
processos — nenhuma curva até 256 faz sentido sem malha gerada para isso. Para o
Poisson do t8Flow o problema não existe (a floresta se reparte na criação), e é
por isso que ele é o veículo do artigo.

**O escalonamento fraco precisa de critério decidido.** O script usa "trabalho
por processo constante", o que em 3D significa `np` octuplicando a cada nível —
logo a lista útil é 1, 8, 64. Com 256 núcleos caberia um quarto ponto se a
malha base for pequena o bastante para o primeiro ponto caber num processo.
