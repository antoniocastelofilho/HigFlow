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
ainda ganhando em 8 nesta máquina, mas a 33% de eficiência. Uma curva de
escalonamento forte até 256 exige malha gerada para isso, não a existente.

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

## 3. A linha de base, medida em 20/10/2026

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

| `np` | laço (10 passos) | [mín, máx] | dispersão | ganho | eficiência |
|---|---|---|---|---|---|
| 1 | 388,0 s | 357 – 440 | 21% | 1,00× | 100% |
| 2 | 263,6 s | 225 – 318 | 35% | 1,47× | 74% |
| 4 | 187,1 s | 184 – 221 | 20% | 2,07× | 52% |
| 8 | 146,7 s | 144 – 150 | 4% | 2,64× | 33% |

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

## 4. O que falta saber do ambiente

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

## 5. Sequência proposta

1. **Gerar malhas maiores** para o caso de escalonamento. O teto de 33
   subdomínios é do caso, não do método, e nenhuma medida até 256 faz sentido
   sem isso.
2. ~~Repetir a linha de base com três corridas por ponto~~ — **feito**, seis
   repetições em duas varreduras (§3). Derrubou duas conclusões.
3. **Poisson 3D do t8Flow** como veículo do artigo: já é paralelo e verificado
   em `np = 1, 2, 4`, e o escalonamento dele não depende do NS.
4. **Antes de qualquer tempo, verificar que o resultado paralelo é idêntico ao
   serial em caso NÃO exato.** No caso exato tudo é zero de máquina e nada
   discrimina — e foi assim que três defeitos de franja sobreviveram verdes
   (ver `t8flow-franja-paralela`).
5. **Escalonamento fraco** exige a família de malhas crescendo com `np`; é
   trabalho de geração de malha, e vale decidir o critério (células por
   processo constante) antes de gerar.

---

## 6. Ressalva de método

Comparar VTK entre `np` diferentes **não é teste válido** em domínio
multibloco: o desempate do WLS depende da ordem de visita das árvores, que muda
com a partição. Para comparar decomposições, instrumentar os campos por faceta
com a coordenada como chave — isso deu 1,0e−10 contra 0,56% do VTK.

---

## 7. O que está preparado, e o que falta confirmar

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
