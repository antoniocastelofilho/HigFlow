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
já saturado em 8 nesta máquina. Uma curva de escalonamento forte até 256 exige
malha gerada para isso, não a malha existente.

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

| `np` | total (12 passos) | arranque | laço (10 passos) | ganho no laço | eficiência |
|---|---|---|---|---|---|
| 1 | 719,0 s | 189,9 s | 440,9 s | 1,00× | 100% |
| 2 | 419,5 s | 146,6 s | 227,5 s | **1,94×** | **97%** |
| 4 | 265,9 s | 76,1 s | 158,1 s | 2,79× | 70% |
| 8 | 222,3 s | 34,7 s | 156,3 s | 2,82× | 35% |

### O que estes números dizem

**97% em dois processos** — contra 41% registrados para o caso de AMR dinâmico
no `relatorio-sistema`, e contra 1,62× do caso 2D viscoelástico. A diferença é
**tamanho**: 162 mil células contra 9 mil. Confirma a previsão de
`higflow-paralelo-estado`: o regime em que a divisão vence a comunicação começa
quando há trabalho que baste por processo.

**O joelho está entre 4 e 8.** De 4 para 8 o laço vai de 158,1 s a 156,3 s —
ganho nulo. Com 33 subdomínios e 8 processos, cada um fica com ~4 subdomínios, e
a partição deixa de dividir trabalho novo. **O limite é a decomposição do caso,
não o número de núcleos.**

**O arranque escala, e essa previsão minha estava errada.** Eu esperava que não
— o `relatorio-sistema` aponta que um processo materializa a malha global antes
de repartir. Medido: 189,9 → 34,7 s, ganho de 5,5×, *melhor* que o laço. Com 33
arquivos, o laço de leitura reparte arquivos entre processos e isso domina. O
impedimento descrito naquele relatório vale para caso de **árvore única** — a
cavidade `example3d_lid_driven` tem um domínio só e cai nele.

### O que estes números NÃO dizem

- **Sem dispersão.** Uma corrida por ponto, dez passos. Para publicar, no mínimo
  três repetições por ponto.
- **O `np=8` é pessimista por motivo de máquina**: oito ranks em oito núcleos
  físicos, sem folga para o sistema. Parte dos 35% é contenção, e eu não separei
  as duas coisas.
- **Nada sobre escalonamento fraco** — todos os pontos são da mesma malha.
- **Nada sobre vários nós.** Tudo aqui é memória compartilhada; a rede do
  cluster não foi exercitada.

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
2. **Repetir a linha de base com três corridas por ponto**, aqui, antes de
   migrar. Barato, e dá a dispersão que falta.
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
