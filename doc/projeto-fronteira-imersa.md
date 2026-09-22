# Fronteira imersa com força direta

Projeto, não implementação. O terreno abaixo foi **medido** no código, não
suposto — cada afirmação sobre o que existe hoje tem o arquivo e a linha.

As três decisões tomadas em 21 de setembro de 2026:

```
corpo          RÍGIDO E FIXO          →  U_corpo = 0
força          DEFASADA               →  calculada com a velocidade do passo anterior
alcance        SÓ O NEWTONIANO        →  hig-flow-step.c, não os doze step-*.c
```

Elas colapsam o problema, e vale dizer exatamente o quanto.

---

## 1. O que as decisões eliminam

**Corpo fixo elimina a migração.** Os marcadores nunca mudam de posição, logo
nunca mudam de rank. Toda a questão de estrutura distribuída — que motivou as
sondas do `DMSwarm` e do `DMPlex` — deixa de ser necessária nesta primeira
versão. Os marcadores são gerados uma vez e ficam.

**Corpo rígido elimina a deformação**, e com ela a necessidade de topologia: sem
área dual variável nem normal recalculada, o peso de cada marcador é constante.
Em 3D isso é o que dispensa o `DMPlex`; em 2D nunca foi preciso.

**Força defasada elimina o estágio novo.** É a diferença entre tocar um ponto e
tocar doze funções de passo — ver a seção 3.

Nenhuma das três é irreversível, e a seção 7 diz o que cada uma custa para
desfazer.

---

## 2. O terreno, medido

**O lado euleriano já existe e está vivo.** Há um campo de força por faceta,
`ns->dpFU[dim]`, com subdomínio próprio (`psfdF[dim]`), preenchido a cada passo
por `higflow_calculate_facet_source_term` (`hig-flow-step.c:415`) e lido de volta
na equação:

```
dpFU[dim]  →  cc.F              hig-flow-discret.c:133
cc.F       →  higflow_source_term    hig-flow-terms.c:29
                                →  lado direito da quantidade de movimento
```

A cadeia foi percorrida inteira, e isso não é zelo excessivo: neste código há
campos preenchidos que ninguém lê. Este é lido.

**A ordem do passo** (`higflow_solver_step`, `hig-flow-step.c:1326`):

```
contorno → fonte (analítica, f(x,t)) → preditor u* → pressão → projeção
```

A fonte nasce **antes** do preditor e só depende de posição e tempo.

**A localização de ponto existe e é barata**: `sfd_get_facet_with_point`
(`domain.h:481`), e o instantâneo a tornou barata para geometria. As cláusulas
C7–C9 do contrato de Mesh a afirmam nos dois backends.

**O `dp_sync` copia dono→franja e sobrescreve** (`pdomain.c:1308`). Não existe
acumulação franja→dono. Isso decide a seção 4.

---

## 3. A formulação, e onde ela entra

Corpo rígido e fixo: a velocidade desejada na fronteira é zero. A força direta
defasada, no marcador `k`:

```
f_k = ( 0 − u(x_k, t^n) ) / Δt
```

com `u(x_k, t^n)` interpolada do campo do passo **anterior** — que é o que
"defasada" significa e o que dispensa esperar o preditor.

**Entra em exatamente um lugar**: logo após a chamada de
`higflow_calculate_facet_source_term` em `higflow_solver_step`, somando em
`dpFU[dim]` antes que o preditor o leia. Uma chamada nova, num arquivo.

O preço, dito sem rodeio: **o não escorregamento não é imposto no passo
corrente.** A velocidade na fronteira não vai a zero exatamente; vai a `O(Δt)`
do zero. É a diferença entre isto e o estágio pós-preditor do Uhlmann, e é a
razão de a seção 6 medir justamente esse resíduo em vez de olhar figura.

---

## 4. A estrutura de dados: as malhas do PETSc

**Decidido: `DMSwarm` para os marcadores, `DMPlex` para a topologia quando ela
for necessária.** Não é dependência nova — o PETSc 3.25.4 já é dependência, e as
sondas de `higflow/fronteira-imersa/sondas/` mediram as duas capacidades nesta
instalação (migração com carga intacta; superfície triangulada distribuída com o
suporte do vértice preservado).

Registro que **eu havia recomendado replicar**, e por que a decisão contrária é
melhor. Meu argumento era de tamanho: nestas malhas `h = 0,01`, um cilindro de
diâmetro 0,05 dá ~160 marcadores, e replicar 240 KB por rank não custa nada.
Está certo — **para a primeira versão**. O que ele ignora é a trajetória: corpo
móvel exige migração e corpo deformável exige topologia distribuída, e chegar
neles com marcadores replicados significa trocar a estrutura depois, com o
método já funcionando em cima dela. Otimizar a primeira versão ao preço da
segunda é o tipo de economia que se paga com juros.

**Mas a escolha ressuscita um problema que a replicação fazia sumir**, e ele
precisa de resposta explícita, porque descobri-lo tarde é caro:

Com marcadores distribuídos por posse euleriana, o espalhamento de um marcador
perto da borda escreve em facetas de **outro** rank. O `dp_sync`
(`pdomain.c:1308`) manda dono→franja e **sobrescreve**: contribuição escrita numa
franja é descartada, calada. Acumulação franja→dono não existe na HiGTree.

**A resposta está na mesma biblioteca.** `PetscSF` é exatamente a primitiva de
comunicação irregular para isto, e `PetscSFReduce` com `MPI_SUM` faz a
acumulação na direção que falta. O grafo se monta a partir do que a HiGTree já
tem — `gid_map`, `local_count`, `total_count` e os vizinhos filtrados de
`_dp_shared` (`pdomain.h:138`) — sem tocar na estrutura dela: o `PetscSF` fica
por fora, lendo o mapa que já existe.

**Sondado, e fecha** (`sondas/sonda-petscsf.c`). Cada rank escreve 1,0 nas
entradas que possui e 10,0 numa entrada de franja do vizinho; depois do
`PetscSFReduce` com `MPI_SUM`, o dono lê **11,0** — a sua contribuição mais a do
vizinho — em np=2 e np=3. Se o `PetscSF` sobrescrevesse em vez de somar, daria
10,0, e o teste distingue os dois casos.

É a primitiva que faltava, e ela vem da mesma biblioteca que já decidimos usar.

## 5. O suporte do marcador, sem varredura

O delta regularizado tem suporte de três células por direção (Roma). Achar as
facetas do suporte **varrendo** seria `O(N_L × N_facetas)`.

Não é preciso: a malha é uniforme e o espaçamento é conhecido, então os centros
das facetas do suporte estão em **deslocamentos conhecidos** a partir do
marcador. Cada um se acha com uma `sfd_get_facet_with_point`. Custo
`O(N_L × 3^DIM)` localizações, cada uma logarítmica.

Isso também fixa a exigência de franja **euleriana**: o suporte se estende uma
célula e meia para cada lado, então um marcador a menos disso da borda do rank
tem suporte atravessando a fronteira. Essa franja é a da HiGTree, e **não** é a
sobreposição do `DMPlexDistribute` — as duas foram conflacionadas numa versão
anterior deste texto. A sobreposição do Plex é da malha lagrangeana e serve a
outra coisa (medida: `sondas/sonda-plex-sobreposicao.c`; sobreposição 1 basta
para área e normal do vértice, e custa +19% de células em np=2 e até +53% em
np=3 numa esfera de 1280 triângulos). A largura de franja euleriana foi **medida**
(`sondas/sonda-franja-euleriana.c`): o alcance para fora de uma célula própria é
exatamente a franja declarada — 1, 2 e 5 células para franja 1, 2 e 5. O HiGFlow
declara 5 (`hig-flow-kernel.c:1465`), então Roma (1,5) e Peskin (2,0) cabem com
folga, e **a franja não é restrição**. Mas a HiGTree usa 1 por omissão: com ela o
suporte não caberia, e a perda apareceria como força um pouco menor junto às
fronteiras de partição — defeito que se lê como "efeito de malha". Com replicação isso é inofensivo na interpolação (a
soma global recupera as parcelas), mas **no espalhamento cada rank precisa ver
as facetas do suporte que ele possui** — o que já é verdade, porque ele as
possui.

---

## 6. Como verificar — e por que não olhando figura

O erro clássico aqui é validar com a esteira de um cilindro: ela *parece* certa
numa faixa larga de implementações erradas. Três oráculos que discriminam:

**1. O resíduo de não escorregamento.** Medir `max_k |u(x_k)|` ao longo do
tempo. Com força defasada ele não vai a zero, vai a `O(Δt)` — então o teste é
**refinar `Δt` e ver o resíduo cair na mesma ordem**. Um erro de sinal ou de
peso não passa nesse teste; a figura da esteira passa.

**2. Parede imersa alinhada com a malha, num canal.** Pôr a fronteira imersa
sobre um plano coordenado dentro de um canal e comparar com Poiseuille no canal
reduzido, que tem **solução analítica**. Alinhada com a malha, o delta degenera
para o caso mais simples e o resultado tem de bater com a solução, não só
parecer plausível.

**3. Conservação da força total.** `Σ_k f_k · w_k` espalhada tem de igualar a
integral da força euleriana. É o que pega dupla contagem — e a dupla contagem
tem um caminho concreto: com replicação, todo rank varre todos os marcadores, e
um marcador exatamente sobre a face entre células de ranks diferentes é contado
duas vezes. É o mesmo desempate que a cláusula C8 exige e que a decisão D3
resolveu para classificação de ponto; aqui ele deixa de ser convenção e vira
conservação.

O caso de suíte seria 2D, uniforme, rápido — e portanto também rodável pelas
fontes do t8code, o que o põe sob a mesma referência dos demais.

---

## 7. O que fica de fora, e o que custa trazer de volta

| deixado de fora | o que custa depois |
|---|---|
| corpo móvel | migração de marcador: as sondas já mostram o caminho (`DMSwarm`) |
| corpo deformável | topologia distribuída (`DMPlex`) e peso recalculado por passo |
| os outros onze `step-*.c` | a chamada nova em cada um; o cálculo em si não muda |
| estágio pós-preditor (Uhlmann) | um estágio novo por função de passo, e iteração |
| geometria vinda de arquivo | ver abaixo |

**A quarta decisão, ainda sua.** Como o corpo entra. Assumo, para a primeira
versão, **gerado no próprio exemplo** a partir de parâmetros analíticos (centro,
raio, número de marcadores) — porque formato de arquivo é compromisso que merece
esperar o método funcionar. Se preferir arquivo desde já, diga: muda a entrada,
não o resto.

---

## 8. Dois casos que o nome esconde

"Fronteira imersa" cobre duas coisas que compartilham a maquinaria e divergem em
quase tudo o mais. A distinção precisa estar no código, não só na cabeça de quem
o escreveu.

```
(a) CONTORNO RÍGIDO CURVO        um sólido dentro do fluido: cilindro, perfil,
                                 parede que a malha cartesiana não representa
(b) INTERFACE ENTRE DOIS FLUIDOS bolha, gota, superfície livre
```

Na literatura são linhagens diferentes: (a) é fronteira imersa com forçamento
direto (Fadlun, Uhlmann); (b) é **rastreamento de frente** (Unverdi–Tryggvason).
Nomear os dois como a mesma coisa é o começo do erro.

### O que é comum — e é o que já está construído

Os dois operadores de transferência, o núcleo regularizado, o par adjunto, e a
acumulação franja→dono. `fi_interpola` e `fi_espalha` servem aos dois casos sem
uma linha de diferença. Isso não é coincidência: a transferência é sobre
*geometria e quadratura*, não sobre física.

### O que difere — quatro coisas, e elas se acumulam

**1. De onde vem a força.**

Em (a) a força é um **multiplicador de Lagrange**: vale o que for preciso para
que `u = U_corpo`. Não tem lei constitutiva, e sua magnitude cresce como `1/Δt`
— é rígida por construção.

Em (b) a força é **constitutiva**: tensão superficial `σ κ n`. Tem magnitude
física própria e traz a sua própria restrição de passo, a da onda capilar,
`Δt ≲ sqrt(ρ h³ / 2πσ)` — que não existe no caso (a).

**2. Os marcadores se movem?**

Em (a), com corpo fixo, nunca. É por isso que a estrutura atual basta: os
marcadores são colocados uma vez e ficam.

Em (b) eles são **advectados** pela velocidade interpolada, `dX/dt = u(X)`. Daí
saem duas exigências que hoje não existem: migração entre ranks quando o
marcador cruza fronteira de partição — medida na `sonda-dmswarm`, ainda não
usada — e **remalhamento**, porque a malha lagrangeana se deforma e o
espaçamento deriva. Quando `ds > h` o suporte do núcleo deixa de cobrir a
interface e ela fica **permeável**: o fluido atravessa. É o defeito mais
característico de rastreamento de frente, e não se anuncia como erro.

**3. Precisa de conectividade?**

Em (a), não. Os pesos são fixos na criação e cada marcador é independente.

Em (b), **sim, e é a diferença que quebra a estrutura atual.** Curvatura exige
vizinhos: em 2D a ordem ao longo da curva, em 3D a triangulação — o `DMPlex`,
com sobreposição 1, que é o que a `sonda-plex-sobreposicao` mediu.

E aqui está a consequência afiada: **distribuir marcadores por posse euleriana
destrói a ordem da curva.** O `fi_cria_curva` de hoje guarda só os marcadores
deste rank, sem ordem e sem ligação com os vizinhos. Para (a) isso é correto e
barato. Para (b) é exatamente o que não pode acontecer.

**4. Acoplamento com o resto da física.**

Em (a), nenhum: a densidade é uniforme e a força entra no campo por faceta.

Em (b), densidade e viscosidade **saltam** através da interface, e a formulação
de um fluido só precisa de uma função indicadora derivada das posições dos
marcadores. Mais que isso: **o HiGFlow já tem multifásico**, por VOF, com
curvatura própria (`vof-*-normal-curvature`, função altura) e um termo de tensão
interfacial já na equação — `higflow_interfacial_tension_term`, que lê `cc.IF` e
divide por `We * ρ`. Rastreamento de frente seria uma representação
**alternativa** de interface, não um acréscimo ao trabalho de fronteira imersa.

Uma armadilha concreta desse acoplamento: `higflow_source_term` **divide `cc.F`
por `cc.dens`** quando `flowtype == MULTIPHASE`. Uma força espalhada em `dpFU`
numa corrida multifásica já sai dividida por ρ — o que é o certo para a forma
não conservativa da quantidade de movimento, e é um erro de fator ρ para quem
não souber.

### O que isto implica para o código

Separar o que transfere do que **decide a força**:

```
fi_interpola / fi_espalha        comuns aos dois casos, prontos
fi_forca_corpo_rigido(c, dt)     caso (a) -- o multiplicador, f = (U - u)/Δt
fi_forca_tensao(c, sigma)        caso (b) -- exige curvatura, exige topologia
fi_move(c, dt)                   caso (b) -- advecção, migração, remalhamento
```

O caso (a) está construído e verificado. O caso (b) **não é uma extensão dele**:
precisa de topologia distribuída, de marcadores que migram, de remalhamento, e
de uma decisão sobre conviver com o VOF que já existe. É projeto próprio.

### A decisão, tomada em 21/09/2026: alternativa ao VOF, e separada

O rastreamento de frente é **alternativa** ao VOF, não complemento, e fica
**separado** dele. Não entra nos módulos `hig-flow-vof-*`, não usa a fração
volumétrica, e não passa pelo `cc.IF`.

Isso tem uma consequência que não é óbvia e que vale antecipar, porque
descobri-la no meio da implementação custaria caro:

**Separado do VOF, ele não herda a densidade.** A formulação de um fluido
precisa de ρ e μ variáveis através da interface, e hoje quem os produz é a
maquinaria do VOF, a partir da fração volumétrica. Um rastreamento de frente
separado tem de derivar os seus **da posição dos marcadores** — uma função
indicadora própria.

E isso decide como a força tem de ser escalada, por causa de um detalhe já
medido: `higflow_source_term` divide `cc.F` por `cc.dens` **apenas quando**
`flowtype == MULTIPHASE`. Então das duas, uma:

- o rastreamento de frente declara um `flowtype` próprio, e aí **ele mesmo**
  divide a força por ρ antes de espalhar; ou
- ele reusa `MULTIPHASE` — e aí herda o divisor, mas também a expectativa de que
  a densidade venha do VOF, que é justamente o que se quis evitar.

A primeira é coerente com "separado". Fica escrita como o que decidir primeiro
quando (b) começar, não como coisa a descobrir depois.

**O que a decisão NÃO muda:** os dois operadores de transferência continuam
comuns aos dois casos. Separar (b) do VOF não o separa da maquinaria que o caso
(a) já construiu e verificou.

---

## 9. O maior `Δt` estável — medido, e a pergunta estava errada

`example2d_SchaeferTurek`, 400 passos para **todos** os valores, np=3. Número
fixo de passos e não tempo fixo, porque instabilidade de realimentação cresce
**por passo**: a tempo fixo cada `Δt` teria um número diferente de
oportunidades de amplificação, e o passo grande pareceria melhor do que é.

Reproduzir com `ci/tools/varre-dt-fronteira-imersa.sh`.

```
Δt        max|u|      veredito    resíduo de não escorregamento
0,001     1,790       estável     4,20e-02
0,002     1,815       estável     8,93e-02
0,005     1,961       estável     2,24e-01
0,01      1,896       estável     3,98e-01
0,02      1,794       estável     6,45e-01
0,05      1,592       estável     1,01e+00
0,10      1,1e+104    DIVERGIU    1,49e+103
0,20      6,2e+95     DIVERGIU    1,06e+93
0,40      2,0e+109    DIVERGIU    2,25e+105
```

### A conclusão é o contrário do que motivou o experimento

**O despenhadeiro está entre 5e-2 e 1e-1, ACIMA do CFL convectivo (~0,026).**
Então a fronteira imersa **não é o gargalo de estabilidade** — o esquema
semi-implícito absorve `Δt` até quase o dobro do CFL. O teto explícito que eu
queria dar ao `higflow_adjust_timestep` não é necessário: o CFL já é mais
restritivo.

**Mas o limite útil está muito antes, e é de acurácia.** A coluna do resíduo
cresce como `O(Δt)` — razões 2,13 / 2,51 quando `Δt` dobra — e depois **satura**,
porque não pode passar da escala de velocidade local:

```
Δt = 1e-3    resíduo  2,3% de max|u|
Δt = 5e-3            12%
Δt = 2e-2            36%
Δt = 5e-2            63%   <- ultimo estavel, e o corpo deixou de ser corpo
```

Em `Δt = 5e-2` a corrida é estável e **o obstáculo é uma sugestão**. É o pior
desfecho possível para quem confia em estabilidade como critério.

**O que nenhum controlador de CFL captura.** Ele deixaria `Δt` crescer até 0,026
e entregaria um cilindro com ~20% de escorregamento sem reclamar de nada. Se o
`higflow_adjust_timestep` for ligado com a fronteira imersa ativa, o limite tem
de vir do **resíduo**, não do CFL — e isso é escolha de quem usa, porque depende
de quanta violação o resultado tolera.

### O que a tabela NÃO estabelece

- A primeira versão do varrimento (até 2e-2) **não achou limite nenhum** e uma
  tabela toda verde parece resposta. O limite só apareceu estendendo até 0,4.
  Varrimento que não encontra a fronteira não mediu o que se pensa.
- Nada sobre o `Cd`. Resíduo pequeno é condição necessária, não suficiente: a
  próxima medida é o arrasto em dois `Δt` e ver se ele **se move**.
- Nada em 3D. O despenhadeiro lá pode estar em outro lugar, e a malha é mais
  grosseira (10 células no diâmetro contra 20).

---

## 10. O método de segunda ordem da literatura — e o que ele não resolve

Levantamento de 21/09/2026, pedido a partir dos trabalhos de Aristeu da Silveira
Neto (UFU) e Alexandre Roma (IME-USP). **As quatro hipóteses que eu havia
levantado estavam todas erradas**, e uma delas foi verificada como
explicitamente não sendo o caso — registro isso porque o valor do levantamento
foi justamente ter pedido fonte em vez de confirmação.

### O método

**Griffith, Hornung, McQueen & Peskin, "An adaptive, formally second order
accurate version of the immersed boundary method", JCP 223 (2007) 10–49.**

Citado nominalmente no único artigo conjunto de Roma **e** Silveira-Neto em
revista — Ceniceros, Roma, Silveira-Neto & Villar, *Commun. Comput. Phys.*
8(1):51–94, 2010 — em contraste com o Roma-Peskin-Berger de 1999, que é de
**primeira** ordem e é de onde vem o núcleo de 3 pontos usado aqui.

Base não adaptativa: **Lai & Peskin, JCP 160 (2000) 705–719**, o artigo canônico
de "IB de segunda ordem".

O **Modelo Físico Virtual** (Lima e Silva, Silveira-Neto & Damasceno, JCP 189,
2003) **não é** o método de segunda ordem: Campregher/Silveira-Neto (2009)
registra que ele *não* impõe o não escorregamento diretamente, e ali
"second-order" se refere ao esquema espacial do solver.

### Segunda ordem em quê — a distinção que muda o plano

Segunda ordem **no esquema**, em tempo e espaço. **Não** na imposição da condição
de contorno. Lai & Peskin, literal: *"the use of the discrete delta function in
the boundary-fluid interaction prevents the immersed boundary method from being
more than first-order accurate"* — a velocidade não tem derivada contínua
através da fronteira. O ganho real é **menos viscosidade numérica**.

O próprio artigo do Roma/Silveira-Neto mediu isso separadamente: **ordem 1 na
vizinhança da interface, tendendo a 2 longe dela.**

Consequência direta para este projeto: a defasagem temporal da força e a primeira
ordem espacial do núcleo são **dois erros independentes**. Passar para um esquema
de segunda ordem mata o erro dominante de hoje, o `O(Δt)`, e **não** faz o resíduo
convergir melhor que `O(h)` com o refino de malha.

### As três receitas, e o custo de cada uma

```
Uhlmann / Fadlun    forca da velocidade PROVISORIA do passo corrente, nao de u^n.
                    Formula identica a' nossa; muda de onde vem a velocidade.
                    Uhlmann usa O MESMO nucleo de 3 pontos do Roma.
                    -> menor custo para eliminar o nosso O(dt)

Lai-Peskin          ponto medio, DUAS avaliacoes da forca por passo, ZERO iteracoes

Multi-direct        ITERA.  E o motivo nao e' a defasagem temporal: e' que
 forcing            espalhar e interpolar nao comutam (H.A != I).
                    20 iteracoes no classico, 5 em variantes aceleradas
```

### Os valores de referência do 2D-1

Verificados em duas fontes independentes — Nabh, tese, Universität Heidelberg,
Preprint 42/98, 1998, p. 74; confirmados pelo Featflow/TU Dortmund e por John &
Matthies, *IJNMF* 37 (2001) 885–903:

```
Cd = 5,57953523384     Cl = 0,010618948146     dp = 0,11752016697
dp medido entre (0,15 ; 0,2) e (0,25 ; 0,2)
```

### DUAS ADVERTÊNCIAS QUE MUDAM O QUE MEDIR

**1. `Cl` não é resolvível aqui.** O deslocamento do cilindro fora do eixo é
`D/20 = 0,005`. O suporte do núcleo de 3 pontos é `3Δx ≈ 0,15 D` — **maior que o
próprio deslocamento**. Com forçamento difuso em resolução moderada, `Cl` serve
de diagnóstico de simetria, **não** de critério de acurácia. Usar `Cd` e `Δp`.

**2. Vinte células por diâmetro é pouco.** **Não há na literatura nenhum `Cd` do
2D-1 obtido com fronteira imersa de forçamento contínuo** — lacuna declarada, não
suposta. O que há é indireto: Peng, Ayala & Wang (arXiv:1906.05445) medem
`D/δx ≈ 37` para 1% no arrasto com o IBM do Uhlmann; Jiang & Liu
(arXiv:1806.09403) medem inclinação de convergência **≈1,0** para forçamento
direto simples contra ≈1,5 para MDF. Com inclinação 1, cada refino pela metade
corta o erro pela metade. Nossos 20/D estão bem abaixo de 1%.

### Um teste falsificável, sem mexer na formulação

A leitura da tabela da seção 9 é que há **dois erros somados**: a rampa linear é
o `O(Δt)` da defasagem, o **patamar** é o piso `O(h)` do `H·A ≠ I`.

Isso se testa refinando a malha **com `Δt` fixo** e vendo se o patamar cai com
`h`. Se cair: a rota do Uhlmann remove a rampa e **não** o patamar; só iterar
ataca o patamar. É o experimento de melhor retorno antes de mexer em formulação.

*(Esta decomposição é análise sobre a nossa medição, não citação: o levantamento
declarou não ter achado artigo que enuncie o escalonamento `O(Δt)` do resíduo do
forçamento defasado.)*

### Procedência a conferir

A referência do Uhlmann 2005 veio com um link de arXiv datado de 2018, o que é
inconsistente. O conteúdo citado é coerente com o artigo, mas o link precisa ser
conferido antes de entrar em bibliografia.
