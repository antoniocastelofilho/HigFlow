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
np=3 numa esfera de 1280 triângulos). A largura de franja euleriana exigida pelo
delta segue **sem medir**. Com replicação isso é inofensivo na interpolação (a
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
