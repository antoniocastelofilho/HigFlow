# Sondas de viabilidade — fronteira imersa

Duas perguntas, dois programas. Nenhum entra em build do projeto: a pergunta aqui
é **o que as bibliotecas já ligadas sabem fazer**, não como integrá-las.

```bash
set -a; . ../../../varsrc; set +a
make && make rodar
```

## O que se estava julgando

A malha lagrangeana de uma fronteira imersa precisa de duas coisas que a
`distributed_property` da HiGTree **não tem**, e não por descuido — ela foi feita
para outra coisa:

1. **Migração.** Um `distributed_property` é indexado por id local da partição
   euleriana, com o `gid_map` fixado no `psd_synced_mapper`. Elemento não muda de
   dono entre reconstruções. Marcador de corpo móvel atravessa fronteira de rank
   durante a simulação.

2. **Conectividade de superfície.** Em 2D a curva tem ordem natural e a
   vizinhança é implícita. Em 3D a superfície é triangulada, a área dual e a
   normal do vértice dependem dos triângulos incidentes, e **não existe ordem
   linear que torne essa vizinhança local**. É o 3D que decide a estrutura.

## O que cada sonda estabelece

### `sonda-dmswarm.c` — migração com carga

Quatro marcadores por rank, com campos `forca` (3 componentes) e `peso`, todos
mandados para o rank seguinte. Mede se chegam e se o conteúdo sobrevive.

```
np=2   rank 0: antes 4, depois 4, carga do rank 1 INTACTA
       rank 1: antes 4, depois 4, carga do rank 0 INTACTA
np=3   os três, INTACTA
```

Usa `DMSWARM_BASIC` e escreve o rank de destino à mão, no campo
`DMSwarmField_rank`. **Isso é de propósito**: a malha de fundo deste projeto não
é um `DM` do PETSc, então o `DMSWARM_PIC` — que localizaria o ponto sozinho —
não se aplica. A localização sai do instantâneo (cláusulas C7–C9 do contrato de
Mesh) e o destino é escrito por nós. A sonda mostra que esse caminho existe.

### `sonda-dmplex.c` — superfície triangulada distribuída

Um octaedro fechado: 8 triângulos, 6 vértices, topologia de dimensão **2** num
espaço de dimensão **3**. Criado inteiro no rank 0 e distribuído.

```
np=2   4 cel / 5 ver | 4 cel / 5 ver            suporte máx. de vértice 4
np=3   3 cel / 5 ver | 3 cel / 6 ver | 2 cel / 4 ver    suporte máx. 4, 4, 3
```

O número que importa é o **suporte do vértice**: os triângulos incidentes
continuam alcançáveis depois da distribuição. É o que área dual e normal exigem,
e é exatamente o ponto em que uma estrutura caseira ficaria cara.

### `sonda-petscsf.c` — acumulação franja→dono

O buraco que a estrutura distribuída reabre. O espalhamento do delta escreve em
facetas de **outro** rank, e o `dp_sync` da HiGTree (`pdomain.c:1308`) manda
dono→franja e **sobrescreve**: a contribuição na franja seria descartada, calada.

Cada rank escreve 1,0 nas entradas próprias e 10,0 numa entrada de franja do
vizinho. Depois do `PetscSFReduce` com `MPI_SUM`:

```
np=2 e np=3   entrada 0 = 11,0  (1 própria + 10 do vizinho)   ACUMULOU
              entrada 1 =  1,0
```

O teste discrimina: se o `PetscSF` sobrescrevesse em vez de somar, daria 10,0. O
grafo se monta do que a HiGTree já tem (`gid_map`, `local_count`, `total_count`
e os vizinhos filtrados de `_dp_shared`), sem tocar na estrutura dela.

### `sonda-plex-sobreposicao.c` — quanto halo cada nível de sobreposição custa

Esfera triangulada, refinada três vezes (1280 triângulos), distribuída com
sobreposição 0, 1 e 2. Fantasmas são as folhas do `pointSF`.

```
np=2   640 próprias    sobrep 0:    0 fantasmas       sobrep 1: +124 (+19%)    sobrep 2: +240 (+38%)
np=3   427 próprias    sobrep 0:    0 fantasmas       sobrep 1: +161..228 (+38..53%)   sobrep 2: +297..441 (+70..103%)
```

**A leitura que importa: sobreposição 1 basta, e sobreposição 0 não.** Já em
sobreposição 0 há vértices fantasmas (64 no rank 0, np=2) — os vértices da
fronteira da partição são compartilhados. Mas os **triângulos incidentes** a esse
vértice estão parte do outro lado, e área dual e normal dependem deles. Um anel
de células resolve; é o que sobreposição 1 entrega.

Sobreposição 2 só se justifica para quem precisa de dois anéis — curvatura por
estêncil mais largo, energia de flexão — o que é assunto de corpo **deformável**,
não do corpo rígido e fixo decidido.

**O custo cresce com o número de ranks**, e é preciso ver por quê: o halo é
proporcional ao perímetro da partição, e o perímetro relativo cresce quando a
superfície é dividida em mais pedaços. Numa esfera de 1280 triângulos em np=3, a
sobreposição 2 quase **dobra** a malha local. Para corpo pequeno partido em
muitos ranks, isso deixa de ser detalhe.

### Uma medida que não mediu nada, registrada de propósito

A primeira execução usou a esfera **sem refinar** — o `DMPlexCreateSphereMesh`
entrega um icosaedro de 20 triângulos. O resultado:

```
sobrep 1 | rank 0 | celulas 20 (10 fantasmas)
sobrep 2 | rank 0 | celulas 20 (10 fantasmas)
```

Parece a descoberta de que o halo satura no nível 1. **Não é.** Com 20 células e
sobreposição 1, cada rank já recebeu a malha inteira; o número não cresce porque
acabou a malha, não porque o halo parou. Medida cujo resultado é o teto do
domínio não mede o que se pensa estar medindo — e essa é a forma mais fácil de
errar aqui, porque o número sai bonito e estável.

## O que as sondas NÃO provam

- **Nada sobre integração.** Elas não tocam HiGTree nem HiGFlow. Que as
  bibliotecas saibam fazer isso não diz que o acoplamento seja barato.
- **Nada sobre o halo do delta regularizado, e a sobreposição do Plex NÃO é esse
  halo.** Conflacionei as duas coisas antes e corrijo aqui: a sobreposição é da
  malha LAGRANGEANA e serve a área e normal do vértice; o suporte do delta é
  EULERIANO, e quem o cobre é a franja da HiGTree mais a acumulação do `PetscSF`.
  A largura de franja euleriana que o delta exige segue sem medir.
- **Nada sobre desempenho.** Oito triângulos e quatro marcadores não medem nada.
- **Nada sobre a escolha** — que já foi tomada, e contra a minha recomendação.
  Eu recomendava replicar, pelo tamanho: ~160 marcadores para um cilindro nestas
  malhas. O argumento está certo para a primeira versão e errado para a
  trajetória: corpo móvel exige migração e corpo deformável exige topologia
  distribuída, e chegar neles com marcadores replicados custa trocar a estrutura
  com o método já funcionando em cima dela. Decidido usar as malhas do PETSc.

## Por que não uma biblioteca nova

`DMSwarm`, `DMPlex` e `PetscSF` vêm do PETSc 3.25.4, que já é dependência. O
Zoltan já é ligado em todo binário (`-ltrilinos_zoltan`, usado pelo `lbal.c`) e
traz `zoltan_dd.h` (diretório distribuído: id global → dono) e `zoltan_comm.h`
(comunicação irregular), que serviriam a uma estrutura própria.

libMesh, MOAB e deal.II resolveriam também, mas cada uma é dependência nova e
grande, com build próprio — e este projeto acabou de pagar esse custo com o
t8code. CGAL é serial e tem licenciamento a considerar. Com duas bibliotecas
capazes já ligadas, dependência nova precisa se justificar.
