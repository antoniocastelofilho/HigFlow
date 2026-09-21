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

## O que as sondas NÃO provam

- **Nada sobre integração.** Elas não tocam HiGTree nem HiGFlow. Que as
  bibliotecas saibam fazer isso não diz que o acoplamento seja barato.
- **Nada sobre o halo do delta regularizado.** O `DMPlexDistribute` foi chamado
  com sobreposição **zero**. O suporte do delta (3 a 4 células) exige sobreposição,
  e isso não foi medido.
- **Nada sobre desempenho.** Oito triângulos e quatro marcadores não medem nada.
- **Nada sobre a escolha.** Nas centenas de marcadores destas malhas (`h = 0,01`;
  um cilindro de diâmetro 0,05 dá ~160 marcadores), **replicar a malha inteira em
  todo rank ainda é mais simples** que distribuir — e faz desaparecer o problema
  de acumulação franja→dono, porque nenhum rank escreve fora do que possui. As
  bibliotecas passam a valer com muitos corpos, superfície bem mais fina que a
  malha, ou corpo deformável.

## Por que não uma biblioteca nova

`DMSwarm`, `DMPlex` e `PetscSF` vêm do PETSc 3.25.4, que já é dependência. O
Zoltan já é ligado em todo binário (`-ltrilinos_zoltan`, usado pelo `lbal.c`) e
traz `zoltan_dd.h` (diretório distribuído: id global → dono) e `zoltan_comm.h`
(comunicação irregular), que serviriam a uma estrutura própria.

libMesh, MOAB e deal.II resolveriam também, mas cada uma é dependência nova e
grande, com build próprio — e este projeto acabou de pagar esse custo com o
t8code. CGAL é serial e tem licenciamento a considerar. Com duas bibliotecas
capazes já ligadas, dependência nova precisa se justificar.
