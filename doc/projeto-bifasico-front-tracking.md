# Escoamento bifásico com front-tracking

Documento vivo. Estado em 24/09/2026: **planejado; nenhuma fase começou.**
Seção separada de propósito, para não sobrecarregar a fronteira imersa rígida
nem o AMR dinâmico — mas construída sobre a maquinaria dos dois.

## A fronteira que não se cruza

Front-tracking é **alternativa ao VOF**, não um complemento dele, e os dois
ficam separados no código. VOF captura a interface por um campo de fração
volumétrica advectado (euleriano); front-tracking a **rastreia** por uma malha
lagrangeana de marcadores conectados que se movem com o fluido. Não se mistura
a reconstrução geométrica de um com a advecção de fração do outro. Essa decisão
é do início da fronteira imersa e continua valendo.

## O que já existe, e reusa quase inteiro

O front-tracking é, estruturalmente, a fronteira imersa com três mudanças (os
marcadores andam, as propriedades saltam, a força é de tensão superficial). As
peças construídas nesta sessão servem quase todas:

| peça | onde | como serve ao front-tracking |
|---|---|---|
| DMSwarm de marcadores | `hig-flow-fronteira-imersa` | a frente É um swarm; já sobrevive à repartição |
| kernel de Roma (espalha/interpola) | idem | espalha a força singular de tensão superficial; interpola u nos marcadores para advectar |
| campo `indice` por marcador | idem | guardado JÁ pensando nisto — a curvatura precisa da conectividade (ver o comentário em `hig-flow-fronteira-imersa.c:167`) |
| transferência por posição | `hig-flow-remalha` | move u/p quando a malha euleriana readapta |
| AMR dinâmico + critério | `hig-flow-kernel` + exemplo | refina onde a ação está: **na interface** — o critério híbrido troca "vorticidade" por "distância à frente" |
| projeção consistente na interface 2:1 | `hig-flow-step` | a pressão salta na interface do fluido; o par montado=aplicado importa ainda mais aqui |
| VTK lagrangeano por quadro | exemplo (`HIGFLOW_VTK_LAGR`) | visualizar a frente que se move |

O que **não** existe, e é o trabalho novo: a advecção da frente, o salto de
propriedade, a força de tensão superficial, e a cirurgia da frente.

## As três mudanças estruturais

### 1. Os marcadores advectam

O coração do método. Onde a fronteira imersa rígida **calcula uma força** para
impor não escorregamento a marcadores fixos, o front-tracking **move** os
marcadores:

    X^{n+1} = X^n + dt · u(X^n)

`fi_interpola` já entrega u(X) — a mesma interpolação. Em vez de alimentar
`fi_forca_corpo_rigido`, alimenta a atualização da posição no DMSwarm. É reúso
direto, com o sinal trocado: interpola para mover, não para resistir.

### 2. As propriedades saltam

ρ e μ mudam através da interface — não há análogo na fronteira imersa rígida
(fluido único). Precisa de uma **função indicadora** I(x) (0 de um lado, 1 do
outro) reconstruída da posição da frente. O caminho de Tryggvason: espalhar o
gradiente da indicadora dos marcadores para a malha (o mesmo kernel de Roma) e
resolver um Poisson ∇²I = ∇·(espalhado) — reúso do solver de pressão que já
existe. ρ(x) e μ(x) então interpolam entre as duas fases por I(x).

**Consequência no solver:** os termos viscoso e de massa deixam de ter
coeficiente constante. É a mudança mais invasiva no passo de tempo, e a razão de
começar newtoniano-bifásico antes de qualquer reologia.

### 3. A força é de tensão superficial

    F = σ · κ · n · δ(x − X)

espalhada pelo kernel de Roma — a mesma máquina que espalha a força do corpo
rígido, com um integrando diferente. O que é novo: κ (curvatura) e n (normal)
vêm da **geometria da frente**, isto é, da conectividade dos marcadores. Em 2D a
frente é uma curva (marcadores em sequência, κ da segunda derivada); em 3D é uma
superfície triangulada (κ do operador de Laplace-Beltrami). O campo `indice` é o
que permite reconstruir essa sequência/malha depois da repartição.

## O ponto de projeto mais duro: a conectividade sob repartição

Na fronteira imersa rígida os marcadores são **desconexos** — cada um acha seu
suporte sozinho, e a repartição não dói. A frente do front-tracking é
**conectada**: em 2D uma polilinha, em 3D uma malha de triângulos. A curvatura
precisa dos vizinhos, e os vizinhos podem estar em outro rank depois do `lbal`.
O DMSwarm guarda pontos, não arestas.

Duas saídas, e é decisão de projeto (ver abaixo):
- **frente replicada**: a topologia da frente (índices e conectividade) vive
  inteira em todos os ranks, e só as posições se distribuem — barato enquanto a
  frente for pequena perto do volume, que é o caso de gotas/bolhas isoladas;
- **frente distribuída com franja**: a frente também se particiona, com uma
  franja de marcadores vizinhos — escala para muitas interfaces, ao custo de um
  protocolo de franja lagrangeano (o análogo do consenso de franja euleriano
  que esta sessão construiu).

Começar replicado é o caminho do primeiro resultado; o distribuído é o que casa
com a motivação de escalabilidade do t8code, para depois.

## A cirurgia da frente

Quando a interface estica, os marcadores se afastam e a curva perde resolução;
quando comprime, se amontoam. O front-tracking re-semeia: insere marcadores onde
o espaçamento passa de um limiar, remove onde fica abaixo. É maquinaria nova,
mas local e barata, e o critério é o mesmo Δs ≈ h que o corpo rígido já usa.

**Mudança de topologia** (gotas coalescendo ou se partindo) é o calcanhar do
front-tracking — ao contrário do VOF, não é automática. Fica **fora do escopo
inicial**: as primeiras fases assumem topologia fixa (uma gota, uma bolha).

## As decisões que são do usuário

1. **Representação da frente**: replicada (primeiro resultado rápido) ou
   distribuída com franja (escala, mas pede protocolo novo).
2. **Reconstrução da indicadora/salto**: Poisson do gradiente espalhado
   (Tryggvason) ou distância assinada à frente.
3. **Formulação**: um fluido com propriedade variável (I(x) interpolando ρ,μ) —
   o caminho padrão — versus dois domínios acoplados.
4. **Escopo de topologia**: confirmar que a versão inicial não trata
   coalescência/quebra.

## Fases, cada uma com o seu portão exato

**B1 — advecção cinemática, campo dado.** A frente advecta num campo de
velocidade **prescrito** (vórtice único, deformação reversível de Rider–Kothe).
Sem acoplamento, sem tensão superficial. *Portão: no teste reversível, a frente
volta à forma inicial em t=T; o erro de área é a medida.* Isola a advecção +
cirurgia de toda a física.

**B2 — gota estática (lei de Laplace).** Gota em repouso, só tensão superficial.
O salto de pressão tem de ser σ/R (2D) ou 2σ/R (3D). *Portão: o salto bate com
Laplace, e as correntes parasitas (spurious currents) — a velocidade espúria que
todo método de tensão superficial gera — ficam abaixo de um teto medido.* É o
teste que valida força + curvatura + salto de propriedade juntos.

**B3 — oscilação de gota.** Gota perturbada oscila na frequência de Lamb.
*Portão: a frequência bate com a teoria; o amortecimento numérico é medido.*
Valida a dinâmica acoplada.

**B4 — bolha subindo (benchmark de Hysing 2009, 2D).** O caso de referência com
número, contra o qual três grupos publicaram. *Portão: forma terminal,
velocidade de subida e circularidade dentro das faixas do benchmark.* Aqui o AMR
dinâmico entra de verdade — a interface se move pelo domínio e a malha a segue,
exatamente o que a F4 exercitou para a esteira.

**B5 — 3D.** Superfície triangulada, Laplace-Beltrami, bolha 3D. Depois de B4.

A ordem é a lição das seções anteriores: cada fase tem um oráculo exato — forma
reversível, lei de Laplace, frequência de Lamb, benchmark publicado — antes de a
próxima introduzir física nova. B1 e B2 não têm dinâmica acoplada; é neles que os
defeitos de máquina (advecção, cirurgia, curvatura) aparecem baratos.

## O que herda de graça, e o que não

Herda: o DMSwarm, o kernel, a transferência por posição, o AMR, a projeção
consistente, o VTK lagrangeano. Herda também as **armadilhas mapeadas**: o
vazamento de árvore do remalhamento (bloqueio real, ver
`higflow-vazamento-remalha`), a posse meio-aberta do marcador na fronteira de
partição, e o consenso de franja euleriano.

Não herda: tudo que depende de a densidade ser uniforme e de os marcadores serem
fixos. O passo de tempo com coeficiente variável (mudança 2) é o que mais se
afasta do que existe, e é o risco técnico central.
