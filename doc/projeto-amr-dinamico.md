# Refinamento adaptativo dinâmico

Documento vivo.  Estado em 24/09/2026: **F1, F2 e F3 concluídas**, com o
ciclo validado em DOIS backends (t8code e MTree, Cd a 0,0085% um do outro).
Com histerese a malha converge monotonicamente; sem, nem o estacionário
converge.  F4 (2D-2, Re=100): a maquina roda ate' o fim e a esteira despende, mas o Strouhal nao sai limpo -- remalha a cada 0,5 injeta picos de forca (11% da serie) e N=100 perturba o ciclo.  Proximo: congelar a malha apos o transiente e medir a cauda fixa.

## Por quê, e por quê agora

O refinamento **estático** resolveu o 2D-1: escoamento estacionário, esteira
parada, uma caixa bem posta basta — medido, o Cd é insensível ao comprimento
da caixa e responde ao refinamento como a malha uniforme, a 8,7× menos células.

Ele deixa de bastar em dois lugares:

- **2D-2 (Re = 100)**: a esteira de von Kármán oscila; nenhuma caixa fixa a
  cobre sem virar o domínio inteiro.
- **3D**: uma caixa grande o bastante custa o que se queria evitar.

E o motivo de fazer sobre o t8code, não sobre a máquina de AMR da HiGTree,
continua o registrado: a adaptação dinâmica da HiGTree está quebrada em np=2, e
o escalonamento paralelo dela é justamente o que motivou integrar o t8code.

## O que já existe, todo verificado

O ciclo de remalhamento precisa de cinco peças, e as cinco existem:

| peça | onde | verificação |
|---|---|---|
| produção de malha refinada | `t8_produz_por_rank_brick_refinado` | `test-refino-caixa-t8code`: sem lasca, 2:1, duas camadas, cobertura = t8code |
| transferência por posição | `hig-flow-remalha` | bit a bit, 2D e 3D, com guarda de vacuidade |
| interpolação conservativa | `rem_interpola_*` | integral exata em célula; constante exato em faceta |
| projeção pós-remalha | `higflow_projecao_remalha` | razão 1e-6 uniforme, 1,3e-3 refinada |
| pressão consistente na interface | montagem composta em `higflow_pressure` | div por fluxo ~1e-6; dourados bit a bit |

O que **não** existe: o critério, o ciclo que liga as peças, e a reconstrução
do domínio no meio da corrida.

## O ciclo

A cada N passos:

1. **Critério** avaliado na malha corrente → marca por célula: refinar,
   engrossar, manter.
2. **Expansão de duas camadas** das marcas de refino (a exigência do MLS — a
   mesma da escada estática) e **histerese**: refina acima de θ_alto, só
   engrossa abaixo de θ_baixo.  Sem histerese a malha bate à toa
   ("flapping") e o custo de remalhar domina.
3. **Marcas → floresta**.  A floresta do t8code vive entre ciclos (o padrão do
   `t8-forest-cache`); as marcas chegam a ela por **encontro marcado por
   posição** — a partição do t8code não é a do HiGFlow, e o mecanismo do
   `hig-flow-remalha` já resolve exatamente isso.
4. `adapt` + `balance` + `partition(set_for_coarsening=1)` + produção das
   caixas pelo produtor consertado.
5. **Reconstrução do domínio** — o núcleo de engenharia, ver abaixo.
6. **Colher → plantar → interpolar → sincronizar → projetar.**  A colheita vem
   **antes** de destruir os domínios velhos (ela não guarda ponteiro para eles
   — propriedade já verificada).
7. **Corpo imerso**: os marcadores são fixos e sobrevivem no DMSwarm; o suporte
   é reencontrado a cada interpolação.  O critério carrega uma **cláusula
   estática**: a vizinhança do corpo fica sempre no nível mais fino, com folga
   maior que o suporte do núcleo — e as guardas existentes
   (`fi_suporte_nivel_trocado == 0`) afirmam isso a cada quadro.

## A reconstrução do domínio

É a peça sem precedente no código: destruir e recriar, no meio da corrida,
`psd`/`psfd`/`sfbi`/mapeadores/propriedades/solvers/CCs, mantendo `ns->par` e o
estado temporal.  Uma função nova, `higflow_remalha(ns)`, feita das mesmas
chamadas do arranque (`higflow_create_domain`, `psfd_compute_sfbi`,
`higflow_initialize_boundaries_yaml`, `higflow_create_solver`) — reusar o
caminho do arranque, não duplicá-lo.

Riscos mapeados:
- as CCs precisam ser reinstanciadas na malha nova (o caminho yaml do arranque
  deve servir; verificar que é re-entrável);
- os solvers lineares mudam de tamanho (recriar, não redimensionar);
- a franja de 5 células vem do `partition_graph` novo;
- custo: a reconstrução tem de custar menos que N passos, senão N sobe.

## As decisões que são do usuário

1. **A física do critério.**  Candidatos: |ω| (vorticidade), ‖∇u‖, ou híbrido
   — zona estática do corpo + sensor de esteira.  E os limiares θ_alto/θ_baixo.
2. **A cadência N** (custo × acompanhamento da esteira).
3. **Engrossar ou não** fora da esteira, e o nível de piso.
4. **O caso graduado no conjunto dourado** — pendência herdada do conserto da
   interface: hoje nenhuma referência exercita a montagem composta.

## Fases, cada uma com o seu portão

**F1 — reconstrução a malha igual.  CONCLUÍDA.**  `higflow_reconstroi_dominio`
(hig-flow-kernel.c) destrói e recria domínios/propriedades/solvers/CCs no meio
da corrida, transfere u e p pela posição e interpola o que faltar.  Medido, com
o corpo imerso ativo (Uhlmann, 5 iterações), rebuild no passo 15 de 30:

- **Identidade no rebuild**: bit a bit em malha uniforme (np=1 e np=6).  Na
  malha t8-refinada os lids permutam — a produção reordena árvores na segunda
  chamada do mesmo processo — e o veredito é em dois níveis: nada sem valor, e
  somas invariantes a permutação iguais a 1e-14 relativo.  Lid é detalhe de
  implementação, não identidade da malha.
- **Continuação**: a diferença com/sem rebuild É a tolerância do solver linear,
  e escala com ela — sum(p) difere 2,5e-3 relativo com rtol=1e-5 e 8e-10 com
  rtol=1e-12 (sete ordens).  A reconstrução não introduz erro próprio.
- **Três armadilhas achadas e curadas no caminho**: o rascunho do adaptador
  Uhlmann era cache preso ao domínio antigo (curado com o *aviso de remalha*,
  gancho chamado depois de colher e antes de destruir — depois, nem o
  dp_destroy é seguro); os oráculos não podem ler lids fora do instantâneo (o
  solve implícito carrega lixo dependente de ambiente nos DOFs sem linha
  montada); e após mudança de layout do struct, a cadeia INTEIRA rebuilda.
- **Vazamento deliberado e documentado**: sd/sfd/árvores/CCs antigos não são
  liberados (posse compartilhada com o balanceador); custo por remalhamento, a
  medir na F2.

**F2 — remalha estático-equivalente.  CONCLUÍDA.**  Cd do ciclo com 9
reconstruções (N=200) = 5,536659 contra 5,536660 da corrida contínua —
diferença de 1e-6.  Custo: 0,9 s por reconstrução (0,4% do tempo total);
vazamento deliberado ~15 MB/ciclo no rank 0 (104 → 239 MB em 9) — apertar
antes de corridas longas de F4.  Nada sem valor em nenhuma das nove.

Como planejada: **F2 —**  O ciclo inteiro rodando no 2D-1, com um
critério que reproduz a caixa estática a cada ciclo: malha igual, campos
transferidos, Cd no fim igual ao da corrida estática.  *Portão: Cd da corrida
com N remalhamentos == Cd da corrida sem nenhum (referência: a remedição de
dois níveis com o conserto, em curso).  Medir o custo por remalhamento aqui.*

**F3 — critério dinâmico no 2D-1.  RODADA; PORTÃO DE CONVERGÊNCIA NÃO FECHOU.**
2000 passos, remalha a cada 200, critério híbrido (VORT1=1,5, VORT2=6, SEM
histerese).  O que passou: guardas zeradas nos 9 ciclos, resíduo do corpo
igual ao estático (1,73e-2), custo 5,2 s/remalha, toda posição preenchida.  O
que não passou: a malha DERIVA (67,9k → 71,1k → 64,9k células, nenhuma
pulada), |ω|max oscila em picos (24 → 60 → 25 → 41) num escoamento que deveria
assentar, Cd = 5,752 (pior que uniforme-20) e Cl = −0,516.  Assinatura de
**flapping**: células cruzando os limiares a cada ciclo, malha mudando,
remalha perturbando, realimentação.  É o dado que a decisão de histerese
(que é do usuário) pedia; um experimento com histerese está programado para
gerar a comparação.

Como planejada: **F3 —**  Escoamento estacionário: a malha deve
**convergir** para uma malha fixa e parar de mudar.  *Portão: número de células
estabiliza; Cd igual ao estático; custo de remalha → 0 depois da convergência.*

**F4 — 2D-2, Re = 100.**  O alvo.  A esteira oscila e a malha a segue.
*Portão: Strouhal e amplitudes de Cd/Cl dentro das faixas do benchmark de
Schäfer–Turek; a malha acompanha a esteira nos VTKs.*

**F5 — 3D.**  Depois do 2D-2, com a geometria do 3D-1Z.

A ordem é a lição da sessão: cada fase tem um oráculo exato antes de a próxima
introduzir física nova.  F1 e F2 não têm física — são identidade e
equivalência — e é neles que os defeitos de máquina aparecem baratos.
