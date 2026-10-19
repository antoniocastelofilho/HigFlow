# Escoamento bifásico por front-tracking, 2D — relatório

Interface rastreada por poligonal de marcadores, acoplada ao
`t8SolverNS-FV 2D` por núcleo de Peskin, sobre a malha de floresta do t8code
com refino adaptativo e remalhamento. Este documento cobre o que foi medido,
o que foi refutado, e o que está em aberto.

**Estado:** a escada B1→B3 está completa. A bolha ascendente de Hysing, caso 1,
**acerta a velocidade máxima dentro da faixa das três referências publicadas**,
em malha adaptativa com menos da metade das células que a uniforme equivalente
precisaria.

---

## 1. A escada, e por que ela tem três degraus

Cada degrau tem um oráculo que **não depende do degrau seguinte**, e é isso que
o torna degrau e não etapa:

| | caso | o oráculo | depende de |
|---|---|---|---|
| **B1** | frente advectada em campo **prescrito** | forma fechada: curvatura de círculo e elipse, reversão de LeVeque | nada — não há solver |
| **B2** | gota **estática** | salto de Laplace `Δp = σ/R`, forma fechada | B1 + o solver, com velocidade ~0 |
| **B3** | bolha ascendente de Hysing | três grandezas do *benchmark*, três referências independentes | tudo |

Se B2 falhasse, o salto de Laplace não diria se o erro é da curvatura, do
núcleo ou da projeção — por isso ele vem depois de B1, onde a curvatura já foi
verificada sozinha.

---

## 2. B1 — a frente em campo prescrito

**A curvatura é exata para o círculo em qualquer resolução.** Ela sai do
circunraio de três marcadores, e o circunraio de três pontos sobre um círculo
*é* o raio. Não é "converge": é identidade. A tolerância, porém, tem de seguir
o condicionamento, que degrada como `1/θ²` com marcadores quase colineares.

**A reversão de LeVeque não converge com o passo quando há cirurgia.** Medido:

| `dt` | sem cirurgia | com cirurgia |
|---|---|---|
| 4e−03 | 7,56e−07 | 1,01e−04 |
| 2e−03 | 9,45e−08 | 1,01e−04 |
| 1e−03 | 1,18e−08 | 1,01e−04 |

Sem cirurgia a razão é **8,0** por bissecção — terceira ordem, porque a reversão
é simétrica no tempo e o termo de segunda ordem cancela. Com cirurgia há um
**piso** independente de `dt`: partir e fundir arestas é irreversível, e o piso
é o preço de manter resolução. Não é defeito; é o custo, medido.

**O que a cirurgia promete não é um fator, é uma cota.** A primeira versão deste
caso afirmava um contraste num instante — "a maior aresta fica 3× menor" — e
isso era um número escolhido, porque o fator depende de quando se olha. O que
ela promete é que a maior aresta fica **limitada** enquanto sem ela cresce
**monotonicamente**. Medido em múltiplos de `L₀`, nos quartos de `T/2`:

```
sem cirurgia   2,02   3,05   3,61   3,77      (cresce)
com cirurgia   1,50   1,49   1,50   1,49      (limitada)
```

A cirurgia **recusa** `l_min ≥ l_max/2`: a frente oscilaria entre dois estados a
cada passo, sem nada acusar.

---

## 3. B2 — a gota estática e o salto de Laplace

Fração volumétrica por **recorte exato** de polígono (Sutherland–Hodgman contra
a caixa da célula): um algoritmo só, sem casos dentro/fora. Força de tensão
distribuída pelo núcleo de Peskin de quatro pontos; a velocidade dos marcadores
vem da **adjunta** da distribuição — mesmo núcleo —, que é o que impede a frente
de ganhar ou perder energia na travessia.

**Medido**, com o termo viscoso completo (ver §4 do relatório do NS):

| nível | erro do salto | razão | correntes espúrias |
|---|---|---|---|
| 4 | 1,024e−04 | — | 5,570e−04 |
| 5 | 2,809e−05 | 3,65 | 3,659e−04 |
| 6 | 1,297e−05 | **2,17** | 2,446e−04 |

**A convergência do salto degrada.** A razão 3,65 é quase segunda ordem; a 2,17
é primeira. Relatórios anteriores deste projeto disseram "segunda ordem" — era
otimismo de olhar só o primeiro par. O que se pode afirmar é que converge, e que
a taxa cai.

**As correntes espúrias não vão a zero em taxa útil**: razões 1,52 e 1,50,
compatíveis com `h^0,6`. São artefato do método, e o exemplo as **imprime** em
vez de as omitir.

**A guarda de vacuidade pagou aqui.** O caso `R = 0,20` no nível 4 falhou com
*"0 células puras dentro"*: a banda excluída da medida tem 2,5 h de cada lado
(largura do núcleo), sobra um disco de raio `R − 2,5h`, e com menos de 1,5 h não
há centro de célula. A condição é `R > 4h`. Sem a guarda, o teste teria
comparado 0,0 com 5,0 e mandado procurar defeito no método. **Corrigi o caso,
não a guarda.**

O núcleo usa o `h` **da célula que contém o marcador** — numa malha adaptada não
há "o" `h`, e o errado espalha a força por largura errada: muda a curvatura
efetiva e desloca o salto sem nada acusar.

---

## 4. B3 — a bolha ascendente de Hysing

Domínio `[0,1]×[0,2]`, `R = 0,25` em `(0,5; 0,5)`, `ρ` 1000/100, `μ` 10/1,
`g = 0,98`, `σ = 24,5` — `Re = 35`, `Eo = 10`.

Malha **igual à do HiGFlow**: base 40×80 (`h = 0,025`) com dois níveis de refino
em volta da frente, `h_min = 0,00625` — quarenta células por raio. Avanço
semi-implícito, projeção incremental, CUBISTA, um marcador por célula fina
(`ds/h = 1`), TSUR a cada 20 passos.

### O resultado

| grandeza | t8Flow | HiGFlow FT | referência |
|---|---|---|---|
| **velocidade máxima** | **0,2418** | 0,241577 | **0,2417 – 0,2421** |
| instante do máximo | 0,9153 | — | ~0,92 |
| circularidade mínima | 0,8983 | 0,900344 | 0,9012 ± 0,0001 |
| instante do mínimo | 1,9053 | — | ~1,90 |
| centro de massa em `t=3` | 1,0766 | 1,078587 | 1,0799 – 1,0817 |
| deriva de massa | **1,3e−03** | 2,44e−03 | 0 |

A velocidade máxima cai **dentro** da faixa, e os dois instantes batem. A
circularidade e o centro de massa erram ~0,3%, contra ~0,1% do HiGFlow; a deriva
de massa é **metade** da dele. 23714 células ao fim, 385 remalhamentos.

**A discordância ENTRE as três referências publicadas é o limite de resolução
deste oráculo.** Bater melhor que ela não seria conquista — seria sinal de que
algo foi ajustado para dar o número.

### As três grandezas saem da fração volumétrica

Isso não é detalhe de implementação: é o que tornará a comparação com o VOF
possível. O VOF **advecta** esse campo; o front-tracking o **escreve** da
geometria da frente. A medição é a mesma nos dois, e qualquer diferença que
sobreviva pertence ao método, não à forma de medir.

---

## 5. O episódio da banda do AMR

É o trecho mais instrutivo deste relatório, e o que mais se aproveita noutros
casos.

### O sintoma: o AMR anticonvergia

Par de controle — **mesma `h_min`, mesmo `dt`, mesmos 2839 passos**:

| | células | `c_min` | `v_max` | `y_c(3)` |
|---|---|---|---|---|
| uniforme `h = 0,0125` | 12800 | 0,8963 | **0,2431** | 1,0764 |
| AMR `h_min = 0,0125` | 3620 | 0,8935 | 0,2262 | 1,0646 |
| AMR `h_min = 0,00625` | 8852 | 0,8966 | **0,2163** | 1,0648 |

O uniforme batia o *benchmark* (+0,4% em `v_max`). O AMR errava 7% **na mesma
resolução de interface**, e **piorava ao refinar**. Isso não é convergência
lenta — é anticonvergência, e aponta defeito e não malha grossa.

**O par de controle existe para isto.** Sem ele, a comparação seria contra a
referência publicada, e 7% ficaria indistinguível de "a malha ainda é grossa".

### A hipótese errada, e por que medi-la valeu

O erro escalava com o **número de remalhamentos** (136 → −6,5%, 391 → −10,5%),
não com `h`. E `c_min` e área quase não se moviam, enquanto `v_max` e `y_c` sim
— a interface é transferida **exata** (a frente não vive na malha), mas o
escoamento é interpolado. Suspeita natural: a transferência dissipa.

Contraste, mesma malha, variando só a frequência de remalhamento:

| `--amr-cada` | remalhamentos | `v_max` | energia perdida |
|---|---|---|---|
| 5 | **411** | 0,2254 | −7,79e−03 |
| 20 | **136** | 0,2254 | −7,77e−03 |
| 100 | **28** | 0,2254 | −7,62e−03 |

Quinze vezes mais remalhamentos, resultado **idêntico em quatro algarismos**. E
a perda é **invariante** — o salto por remalhamento escala ao inverso da
contagem, porque o que se perde é o desencontro acumulado *entre* remalhamentos,
não um custo fixo de cada um. Além disso 0,78% é pequeno demais para explicar
7%.

**Hipótese refutada.** Se eu tivesse "consertado" a transferência, teria mexido
em código são e o defeito continuaria.

### A causa: largura em contagem de células

`amr_limiares` derivava `banda_fina = cel_min·h_fino + margem` — tudo
proporcional a `h`. Logo a largura **física** da região resolvida encolhe com o
refino: **0,50 R** em `h_min = 0,0125` e **0,25 R** em 0,00625. A camada em
volta de uma bolha a `Re = 35` tem espessura própria, que não diminui porque a
malha refinou.

**Era resultado conhecido** — registrado na seção companheira — e este código
fazia o contrário. Ter anotado não substituiu conferir o código contra a
anotação.

### A correção, e a largura medida

`banda_fina = max(largura_física, cel_min·h_fino) + margem`. Duas exigências de
naturezas diferentes: a contagem de células é da **discretização** e escala com
`h`; a largura é do **problema** e não escala.

Quanto o escoamento pede é **medição, não escolha**. Varredura em
`h_min = 0,0125`, contra o uniforme na mesma resolução:

| `L` | banda | células | `v_max` | `y_c(3)` |
|---|---|---|---|---|
| 0,05 = 0,2 R | 0,50 R | 3704 | 0,2262 | 1,0646 |
| 0,125 | 0,75 R | 4886 | 0,2395 | 1,0694 |
| **0,25 = R** | 1,25 R | 6545 | **0,2429** | 1,0751 |
| uniforme | — | 12800 | 0,2431 | 1,0764 |

Monótona. Em `L = R` o adaptativo reproduz o uniforme a **0,08% com metade das
células**, que é o que o AMR tinha de entregar e não entregava. O padrão do
exemplo passou a ser `R`.

Eu havia posto `0,2 R` primeiro, por corresponder às oito células finas que o
HiGFlow resolveu na corrida dos números publicados. **Contar células da corrida
dele não substituiu medir a largura que o escoamento pede**, e as duas não
coincidiram.

### E a energia confirmou, depois

Com a banda larga a perda na transferência cai de −7,8e−03 para **−4,4e−04**:
ela era sintoma, não causa — como o contraste já havia indicado.

---

## 6. A simetria, e o medidor que faltava

O caso é simétrico em `x = 0,5`: dado inicial, contorno e gravidade. Qualquer
assimetria é numérica.

**A semente era a paridade dos marcadores.** Eles ficam em `θ_k = 2πk/n`, e a
reflexão leva `k → n/2 − k`: com `n` **ímpar** esse índice não existe e o
conjunto inicial não é simétrico. O nível 5 dava `nmarc = 201`. Medido: o centro
de área saiu de 0 e chegou a **+2,0e−02** em `t = 3` — 8% do raio, crescendo
monotonicamente. Com contagem par: **−8e−08**.

**Só apareceu porque alguém olhou uma animação.** Nenhum medidor dizia nada. Hoje
`x−0,5` é coluna da tabela de acompanhamento, ao lado de CFL, divergência,
iterações do Poisson, contagem de marcadores e área — todos impressos **durante**
a corrida, porque uma corrida de milhares de passos que degenera no meio tem de
ser vista *enquanto* degenera.

Os remalhamentos também são contados e ditos: conta **zero** com `--amr` ligado
significa malha estática com a banda no lugar errado, que produz resultado
plausível.

---

## 7. O ponto de salvamento

A corrida casada leva ~9,5 h, e esta máquina já drenou a bateria na tomada com um
solver num núcleo. Sem retomada, uma interrupção custa a corrida inteira.

**Ele é pequeno por desenho:** a malha é **derivada** da frente, então não se
grava malha nenhuma — bastam os marcadores, `u`, `p`, o relógio, o passo e os
contadores. 162 KB para 6572 células.

**A única sutileza real foi a tabela de refino.** A malha *em uso* não é
`tabela_da_frente(frente)` no instante do salvamento: o remalhamento só acontece
quando a tabela **muda**, logo a malha é a da última mudança e a frente andou
desde então. Recalcular a tabela na retomada daria malha diferente, e `u` e `p`
— gravados por índice — cairiam nos lugares errados **sem nada acusar**.

Critério de aceitação: reproduzir a corrida direta **bit a bit**, com
remalhamentos no meio. Salvando no passo 450 (não múltiplo de `--amr-cada=20`):

```
                    direta        retomada
medida em t=0,6     0,218515      0,218515      (e as iteracoes do Poisson: 108)
remalhamentos       22            22
energia acumulada   -5,6260e-05   -5,6260e-05
```

As iterações do Poisson são o indicador mais sensível: o KSP só gasta a mesma
contagem partindo do mesmo palpite com a mesma matriz. E a energia acumulada só
fecha se cada remalhamento tiver ocorrido no mesmo instante com o mesmo campo —
é ela que afirma que **gravar a tabela era necessário**.

---

## 8. Em aberto

**Os 0,3% residuais.** `c_min` e `y_c` erram ~0,3% contra ~0,1% do HiGFlow, **na
mesma direção e magnitude** — o que sugere causa única. Não há hipótese
confirmada. O candidato que eu investigaria primeiro é a densidade de marcadores
combinada com o TSUR, que praticamente não atuou nesta corrida: 4 marcadores
movidos em 8030 passos, contra 30 nas corridas mais grossas.

**O TSUR com `rlx = 0` não remove ondulação em 2D**, ao contrário do padrão do
HiGFlow: cada marcador é tocado por duas arestas e os deslocamentos acumulam, de
modo que `rlx = 0` passa 2× e o zigue-zague troca de sinal mantendo amplitude. O
fator de amortecimento é `|2·rlx − 1|`, e **0,5 é a conta em 2D**. Está medido e
é divergência deliberada em relação ao HiGFlow.

**Serial.** O caso bifásico recusa `np > 1` explicitamente.

**A comparação com o VOF** espera o VOF, cujo núcleo geométrico está feito e
verificado (ver o contrato, cláusulas V1–V7).
