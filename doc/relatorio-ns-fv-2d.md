# t8SolverNS-FV 2D — primeiro relatório

Navier–Stokes incompressível newtoniano, volumes finitos, arranjo deslocado
(MAC), sobre a malha de floresta do t8code. Este documento cobre o que existe,
o que foi medido, o que está recusado por decisão e o que está em aberto.

**Estado:** o canal de Poiseuille roda e converge; a divergência é zero de
máquina; `u`, `v` e `p` convergem com ordem medida. Serial apenas, malha
uniforme apenas — as duas limitações são recusadas explicitamente pelo código,
não descobertas pelo usuário.

---

## 1. As estratégias, e por que são as do HiGFlow

O objetivo declarado é comparar os dois sistemas **no mesmo caso**. Um método
diferente faria qualquer diferença medida ser inconclusiva: não se saberia se
veio da malha, que é o que se quer medir, ou do método, que não é.

| eixo | escolha | observação |
|---|---|---|
| acoplamento p–u | segregado, por projeção | não incremental e incremental |
| tempo | Euler explícito | semi-implícito desenhado, **não implementado** |
| convectivo | central, upwind 1ª, CUBISTA, QUICK | padrão CUBISTA, o do caso de referência |
| arranjo | MAC deslocado | escalares/tensores nos centros, vetores nas facetas |
| forma do convectivo | **divergente**, não advectiva | conserva quantidade de movimento por construção |
| entrada | YAML, esquema do HiGFlow | o mesmo arquivo roda nos dois |
| saída | VTK, um escritor por registro de campos | contra ~20 variantes no HiGFlow |

A forma divergente `∑ ∂(u_e u_d)/∂x_e` e a advectiva `∑ u_e ∂u_d/∂x_e`
coincidem no contínuo porque a divergência é nula, e **não** coincidem no
discreto: a divergente conserva quantidade de movimento por construção, a
advectiva só a conserva quando a divergência discreta já é nula — o que no meio
do passo ela não é.

### A condição do incremento não é a da pressão

No Poisson da projeção, `∇²δp = ∇·u*/dt`, a condição de contorno é **derivada**
da condição da pressão, não igual a ela:

- onde a pressão tem **Neumann**, o incremento tem Neumann homogêneo;
- onde a pressão tem **Dirichlet** com valor `g`, o incremento tem
  `δp_parede = g` na não incremental (ali `p = δp`) e `δp_parede = g − p^n_parede`
  na incremental.

Esse último ponto custou um nível de pressão inteiro, e está na seção 4.

---

## 2. O caso: canal de Poiseuille

Geometria e números do `example2d_Newt` do HiGFlow, para que a comparação seja o
mesmo caso e não dois parecidos: `x ∈ [0,8]`, `y ∈ [−1,1]`, Re = 1.

**Condições:** entrada com perfil parabólico de Dirichlet e pressão de Neumann;
saída com pressão Dirichlet nula e gradiente normal nulo nas duas componentes
da velocidade; duas paredes com velocidade nula e pressão de Neumann. As
condições são funções do usuário (`std::function` do ponto da parede), como no
HiGFlow.

**Solução analítica**, com a adimensionalização `∂u/∂t + ∇·(uu) = −∇p + Re⁻¹∇²u`:

```
u(y) = (3/2)(1 − y²),    v = 0,    dp/dx = −2 u_max / Re = −3
```

O domínio é coberto por **quatro blocos quadrados** de lado 2, e não por um
bloco 8×2: a folha do t8code é quadrada, e um bloco 8×2 com refino uniforme
daria célula de proporção 4:1 — a comparação deixaria de ser do mesmo caso.

### O que este caso NÃO verifica

O termo convectivo **se anula identicamente**: `u` não depende de `x` e `v ≡ 0`.
Central, upwind de primeira e CUBISTA dão o mesmo resultado aqui. Um caso que
não distingue os três não está exercitando nenhum, e até existir um que
distinga, *"CUBISTA está correto"* não é afirmação verificada.

---

## 3. Resultados

Euler explícito, `dt = 5·10⁻⁴`, critério de estacionariedade
`max|u^{n+1} − u^n|/dt < 10⁻¹¹`. Erro na norma do máximo contra a analítica.

| nível | h | células | `u` | `v` | `p` incremental | `p` não incr. | div. máxima |
|---|---|---|---|---|---|---|---|
| 3 | 0,25   | 256  | 1,4926e−02 | 6,1378e−03 | 6,6604e−01 | 1,0325e+00 | 6,1e−11 |
| 4 | 0,125  | 1024 | 4,7708e−03 | 1,8586e−03 | 2,2378e−01 | 4,1019e−01 | 3,9e−12 |
| 5 | 0,0625 | 4096 | 1,3280e−03 | 5,1731e−04 | 8,3419e−02 | 1,7703e−01 | 1,6e−13 |

**Ordens medidas:**

| grandeza | 3→4 | 4→5 |
|---|---|---|
| `u` | 1,65 | 1,85 |
| `v` | 1,72 | 1,85 |
| `p` incremental | 1,57 | 1,42 |
| `p` não incremental | 1,33 | 1,21 |

Três leituras, em ordem de importância.

**A divergência é zero de máquina** — 6,1e−11 a 1,6e−13, caindo com o refino. É
o oráculo mais barato e o único que não precisa de solução de referência: o
método de projeção *existe* para zerar isto. Ele está fazendo o que promete.

**A incremental é consistentemente melhor que a não incremental** na pressão,
por um fator de 1,5 a 2,1, e com ordem maior. É exatamente o que a teoria prevê
— a incremental é de segunda ordem na pressão e a não incremental de primeira —
e o par é um oráculo que nenhuma das duas produz sozinha: se as duas dessem o
mesmo, uma delas estaria errada.

**`v` não é zero, e isso é esperado.** A analítica tem `v ≡ 0`; o discreto dá
`v = O(h^1,85)`. A razão está na seção 5.

### Uma afirmação minha que a medida refutou

Escrevi, ao desenhar o caso, que o erro seria de **máquina** e não uma taxa: o
perfil é quadrático, o laplaciano de três pontos é exato para quadráticas, o
gradiente de pressão é constante e o convectivo se anula. A medida diz ordem
1,85, não exatidão. O raciocínio ignorava o **fechamento de parede** — ver a
seção 5. Fica registrado porque o argumento era convincente e estava errado.

---

## 4. O defeito que o caso encontrou: o nível da pressão

Com Dirichlet **homogêneo** no incremento — que é a escolha óbvia e é o que
muitos textos descrevem —, a pressão saía assim:

```
p calculada  [6,2988e+01, 8,5729e+01]
p exata      [3,7500e-01, 2,3625e+01]
```

Erro 62,6, **maior que a própria solução**. E o diagnóstico está na faixa, não
no erro: a faixa calculada é 22,74 contra 23,25 da exata — **o gradiente está
certo** (2%, compatível com segunda ordem) e o campo inteiro está deslocado por
uma constante.

**A causa.** No ponto fixo, `δp ≡ 0` resolve a equação *qualquer que seja o
nível de `p`*. Dirichlet homogêneo prende a **variação**, não a pressão. O
primeiro passo — partida impulsiva, campo de pressão da ordem de 4·10³ —
deixava `p` em 62,99 perto da saída, e o nível ficava lá para sempre.

A não incremental não sofria disso **porque recalcula `p` do zero a cada
passo**, e esse contraste foi o que localizou a causa: mesmo matriz, mesma
condição, comportamentos opostos.

**O conserto** é a condição correta, `δp_parede = g − p^n_parede`, com
`p^n_parede` por extrapolação **linear** das duas células mais próximas
(`p_parede = (3p_c − p_cc)/2`; a parede está a `h/2` da primeira e `3h/2` da
segunda). Extrapolar a ordem zero forçaria `p_c → g` e introduziria erro de
`h/2` vezes o gradiente — primeira ordem disfarçada de detalhe.

Depois: erro de **62,6 para 0,666**, 94×, e a incremental passou a ser melhor
que a não incremental, como deve.

---

## 5. Por que o canal não é exato: o fechamento de parede

A velocidade tangencial numa parede é fechada por **reflexão**,
`u_fantasma = 2u_parede − u_dentro`, que acerta o valor *na parede* e é de
segunda ordem. Mas ela **não é exata para quadráticas**.

Com `u = U(1−y²)`, parede em `y = −1`, primeira faceta a `h/2`:

| | valor em `y = −1 − h/2` |
|---|---|
| exato | `−U(h + h²/4)` |
| reflexão | `−U(h − h²/4)` |
| diferença | `U·h²/2` |

No laplaciano essa diferença vira `U/2`, isto é, o laplaciano discreto na
primeira fileira dá `−1,5U` em vez de `−2U`. O efeito na solução é `O(h²)`,
que é o que se mede — mas **a exatidão se perde na fileira junto à parede**.

E é daí que vem o `v ≠ 0`: o perfil interior que o esquema sustenta não é
exatamente o da entrada, de modo que há uma região de desenvolvimento, e numa
região de desenvolvimento `v` não é nulo.

**Um fechamento exato existe**: ajustar uma quadrática por
`(parede, 1ª faceta, 2ª faceta)` dá `u_fantasma = (8/3)u_parede − 2u_1 + (1/3)u_2`.
Com ele o canal passaria a ser exato. **Não foi adotado** porque muda o método
em relação ao HiGFlow, e a comparação é o objetivo declarado. Fica como opção a
oferecer, e como contraste de verificação.

---

## 6. Recusado por decisão, e dito pelo código

Nenhuma das duas produz número: as duas lançam com a razão.

**Paralelo.** Faltam a troca de franja dos campos — o estêncil lê célula e
faceta do vizinho, e os campos cobrem só os ids locais — e o percurso por
blocos, que hoje supõe todas as árvores locais. Com a franja ausente o estêncil
leria lixo perto da fronteira de partição e a corrida daria número.

**Malha adaptada.** A transição 2:1 numa faceta de *velocidade* não é o mesmo
problema que numa faceta de fluxo escalar: a componente normal se restringe por
conservação de vazão e a tangencial por interpolação, e aplicar a mesma regra às
duas produz divergência que não fecha — defeito que aparece como o solver de
pressão não convergindo, longe da causa.

Isso responde ao pedido de **refinar junto à parede**: hoje é recusado.

```
malha REFINADA junto as paredes: faixa de 0.3, nivel 4
celulas: 448, em transicao 2:1: 192
NavierStokes: a malha tem 192 folhas em transicao 2:1, e este solver ainda nao
a trata.  A restricao da componente normal e a da tangencial sao diferentes, e
aplicar a mesma as duas produz divergencia que nao fecha.
```

A opção `--refina-parede=E` existe no exemplo e hoje serve para **mostrar a
recusa**, não para medir. Medir o efeito do refino junto à parede exige a
restrição 2:1 das facetas, que é o próximo item de peso.

---

## 7. Em aberto

**O Poisson não converge em `dt = 5·10⁻⁴`, nível 4, incremental.** Para no passo
1991 e o solver **diz** — não produz número. Com `10⁻³`, `4·10⁻⁴`, `2·10⁻⁴` e
`10⁻⁴` converge, e os dois últimos dão o mesmo erro (2,2378e−01). Não é limite
de estabilidade: o limite difusivo do Euler explícito no nível 4 é
`Re·h²/4 ≈ 3,9·10⁻³`, e `5·10⁻⁴` está dez vezes abaixo. Hipótese não confirmada:
a tolerância do KSP é **relativa**, e com o lado direito indo a zero perto do
estacionário o critério relativo fica difícil de atingir. Pede piso absoluto.

**Euler semi-implícito** está no desenho e não implementado. Exige um sistema de
Helmholtz por componente, sobre as facetas.

**Métricas em YAML** — o canal de resultado, separado do VTK — ainda não foram
escritas. É o que permitirá comparar diretamente com `ci/reference/` do HiGFlow.

**Rebaixamentos do convectivo:** 112 no nível 3. Não é defeito — junto à parede
o ponto a montante do esquema de segunda ordem cai fora do domínio e não há o
que fazer —, mas é contado e impresso, porque um esquema que se anuncia de
segunda ordem e cai para primeira em parte do domínio entrega menos do que diz.

---

## 8. Sobre o VTK

O escritor é **um só**, dirigido por registro de campos: escalar no centro,
vetorial nas facetas (interpolado ao centro), tensorial, e fronteiras de
interface em arquivo próprio. Um arquivo por processo.

**VTK é visualização, não resultado**, e isto não é opinião: o próprio HiGFlow o
registra no cabeçalho do IO dele, depois de medir — `np=4` pode diferir no VTK
com o estado do solver bit a bit idêntico (velocidades concordando em 0,0 sobre
15.920 facetas enquanto métricas tiradas do VTK mudavam por um fator de 115). A
razão é a interpolação faceta→centro, que não é a mesma com e sem vizinho local.

Por isso os números deste relatório saem dos campos em memória, com precisão
cheia, e **não** do arquivo. O `id_global` vai em todo arquivo para permitir
casar uma célula do VTK com a mesma célula de um despejo de precisão cheia.
