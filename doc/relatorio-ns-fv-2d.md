# t8SolverNS-FV 2D — relatório

Navier–Stokes incompressível newtoniano, volumes finitos, arranjo deslocado
(MAC), sobre a malha de floresta do t8code. Este documento cobre o que existe,
o que foi medido, o que está recusado por decisão e o que está em aberto.

**Estado:** canal de Poiseuille e caixa manufaturada convergem com ordem
medida; a divergência é zero de máquina. Malha **adaptada** com transição 2:1 e
**remalhamento** dinâmico. Avanço explícito **e** semi-implícito. O termo
viscoso é o **tensor completo**. Serial apenas — a limitação é recusada
explicitamente pelo código, não descoberta pelo usuário.

O solver sustenta o caso bifásico por front-tracking, que tem relatório próprio
em [`relatorio-bifasico-ft-2d.md`](relatorio-bifasico-ft-2d.md) e acerta a
bolha de Hysing dentro da faixa das referências publicadas.

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

**Paralelo.** Lança com a razão, e não produz número. Faltam a troca de franja
dos campos — o estêncil lê célula e faceta do vizinho, e os campos cobrem só os
ids locais — e o percurso por blocos, que hoje supõe todas as árvores locais.
Com a franja ausente o estêncil leria lixo perto da fronteira de partição e a
corrida daria número.

**Malha adaptada — recusa levantada.** Este relatório registrava a transição
2:1 como recusada, com a razão: numa faceta de *velocidade* a componente normal
se restringe por conservação de vazão e a tangencial por interpolação, e
aplicar a mesma regra às duas produz divergência que não fecha — defeito que
aparece como o solver de pressão não convergindo, longe da causa.

A restrição foi feita (§6e) e a recusa saiu da construção. O que fica do
registro é o **modo de falha**, que continua valendo como aviso: divergência
que não fecha se manifesta num sintoma distante da causa.

A opção `--refina-parede=E` do canal, que existia para mostrar a recusa, hoje
mede: **refinar junto à parede dá a precisão da malha nível 4 com 44% das
células**, o que confirma que o erro do canal é governado pelo fechamento de
parede (§5).

---

## 6b. A caixa com tampa deslizante, por solução manufaturada

O caso que o canal não dava: aqui as duas componentes variam nas duas direções
e **o termo convectivo é exercitado**.

### A solução, e por que é essa

Construída por **função de corrente**, `u = ∂ψ/∂y`, `v = −∂ψ/∂x`, de modo que a
divergência é nula *por construção* e não por cancelamento que alguém conferiu
uma vez. Numa manufaturada incompressível essa é a propriedade que não pode
falhar: se `div u ≠ 0`, a projeção tenta corrigir o que é a própria solução.

```
ψ(x,y,t) = f(t) G(x) H(y),   G = sin²(πx),   H = sin(πy) − ½sin(2πy)
f(t) = sin(2πt)     (f(0) = 0: o repouso é condição inicial EXATA)
p(x,y) = cos(πx) cos(πy)
```

**Trigonométrica, não polinomial**, de propósito: um polinômio de grau baixo
pode ser reproduzido *exatamente* por um esquema de segunda ordem, e aí a ordem
medida é a da exatidão acidental, não a do esquema.

Cada raiz de `G` e `H` tem um propósito — `G(0)=G(1)=0` anula `u` nas laterais,
`G'(0)=G'(1)=0` anula `v` nelas, `H(0)=H(1)=0` anula `v` no fundo e na tampa,
`H'(0)=0` anula `u` no fundo, e `H'(1) = −2π ≠ 0` deixa `u` na tampa. Logo
`u(x,1,t) = −2π f(t) sin²(πx)` e **só a tampa depende do tempo**.

A pressão também não é livre: o método de projeção impõe `∂δp/∂n = 0` nas
paredes de Neumann, e uma pressão manufaturada com derivada normal não nula ali
seria incompatível com o próprio esquema. `cos(πx)cos(πy)` tem derivada normal
nula nos quatro lados e média zero.

### O Poisson fica singular, e o caso resolve com um ponto

Com velocidade de Dirichlet nas quatro paredes, a condição natural da pressão é
Neumann em todas, e a constante fica no núcleo do operador. O módulo **recusa**
esse caso em vez de resolver — o PETSc com matriz singular converge para um
elemento qualquer do subespaço, e isso dá número.

A saída é Dirichlet em **um ponto de parede**, onde a manufaturada dá o valor. A
região é `x < h` no fundo, para ser exatamente uma faceta em qualquer nível;
uma faixa de largura fixa cobriria mais facetas a cada refino e o caso deixaria
de ser o mesmo.

### Resultados — CUBISTA, `dt = h²/8`, `T = 1/16`

| nível | h | `u` | `v` | `p` | div. máx. | compat. |
|---|---|---|---|---|---|---|
| 3 | 0,125  | 9,9661e−02 | 4,8182e−02 | 1,6469e+00 | 1,3e−13 | 3,2e−14 |
| 4 | 0,0625 | 2,7059e−02 | 1,2173e−02 | 4,4805e−01 | 2,1e−13 | — |
| 5 | 0,03125| 7,0255e−03 | 3,0479e−03 | 1,1381e−01 | 4,8e−13 | — |

| grandeza | 3→4 | 4→5 |
|---|---|---|
| `u` | 1,88 | 1,95 |
| `v` | 1,99 | 2,00 |
| `p` | 1,88 | 1,98 |

**Segunda ordem nas três**, incluindo a pressão. A divergência é zero de
máquina e a compatibilidade — a soma dos fluxos, que numa caixa fechada tem de
ser nula — dá 3,2e−14.

### O caso discrimina, e o canal não discriminava

Nível 4, `T = 1/16`, variando só o esquema convectivo:

| esquema | `u` | `v` | `p` |
|---|---|---|---|
| central  | 2,7282e−02 | 1,1275e−02 | 4,2065e−01 |
| **upwind1** | **5,4162e−02** | **6,3651e−02** | **2,3733e+00** |
| cubista  | 2,7059e−02 | 1,2173e−02 | 4,4805e−01 |
| quick    | 2,7648e−02 | 1,2165e−02 | 4,4223e−01 |

O de primeira ordem é 2× pior em `u`, 5× em `v` e 5× em `p`; os três de segunda
concordam entre si. **É a verificação de que o termo convectivo está sendo
exercitado** — no canal os quatro davam o mesmo número.

### Duas coisas que a medida ensinou

**O upwind de primeira ordem estagna, e `p` até cresce com o refino.** Em
`T = 1/16`: `u` 1,05e−01 → 5,42e−02 → 4,64e−02 (ordens 0,95 e 0,22) e `p` 2,38
→ 2,37 → 2,65. Erro que cresce com refino é bandeira vermelha, e a causa foi
medida: é **acúmulo no tempo**. Em `T = 1/64` o mesmo upwind converge
normalmente (`u` 1,66e−02 → 5,58e−03 → 1,76e−03). A difusão numérica de
primeira ordem amortece a solução progressivamente, e o ganho do refino é
compensado pelo erro acumulado.

**O tempo de simulação decide se o caso discrimina.** Em `T = 1/64` upwind e
CUBISTA dão praticamente o mesmo (1,66e−02 contra 1,69e−02 em `u`); em
`T = 1/16` separam por 6× em `u` e 23× em `p`. Uma corrida curta teria
concluído *"todos os esquemas concordam"* — exatamente a cegueira que o canal
tinha, por outro motivo.

### E uma armadilha que virou recusa

`T` tem de ser múltiplo inteiro de `dt`. Com `dt ~ h²`, um `T` qualquer dá
número de passos fracionário e **cada nível para num instante diferente**; a
solução exata depende do tempo, de modo que a comparação mistura erro de
discretização com diferença de instante, e a ordem medida sai errada sem que
nada acuse. Medido: com `T = 0,05` o nível 3 parava em 0,050781. O exemplo
agora recusa e diz qual `T` serve.

---

## 6c. O termo viscoso: o tensor completo, e o patamar que o denunciou

Este relatório dizia, numa versão anterior, que o termo viscoso era
`div(μ∇u)` **e que essa era a forma do HiGFlow**. As duas afirmações estavam
erradas, e a segunda era erro de leitura do código dele. Abrindo,

```
div(μ ∇uᵀ)ᵢ = (∂ⱼμ)(∂ᵢuⱼ) + μ ∂ᵢ(∂ⱼuⱼ)
```

o HiGFlow omite só o **segundo**, que é `μ ∇(∇·u)` e some num campo solenoidal.
O primeiro está lá, em `hig-flow-discret.c`, somado em
`higflow_extra_diffusive_term` e dividido por `Re` e por `ρ`.

Descartar o transposto **inteiro** joga fora a tensão cisalhante através da
interface — e com razão 10 de viscosidade, `∇μ` ali é da ordem de `9/h`.

**O que custou**, na bolha de Hysing caso 1, nível 5:

| | sem o transposto | com | referência |
|---|---|---|---|
| circularidade mínima | 0,8361 | **0,8852** | 0,9012 ± 0,0001 |
| velocidade máxima | +5,2% | **+1,6%** | — |
| área / área inicial | 1,0249 | **1,0183** | 1 |

A circularidade errava **650 vezes a dispersão entre as três referências
publicadas**.

**O que denunciou não foi o valor, foi o patamar sob refino.** Nível 5 dava
+5,2% em `v_max`; nível 6 ainda +4,3%. Segunda ordem teria dividido por quatro.
*Patamar sob refino é assinatura de termo faltando, não de malha grossa* — e foi
isso que separou "resolução" de "defeito" sem precisar de mais corridas.

Com `μ` uniforme o termo é **identicamente nulo**, de modo que canal e caixa
manufaturada reproduzem bit a bit o que davam antes.

---

## 6d. Euler semi-implícito

`(I − (dt/ρ) L_μ) u* = uⁿ + dt R`, uma componente por vez. Implícito **só**
`div(μ∇u)`, como no HiGFlow. Convectivo, gravidade, tensão, fonte e o
**transposto** ficam explícitos — este último porque depende da derivada
cruzada `∂_d u_e`, e pô-lo na matriz acoplaria as duas componentes num sistema
só.

**O contorno entra na matriz, não no lado direito.** O fechamento é afim na
própria incógnita (`A u_f + B`): Dirichlet dá `A = −1`, Neumann `A = +1`. Jogar
o fantasma inteiro à direita trata um Neumann como Dirichlet com valor velho —
converge, dá número, e a condição passa a valer outra coisa. Nos lados de
deslizamento livre da bolha de Hysing a condição tangencial **é** Neumann.

**O ganho neste caso é zero, e isso está dito no código.** Em `h = 0,00625` o
limite viscoso é 9,77e−4 e o capilar 9,34e−4: o capilar já mandava, e
implicitar a difusão não o remove — é onda capilar, não difusão. O
semi-implícito entra por **correspondência** com o HiGFlow. Em `h` mais fino a
ordem se inverte (viscoso cai com `h²`, capilar com `h^{3/2}`). Custo medido:
16,4 s/passo contra 17,1 do explícito — os dois sistemas de velocidade são
baratos perto do Poisson.

**Duas implementações do operador, de propósito.** O implícito precisa dele como
matriz, o explícito o avalia como número. Colapsá-las trocaria a ordem das somas
de ponto flutuante, e os casos monofásicos deixariam de reproduzir bit a bit —
critério que já pegou duas refatorações mentindo. A defesa é o teste
(`difusao_explicita` contra `difusao_por_pesos`), como em `amr_tabela` contra
`amr_tabela_exaustiva`.

**É serial, e isso está dito e não suposto:** `DominioFacetas` numera as facetas
localmente, sem deslocamento global, então o índice da faceta é a linha da
matriz só em `np = 1`.

---

## 6e. Malha adaptada: transição 2:1 e remalhamento

**A posse da faceta dependia da orientação do refino.** A regra antiga — "menos
é sempre da célula, mais só no contorno" — dava *duas* facetas finas numa
orientação e *uma grossa* na espelhada, esta última com medida grossa cobrindo
duas faces finas. Regra nova: **a faceta pertence ao lado mais fino**, nas duas
orientações.

O invariante que faltava: para cada célula e cada lado, a soma das medidas das
facetas que cobrem a face é o lado transversal da célula. Mata os dois defeitos
que passam pela posse única — face sem faceta dá zero, faceta grossa dá o dobro.

**O critério de aceitação foi reprodução bit a bit dos casos uniformes**, não
"os números parecem razoáveis": o caminho geral reduz ao particular, logo os
uniformes verificam o código adaptado. Duas vezes isso pegou erro — a
velocidade na quina vinha de uma faceta só (caixa +20%/+37%, e o canal mal se
mexeu, porque nele o convectivo se anula), e uma troca de esquema do ponto a
montante disfarçada de refatoração.

**O remalhamento é um "reset".** Tudo o que depende da malha vive num bloco que
se joga fora e se refaz — malha, facetas, campos, contornos e solver. Refaz-se
**da base**, e não adapta no lugar, porque `refina_por_criterio` só sabe pedir
mais nível: como a frente se move, a banda fina que ela deixa para trás tem de
voltar a ser grossa.

No caso bifásico só `p` e `u` atravessam; `frac`, `ρ`, `μ` e a força são
reescritos da geometria da frente. A velocidade transferida **não** é
discretamente solenoidal, e não precisa de passe próprio: o Poisson do passo
seguinte a limpa.

**Armadilha:** `dt` sai do `h` **mínimo**, não o da base — refinar um pedaço
aperta o passo por `4^níveis` em todo o domínio, e errar isso diverge calado,
com NaN centenas de passos depois. E o pino de Dirichlet da pressão usa o `h`
**da parede**, não o mínimo: a região tem de conter exatamente um centro de
faceta de parede, e `x < h` com o `h` errado pega zero ou dois.

**A base multibloco é parâmetro** porque `h = 1/2^nível` só alcança potências de
dois, e o `h_min` do HiGFlow na bolha é `1/160`. `--base=5 --nivel=3 --amr=2`
reproduz a malha dele. Oráculo do parâmetro: mesma `h` por decomposição de
blocos diferente concorda em 4e−06.

**A banda de refino em contagem de células faz o estudo anticonvergir** — o
episódio inteiro está em [`relatorio-bifasico-ft-2d.md`](relatorio-bifasico-ft-2d.md), §5,
porque foi lá que apareceu.

---

## 6f. A adjunção do par divergência/gradiente

Divergência e gradiente têm casos próprios — exatidão em campo linear,
fechamento da projeção —, e os dois podem passar com o **par** desencontrado. A
adjunção é o que diz se a troca discreta de energia entre velocidade e pressão
se conserva, e nenhum teste de um operador só a vê:

```
<Du, p>_V = -<u, Gp>_W ,   com   W_f = medida_f x dist_f
```

O peso não é escolha: sai de igualar os coeficientes de `u_f p_i` nas duas
formas bilineares, e com outro peso a igualdade não sai nem no caso uniforme.

| malha | desvio relativo |
|---|---|
| uniforme | **4,5e−17** |
| adaptada, nível 4 | 8,47e−06 |
| nível 5 | 1,78e−06 |
| nível 6 | 4,18e−07 |

Razões **4,75** e **4,26**: segunda ordem. Na adaptada não vale exatamente, e
era previsível — o gradiente leva correção transversal na transição 2:1 e a
divergência não tem termo que a case. Sem a correção o Poisson tinha ordem
**negativa**, então ela não é negociável.

**O que importa é a ordem, não o valor.** Um valor pequeno num nível só não
distinguiria "perturbação que some" de "defeito pequeno que fica" — por isso a
asserção é a razão entre níveis.

No método segregado isso é **invisível**: o Poisson da projeção é montado
literalmente como `D G`, consistente por construção. A medida existe porque a
propriedade caracteriza a discretização, e porque num sistema acoplado seria
dela que sairia a escolha do Krylov (GMRES e não MINRES, com o complemento de
Schur montado como se simétrico, válido até `O(h²)`). **O solver acoplado foi
planejado e posto de lado**; a medida sobreviveu ao plano.

---

## 7. Em aberto

**O Poisson não converge em `dt = 5·10⁻⁴`, nível 4, incremental.** Para no passo
1991 e o solver **diz** — não produz número. Com `10⁻³`, `4·10⁻⁴`, `2·10⁻⁴` e
`10⁻⁴` converge, e os dois últimos dão o mesmo erro (2,2378e−01). Não é limite
de estabilidade: o limite difusivo do Euler explícito no nível 4 é
`Re·h²/4 ≈ 3,9·10⁻³`, e `5·10⁻⁴` está dez vezes abaixo. Hipótese não confirmada:
a tolerância do KSP é **relativa**, e com o lado direito indo a zero perto do
estacionário o critério relativo fica difícil de atingir. Pede piso absoluto.

**Paralelo.** O NS é serial. `DominioFacetas` numera as facetas localmente, sem
deslocamento global — o semi-implícito precisa de numeração global antes de
rodar em `np > 1`, e os campos precisam de troca de franja. O lado do Poisson
3D já é paralelo e verificado; é por lá que o escalonamento começa.

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
