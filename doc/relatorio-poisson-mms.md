# Poisson com soluções manufaturadas — verificação do t8Flow 2D

Relatório de **12/10/2026**. Primeiro exemplo do t8Flow, serial, e primeiro caso
da suíte que resolve uma equação.

    exemplo   t8flow/exemplos/poisson-mms.cxx   (parametrizado, sem recompilar)
    teste     t8flow/tests/2d/test-poisson2d.cxx   22 casos
    modulos   src/2d/poisson2d, manufaturada2d, contorno2d, campos2d, malha2d

## 1. O que o método faz, e onde está a força dele

Escolhe-se `u` conhecido, **deriva-se** o termo de fonte `f = ∇²u` e os dados de
contorno a partir dele, e mede-se o erro contra o `u` escolhido. A força está em
`f` e nos dados serem **exatos**, derivados à mão: se fossem calculados
numericamente, o oráculo mediria a consistência entre dois erros meus, não o
esquema.

Por isso cada solução do catálogo (`manufaturada2d.cxx`) carrega as suas derivadas
com a derivação escrita ao lado, e não há diferenciação automática em lugar nenhum.

## 2. Como o operador foi montado

**Balanço de fluxos**, que é o que o arranjo deslocado pede. Pela forma integral,
`∫_V ∇²u = Σ_faces ∂u/∂n · |face|`, então a linha de cada célula é a soma dos
fluxos pelas **suas** facetas igualada a `f(c)·área`.

Não é o laplaciano pontual montado no centro — e essa distinção não é estilo. A
cláusula C5 do contrato deixou medido que **no centro da célula o estêncil MLS é a
identidade**: peso 1,0 na própria célula, 6,4e-11 nas outras, `rhs` ~0. Montar
operador ali daria uma identidade, não um laplaciano.

**O fechamento na parede é de segunda ordem**, com três pontos — a parede (valor
`g`), o centro da célula (a `h/2`) e o centro do vizinho (a `3h/2`):

    du/ds(0) = -8g/(3h) + 3u_c/h - u_viz/(3h)

A soma dos coeficientes é `-8/3 + 3 - 1/3 = 0`, como tem de ser: a derivada de uma
constante é zero. Essa conferência pega erro de sinal e de escala.

**Neumann é o caso trivial aqui**, e vale dizer por quê: no balanço de fluxos, o
que a condição dá é exatamente `∂u/∂n` — o próprio fluxo. Nenhuma aproximação,
nenhum coeficiente novo na linha. É o lado em que o arranjo deslocado paga o que a
forma pontual cobraria com célula fantasma.

## 3. Os resultados

Domínio `[0,1]²`, malha uniforme, solver PETSc com tolerância 1e-14. As colunas de
ordem são `log₂(e_h / e_{h/2})`.

### 3.1 Todas as paredes Dirichlet

| solução | erro inf (nível 5) | erro L2 | ordem inf / L2 |
|---|---|---|---|
| **grau 2** | **3,91e-13** | **1,69e-13** | **exato** |
| grau 3 | 2,22e-05 | 1,19e-05 | 2,95 / 2,99 |
| trigonométrica | 3,78e-03 | 1,53e-03 | 1,97 / 1,98 |

**O grau 2 dá erro de máquina nas duas normas** — o que o pedido pedia. E o
resíduo do solver na mesma corrida é 4,4e-15: o erro que sobra é dele, não do
esquema.

### 3.2 Duas paredes Dirichlet, duas Neumann

| solução | erro inf | erro L2 | ordem |
|---|---|---|---|
| **grau 2** | **1,07e-12** | **7,21e-13** | **exato** |

### 3.3 Uma parede Dirichlet, três Neumann

É o caso de contorno mais exigente que o problema admite: com Neumann nas
**quatro** ele é singular — a solução fica a menos de uma constante, o solver não
tem o que fixar, e medir erro contra a exata mediria a constante.

| solução | erro inf | erro L2 | ordem inf / L2 |
|---|---|---|---|
| **grau 2** | **3,38e-12** | **2,09e-12** | **exato** |
| grau 3 | 5,67e-04 | 2,67e-04 | 1,91 / 1,90 |
| trigonométrica | 4,95e-03 | 2,60e-03 | 1,98 / 1,98 |

### 3.4 O contraste do fechamento: grau 2 com dois pontos

| nível | h | erro inf | ordem | erro L2 | ordem |
|---|---|---|---|---|---|
| 3 | 1,25e-01 | 1,57e-02 | — | 8,28e-03 | — |
| 4 | 6,25e-02 | 4,42e-03 | 1,82 | 2,22e-03 | 1,90 |
| 5 | 3,13e-02 | 1,16e-03 | 1,93 | 5,66e-04 | 1,97 |

**Dez ordens de grandeza** acima do fechamento de três pontos (1,16e-03 contra
3,91e-13), convergindo em ~2. É o fechamento clássico de volumes finitos: correto,
de segunda ordem global, e **não** exato.

## 4. Por que o grau 2 é exato, e o grau 3 não

A explicação é contra-intuitiva e vale guardar, porque ela diz exatamente **até
onde** a exatidão vale.

O dado de Neumann é o fluxo exato **no centro** da faceta, mas o balanço de
volumes finitos precisa do fluxo **médio** sobre ela. A quadratura do ponto médio
é exata para integrando **linear**:

| | `∂u/∂n` ao longo da faceta | resultado |
|---|---|---|
| grau 2 | varia **linearmente** | ponto médio = média → **exato** |
| grau 3 | varia quadraticamente | erro O(h²) → ordem 2 |

O mesmo vale para a fonte: `f(c)·área` aproxima `∫f dV`, e é exato enquanto `f`
for linear — para grau 2, `f` é constante.

Então a exatidão para o grau 2 **não é sorte nem coincidência de coeficientes**: é
a quadratura do ponto médio sendo exata nos dois lugares, no fluxo e na fonte. E é
por isso que o grau 3 cai para ordem 2 com Neumann (5,67e-04, ordem 1,91) mas
chega a ordem ~3 com Dirichlet em tudo (2,22e-05, ordem 2,95): o caminho de
Dirichlet não passa pela quadratura do fluxo dado.

## 5. O que a suíte afirma

13 casos em `test-poisson2d.cxx`. Os que valem destacar não são os de exatidão —
são os que **impedem o verde de significar menos do que parece**:

- **o grau 2 com dois pontos tem de errar muito.** Se ele também desse erro de
  máquina, o caso de exatidão não distinguiria esquema exato de sistema montado
  como identidade — que também daria 1e-13;
- **a exatidão não pode depender do refino.** Se o erro caísse como h², seria
  convergência e não exatidão, e o caso estaria passando por a malha ser fina. O
  caso exige que a razão entre os níveis 3 e 5 fique entre 0,25 e 4;
- **a ordem só é medida onde há erro para medir.** Com esquema exato, `log₂` de
  ruído sobre ruído é um número aleatório que pareceria uma taxa. A tabela imprime
  "exato" e o teste exige erro acima de 1e-6 antes de calcular ordem;
- **o que a montagem não sabe, ela recusa:** transição 2:1 (na face grossa o fluxo
  é a soma dos dois lados finos), Robin (mistura valor e derivada na mesma faceta),
  fonte não instalada, contorno vazio, nome de solução inexistente — este último
  com a lista dos válidos na mensagem, porque errar o nome e receber a solução
  padrão em silêncio seria medir uma coisa e reportar outra.

## 6. Quadtree adaptada: a transição 2:1

Acrescentado em **12/10/2026**, com malha gerada por refino aleatório
reprodutível (40% sorteado, 106 de 136 folhas em transição).

### 6.1 O esquema

Para cada vizinho, **a mesma fórmula** da faceta uniforme: a face compartilhada
tem o tamanho do menor dos dois, e a distância é a soma dos meios-passos. Uma
fórmula e não três (uniforme, grossa→fina, fina→grossa), porque três divergem — e
divergir ali quebra a **conservação** sem quebrar nada visível. A igualdade
conferida:

    grossa->fina:  larg_C/2 / (3h_C/4)  =  (2/3) larg_C/h_C
    fina->grossa:  larg_F   / (3h_F/2)  =  (2/3) larg_C/h_C

Logo o que sai de uma entra na outra com o mesmo valor, e a matriz fica simétrica.

**Mais a correção transversal, que não é opcional.** Na transição os centros não
estão alinhados com a normal — o centro da fina está deslocado um quarto de célula
—, e para a célula **fina** esse erro não cancela, porque ela vê uma única face
contra a grossa. A correção desloca o valor do lado grosso para a linha normal do
centro da sub-face, usando o gradiente transversal.

### 6.2 Os resultados

O oráculo principal aqui **não** é o erro da solução: é o **resíduo do operador**,
`max |A·u_exata − b| / área`, sem passar pelo solver. O erro da solução mistura
consistência, condicionamento e tolerância, e quando é grande não se sabe qual dos
três.

| | malha uniforme | quadtree adaptada |
|---|---|---|
| resíduo do operador, grau 1 | 1,14e-13 | **4,55e-13** |
| resíduo do operador, grau 1 + 3 Neumann | — | **1,82e-12** |
| resíduo do operador, grau 2 | 2,27e-13 | **4,50e+00** |

| solução (adaptada) | erro inf | ordem inf | erro L2 | ordem L2 |
|---|---|---|---|---|
| grau 1 | 6,7e-16 → 3,9e-12 | **exato** | idem | **exato** |
| grau 2 | 1,63e-02 → 2,96e-04 | 2,08 / 1,78 / 1,92 | | 2,16 / 1,94 / 1,98 |
| trigonométrica | 4,80e-02 → 1,09e-03 | 1,63 / 2,11 / 1,73 | | 2,07 / 2,10 / 1,80 |

**O operador é exato para campo linear** na malha adaptada, inclusive com três
paredes de Neumann. **O grau 2 deixa de ser exato** — a exatidão era propriedade da
malha uniforme, exatamente como a análise da quadratura previa: na transição, a
distância entre os centros (3h/4) não é a distância do centro à face (h/2), e o
termo de curvatura não cancela. A convergência fica em ~2.

### 6.3 Os três defeitos, e o que os achou

Nenhum era de formulação — todos de sinal ou de caminho esquecido. Guardo os três
porque o que os encontrou foi sempre o mesmo oráculo:

1. **a correção transversal com o sinal dependendo de quem era o grosso.** O erro
   **dobrou** em vez de cair — 7,8e-02 contra 3,9e-02 sem correção nenhuma. Erro que
   dobra é a assinatura de uma correção invertida, e foi o que levou a refazer a
   álgebra em vez de ajustar números. A conta certa não depende de quem é o grosso;
2. **`s2` medido do centro em vez da parede** no fechamento genérico. Em malha
   uniforme deu `h` onde o certo é `3h/2`, e o caso de regressão do grau 2 pegou na
   hora: 1,97e-01 em vez de 3e-13;
3. **com parede de Neumann no lado transversal, a correção era abandonada** — e o
   resíduo voltava a **53,3**, exatamente o valor de quando não havia correção
   nenhuma. O número repetido foi o que apontou a causa: não era outro defeito, era
   o **mesmo** reaparecendo por um caminho que eu havia deixado aberto. Ali o dado
   de Neumann *é* `∂u/∂t`, e usá-lo é exato — com `∂u/∂t = sgn·g`, não `g`.

E a guarda do contrato pegou a mudança sozinha: a cláusula P4 ficou apontando para
o caso de recusa de 2:1, que deixou de existir, e saiu como **CASO AUSENTE**.

## 7. Diferenças compactas de quarta ordem (malha uniforme)

Duas variantes do **mesmo** esquema de nove pontos, com lados direitos diferentes:

| | lado direito | precisa de |
|---|---|---|
| **Mehrstellen** | `f + (h²/12)∇²f` | `∇²f` analítico |
| **Lele** (implícito na segunda derivada) | `A_x A_y f` — `f` nos nove pontos | só `f` |

O esquema de Lele parte de `(1/10)u''ᵢ₋₁ + u''ᵢ + (1/10)u''ᵢ₊₁ = (6/5)δ²u/h²`. Aplicando
o operador tridiagonal nas duas direções a `∇²u = f`, o lado esquerdo dá os mesmos
nove pontos e o direito fica `A_x A_y f`.

**E os dois são o mesmo esquema**, conferido na álgebra: `A_x f/(6/5)` expande em
`f + (h²/12)f_xx + O(h⁴)`, porque `(1/10)/(6/5) = 1/12`. Os estênceis em `u`
diferem pelo fator `(6/5)² = 1,44` — medido contra Mehrstellen: `−4,8/−3,333 =
1,44` e `0,96/0,667 = 1,44`.

### 7.1 Os resultados

| solução | erro inf (n5) | ordem inf | ordem L2 |
|---|---|---|---|
| grau 1, 2, 3 | 1e-14 a 4e-13 | **exato** | **exato** |
| grau 4 | 2,59e-06 | 3,95 | 3,98 |
| grau 6 | 1,22e-05 | 3,86 | 3,92 |
| trigonométrica | 5,44e-05 | 3,86 | 3,93 |

E no mesmo nível, contra o esquema de segunda ordem: **3,78e-03 → 5,44e-05, 70×**.

### 7.2 O limite é o fechamento, não o interior

O estêncil de nove pontos tem erro proporcional a derivadas de **sexta** ordem,
logo seria exato até grau 5. Mas o fechamento de parede usa quatro pontos — a
parede e três células — e é exato até grau **3**. O corte medido é nítido:

    grau1 1,44e-14   grau2 1,12e-13   grau3 3,58e-13   grau4 3,99e-05

A ordem de convergência não sofre com isso (o erro do fechamento também é O(h⁴));
o que se perde é a **exatidão** em grau 4 e 5. Subir o fechamento para cinco
pontos é trabalho isolado — a fórmula de Lagrange é a mesma com um nó a mais.

### 7.3 A verificação entre as duas variantes

Ter as duas é o que permite verificar a álgebra em vez de confiar nela:

- em polinômio de grau ≤ 4 dão valores **idênticos** (ali o termo de O(h⁴) do lado
  direito envolve derivadas quartas de `f`, que são nulas);
- na trigonométrica diferem **1,8%**, que é a ordem do próprio truncamento.

E há valor prático na diferença: `∇²f` existe na manufatura e **não** existe em
problema real. O esquema de Lele usa apenas `f`.

### 7.4 O que o compacto recusa

Malha adaptada (o estêncil usa vizinhos **diagonais**, que na transição 2:1 não
estão definidos sem reconstrução própria), mais de um bloco, Neumann (o fechamento
é de Dirichlet; usá-lo para Neumann daria quarta ordem no interior e segunda na
parede, sem aviso), e Mehrstellen **sem** `∇²f` — este porque rodar com zero daria
segunda ordem e pareceria quarta falhando, com o número saindo plausível.

## 8. Derivada MLS na faceta e quadratura de Gauss

O esquema de fluxo de segunda ordem perde a exatidão do grau 2 na malha adaptada
(seção 6), e o compacto de quarta ordem não serve na malha adaptada por construção
(seção 7.4). O esquema `mlsgauss` fecha esse vão: é o único do sistema que passa de
segunda ordem **em malha adaptada**.

### 8.1 As duas peças, e por que as duas

O fluxo pela faceta é `∫ ∂u/∂n ds`. O esquema anterior aproximava os dois fatores em
primeira/segunda ordem, e levantar um só não passa de segunda ordem:

| peça | antes | agora |
|---|---|---|
| derivada no ponto | diferença de dois centros | ajuste **MLS de ordem 3 ou 4** |
| quadratura da faceta | ponto médio (exata p/ linear) | **Gauss 2 pontos** (exata p/ cúbico) |
| quadratura do volume | `f(c)·A` (exata p/ linear) | **Gauss 2×2** (exata p/ bicúbico) |

A terceira linha é o achado que eu não esperava, e está medido na seção 8.3.

### 8.2 Conservação, por construção e não por acerto

O suporte do estêncil depende apenas do **ponto** de quadratura, não da célula que
pergunta. Logo as duas células que compartilham uma faceta obtêm o mesmo fluxo, sem
ninguém combinar nada — inclusive na transição 2:1, onde a célula grossa parte a
faceta em duas sub-facetas e a fina não parte a sua.

A cláusula D7 afirma isso pelo caminho onde pode falhar: o oráculo refaz a geometria
das sub-facetas por conta própria e compara os estênceis peso a peso. Oito
sub-facetas de transição, diferença **0,0e+00** no ponto e no peso. Isso importa
porque um esquema **não conservativo ainda converge** — nenhum dos testes de
convergência pegaria um erro de indexação na sub-faceta.

### 8.3 Os resultados

Malha uniforme, todas as paredes Dirichlet:

| solução | fluxo2 | mlsgauss ordem 3 | mlsgauss ordem 4 |
|---|---|---|---|
| grau 2 | exato | exato | exato |
| grau 3 | 1,71e-04 | **exato** (7,99e-15) | exato |
| grau 4 | ordem 2 | 5,50e-06, ordem 3,7 | **exato** (5,91e-14) |
| trig | ordem 2 | 3,30e-04, ordem 2,8 | 1,11e-04, **ordem 3,8** |

Malha adaptada (`--aleatorio=0.4`), que é a razão de existir do esquema:

| solução | fluxo2 | mlsgauss ordem 3 |
|---|---|---|
| grau 1 | exato | exato |
| grau 2 | 3,86e-03 (ordem 2) | **exato** (2,59e-13) |
| grau 3 | ordem 2 | **exato** (4,86e-14) |
| grau 4 | ordem 2 | 4,08e-06, ordem 3,6 |
| trig | ordem 2 | 1,76e-04, **ordem 3,1** |

**O grau de exatidão é a ordem do ajuste MLS, em malha uniforme e em adaptada.** É a
propriedade que o esquema de fluxo não tem.

### 8.4 A quadratura da fonte era o limite, não o fluxo

Com o fluxo da faceta já exato e `∫f dA` ainda na regra do ponto médio, o esquema
**travava em segunda ordem**. Grau 4 com ajuste de ordem 4, malha uniforme — a única
diferença entre as duas colunas é a quadratura do volume:

| nível | ponto médio | Gauss 2×2 |
|---|---|---|
| 3 | 1,49e-03 | 1,51e-14 |
| 4 | 3,81e-04 (ordem 1,97) | 2,52e-14 |
| 5 | 9,58e-05 (ordem 1,99) | 5,91e-14 |

A regra do ponto médio é um erro `O(h²)` **relativo** em `∫f dA`; no esquema de fluxo
de segunda ordem ele se escondia atrás do erro do próprio fluxo. O interruptor
`--fonte-ponto-medio` existe para esta medida: verde que nunca foi vermelho não prova
nada.

### 8.5 A ordem 4 em malha adaptada rebaixa, e agora isso aparece

Medido: em malha adaptada com ordem 4, **8 dos 5172** ajustes não sustentam a ordem.
A causa é o suporte, não o condicionamento: o crescimento por anéis entrega **12 a 16
pontos**, contra 15 monômios que a ordem 4 pede.

O que o código fazia antes era pior que o problema: `resolve_pequeno` só recusava
pivô abaixo de `1e-300` — a singularidade exata — e devolvia um **vetor de zeros**.
O efeito era uma linha da matriz silenciosamente errada:

| nível | erro inf (antes do conserto) |
|---|---|
| 3 | 2,36e-10 |
| 4 | **3,32e-01** |
| 5 | 2,35e-10 |

Um nível com erro de 3,3e-01 cercado por 1e-10 nos vizinhos. Dois consertos:

1. **a base é adimensionalizada** pelo raio do suporte. Sem isso a matriz normal da
   ordem 4 tem condicionamento `O(h⁻⁸)` — com `h = 3e-2` isso passa de `1e12`;
2. **o rebaixamento de ordem é registrado.** Quando o suporte não sustenta a ordem
   pedida, o ajuste cai para a que ele sustenta e `rebaixamentos()` conta; o exemplo
   imprime um AVISO. Zero é o que se quer poder afirmar, e as cláusulas D4–D6 o
   afirmam explicitamente.

Com o rebaixamento, o erro de 3,32e-01 cai para 1,25e-10 — mas a **exatidão do grau 4
se perde** na malha adaptada, e deve se perder: oito pontos de quadratura usam ordem
menor. **A ordem utilizável em malha adaptada é 3**, e levá-la a 4 é um trabalho de
seleção de suporte — mais pontos escolhidos por qualidade geométrica, não por anéis —
não de quadratura nem de fechamento.

## 9. Limites declarados

- **a exatidão do grau 2 vale só na malha uniforme** *no esquema de fluxo*. Na
  adaptada ele é de segunda ordem, exato apenas para campo linear; o `mlsgauss`
  recupera a exatidão até o grau 3 (seção 8);
- **a ordem 4 do `mlsgauss` não sobrevive à malha adaptada.** Rebaixa em 8 de 5172
  ajustes, por falta de pontos no suporte, e o rebaixamento é contado e avisado;
- **Robin não está montado.** Dirichlet e Neumann sim. Robin mistura valor e
  derivada na mesma faceta: a linha ganha coeficiente que depende de α, β e h, e
  montá-lo como se fosse Neumann daria resultado plausível e errado;
- **serial.** `np=1`. Em paralelo entram a franja e a repartição, que são os itens
  seguintes da lista;
- **o piso do erro é o solver.** Pedir menos de 1e-13 seria exigir do solver linear
  algo que ele não promete com tolerância 1e-14.

## 10. Como reproduzir

```bash
cd t8flow && make 2d
mpirun -n 1 ./exemplos/poisson-mms --solucao=grau2 --contorno=dirichlet --niveis=3,4,5
mpirun -n 1 ./exemplos/poisson-mms --solucao=trig --aleatorio=0.4 --niveis=3,4,5,6
mpirun -n 1 ./exemplos/poisson-mms --esquema=mlsgauss --mls-ordem=4 --niveis=3,4,5 \
       --solucao=grau4                          # exato: derivada MLS + Gauss
mpirun -n 1 ./exemplos/poisson-mms --esquema=mlsgauss --aleatorio=0.4 \
       --solucao=grau3 --niveis=3,4,5            # exato em malha ADAPTADA
mpirun -n 1 ./exemplos/poisson-mms --ajuda      # lista as solucoes do catalogo
```

Parâmetros: `--solucao` (grau1, grau2, grau3, trig), `--contorno` (dirichlet,
misto, umdirichlet), `--fechamento` (3 ou 2), `--niveis`, `--tol`,
`--aleatorio=FRAÇÃO` (quadtree adaptada), `--semente`, `--esquema` (fluxo2,
compacto4, lele4, mlsgauss), `--mls-ordem`, `--mls-pontos`,
`--fonte-ponto-medio`. A configuração
inteira vai impressa no alto da saída — relatório sem a configuração é número sem
procedência, e neste projeto um build sem o que se queria medir já produziu duas
corridas "estáveis" que eram outra coisa.
