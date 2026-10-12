# Poisson com soluções manufaturadas — verificação do t8Flow 2D

Relatório de **12/10/2026**. Primeiro exemplo do t8Flow, serial, e primeiro caso
da suíte que resolve uma equação.

    exemplo   t8flow/exemplos/poisson-mms.cxx   (parametrizado, sem recompilar)
    teste     t8flow/tests/2d/test-poisson2d.cxx   13 casos
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

## 6. Limites declarados

- **malha uniforme.** Com transição 2:1 uma faceta passa a ter dois vizinhos e a
  linha deixa de ser a mesma; a montagem recusa em vez de somar como se fossem um.
  É o próximo obstáculo, e já estava apontado antes deste trabalho;
- **Robin não está montado.** Dirichlet e Neumann sim. Robin mistura valor e
  derivada na mesma faceta: a linha ganha coeficiente que depende de α, β e h, e
  montá-lo como se fosse Neumann daria resultado plausível e errado;
- **serial.** `np=1`. Em paralelo entram a franja e a repartição, que são os itens
  seguintes da lista;
- **o piso do erro é o solver.** Pedir menos de 1e-13 seria exigir do solver linear
  algo que ele não promete com tolerância 1e-14.

## 7. Como reproduzir

```bash
cd t8flow && make 2d
mpirun -n 1 ./exemplos/poisson-mms --solucao=grau2 --contorno=dirichlet --niveis=3,4,5
mpirun -n 1 ./exemplos/poisson-mms --ajuda      # lista as solucoes do catalogo
```

Parâmetros: `--solucao` (grau1, grau2, grau3, trig), `--contorno` (dirichlet,
misto, umdirichlet), `--fechamento` (3 ou 2), `--niveis`, `--tol`. A configuração
inteira vai impressa no alto da saída — relatório sem a configuração é número sem
procedência, e neste projeto um build sem o que se queria medir já produziu duas
corridas "estáveis" que eram outra coisa.
