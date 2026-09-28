# Relatório da corrida autônoma — front-tracking (B1 → B2)

Log de experimentos, um bloco por passo. Escrito durante uma ausência do usuário,
com direção: ir até o acoplamento com o solver (teste de Laplace), commitar cada
checkpoint verificado, e relatar cada passo com os experimentos. Sessão iniciada
em 2026-09-26.

Princípio herdado da sessão: **medir, não presumir**. Cada passo declara o
ORÁCULO (o que provaria que está certo) antes do número, e relata o que de fato
aconteceu — inclusive falhas.

Decisões de formulação adotadas por padrão nesta ausência (revisáveis):
- **Propriedades iguais** (ρ, μ iguais nas duas fases) no teste de Laplace, para
  ISOLAR a tensão superficial de qualquer efeito de salto de propriedade
  (decisão #2/#3 do projeto adiada de propósito).
- **Um fluido com força singular** σκn·δ espalhada na malha (não dois domínios
  acoplados).

---

## Passo 1 — checkpoint do B1

**Oráculo:** o arreio de Rider-Kothe já verificado na sessão anterior compila e
roda a partir da árvore do repositório com os mesmos números.

**Feito:** commit `6bf7746` — módulo `hig-flow-front-tracking.{c,h}` + arreio
`higflow/tests/front-tracking/{teste-b1.c,Makefile}`.

**Resultado (reprodução):**

    COM cirurgia    nmarc   dt       cons_max    reversao   n_final
                       64  0.020    1,71e-3     1,02e-3      166
                      128  0.010    4,05e-4     1,69e-4      334
                      256  0.005    9,90e-5     3,58e-5      654

Advecção RK2 de 2ª ordem; cirurgia guardada conserva área a 6e-5 no ida-e-volta.

## Passo 2 — portão do B1 (oráculo de CI)

**Oráculo:** o arreio deve devolver código de saída — 0 se, no mais fino (256,
dt=0,005) com cirurgia, `cons_max < 2e-4` E `reversao < 1e-4`; 1 caso contrário.
Limiares com margem folgada sobre o medido (9,9e-5 / 3,6e-5), pegando regressão
de sinal, escala ou ordem.

**Feito:** `teste-b1.c` agora imprime a linha PORTAO e retorna o código; removida
a escrita de VTK com caminho quebrado. Doc do projeto atualizado (B1 feita).

**Resultado:** `PORTAO B1 ... ==> PASSOU`, `EXIT=0`.

## Passo 3 — curvatura κ e normal (física do B2)

**Oráculo:** curvatura por círculo osculador de 3 marcadores consecutivos,
contra dois casos analíticos — círculo (κ=1/R em todo marcador) e elipse
(κ(t)=ab/(a²sin²t+b²cos²t)^{3/2}). Erro relativo deve cair com o refino.

**Feito:** `ft_vetor_curvatura` (κ·n, aponta para o centro osculador — o lado
côncavo — com módulo 1/R, independente da orientação da curva) e `ft_curvatura`
(escalar) no módulo. Arreio `teste-curvatura.c`; `make check` roda os dois
arreios como oráculos.

**Resultado:**

    circulo   nmarc  err_rel_max        elipse   nvert  err_rel_max
                32   5,3e-15                       128   1,8e-3
                64   4,4e-14                       256   4,5e-4
               128   1,7e-13                       512   1,1e-4
               256   7,3e-13                      1024   2,8e-5

Círculo EXATO (precisão de máquina — 3 pontos num círculo definem o próprio
círculo). Elipse converge em 2ª ordem. Portão PASSOU.

**Caveat descoberto e documentado** (custou uma iteração do arreio): o método de
3 pontos lê **κ=0 em tripla colinear**. A primeira versão do teste de elipse
usava `ft_cria_curva` com reamostragem, que põe marcadores colineares sobre as
cordas → κ=0 onde a elipse tem κ≠0 → erro de 100%. Não é defeito da curvatura (o
tip dava 19,86 vs 20,0 correto), é da inicialização. Consequência real para o
B2: **marcador recém-inserido pela cirurgia (ponto médio, sobre a corda) tem κ=0
até ser advectado** para fora da reta. No B2 estático a frente quase não se move
e a cirurgia raramente dispara, então a curvatura é medida na distribuição
inicial bem-posta — mas fica o registro para quando a interface se mover (B3+).

## Passo 4 — força de tensão superficial σκn (kernel de Roma)

**Oráculo:** dois testes que o arranjo standalone permite, ambos no nível de
máquina para uma gota em repouso:
- **conservação** — Σ_grade F·h² = σ·(integral do vetor curvatura), porque o
  kernel de Roma soma 1 (partição da unidade). Testa o espalhamento.
- **equilíbrio** — σ·∮κ⃗ds → 0 para curva fechada. Testa que a curvatura fecha
  em zero ao redor do círculo (força líquida nula = Laplace).

**Feito:** `ft_delta_roma` (mesmo kernel do corpo rígido, repetido para o módulo
ficar standalone), `ft_espalha_tensao` (espalha σκ⃗ numa malha uniforme centrada
na célula, peso ds_k por marcador), `ft_integral_curvatura`. Arreio
`teste-tensao.c` no `make check`.

**Resultado (círculo R=0,15, σ=0,7, malha 64²):**

    nmarc    conserv (x,y)      equilib |F|
       64    3e-16, 1e-15        2,9e-14
      128    1e-16, 3e-16        1,6e-13
      256    6e-16, 9e-16        2,8e-13
      512    5e-16, 8e-16        1,2e-12

Conservação no nível de máquina (~1e-16) — o espalhamento de Roma é exato.
Equilíbrio ~1e-13 — a integral do vetor curvatura no círculo é zero (por simetria
do círculo, é máquina-zero, não convergência; é a resposta forte para o alvo de
Laplace). Portão PASSOU.

**Nota de honestidade sobre o oráculo de equilíbrio:** ser máquina-zero vem da
simetria do círculo. Uma curva fechada assimétrica daria equilíbrio convergindo
de um valor finito — teste mais exigente, mas o círculo É o caso de Laplace do
B2, então máquina-zero é o resultado certo aqui. O teste FORTE deste passo é a
conservação (1e-16), que não depende de simetria.

## Passo 5 — acoplamento com o solver (em progresso)

**Arquitetura escolhida** (menor risco às suítes): o núcleo
`hig-flow-front-tracking.{c,h}` fica puro/standalone; o adaptador
`examples-common/front-tracking.c` faz a ponte com o solver e é compilado
**por-exemplo** (2D), fora da lib — evita o risco de DIM=3 e de regressão nas
suítes. Única mudança em código compartilhado: **expor `_suporte`** do corpo
rígido como `fi_suporte_facetas` (wrapper puro, comportamento idêntico), para o
adaptador reusar a busca de faceta escalonada sem duplicar a parte delicada.

**Ponto de acoplamento:** o gancho `fronteira_imersa_aplica` roda antes do
preditor (hig-flow-step.c:1587), depois de `facet_source_term` setar `dpFU`
(1571) — então somar σκ⃗ em `dpFU` é a força de corpo que a equação de momento
lê. Mesmo gancho do corpo rígido.

**Decisão de formulação (documentada para revisão):** um fluido só com força
singular na interface, propriedades iguais — isola a tensão superficial. SERIAL
por enquanto (sem reduce de franja; np=1 não tem franja).

**Feito e verificado (passos 5a-5c):**
- `fi_suporte_facetas`/`fi_suporte_capacidade` expostos; corpo rígido recompila.
- `ft_espalha_tensao_solver` no adaptador: espalha σκ⃗ em `dpFU` escalonado,
  peso ds_k, reusando `fi_suporte_facetas`. `ft_forcas_tensao` no núcleo entrega
  (posição, σκ⃗, ds) por marcador.
- Exemplo `example2d_FrontTracking` (gota R=0,25, σ=1 → Laplace Δp=4,0) constrói
  e roda.

**Oráculo (conservação do espalhamento no domínio REAL do solver, 2 passos):**

    FT conserva dim=0: grade=-2,80e-13  marcadores=-2,82e-13  erro=2,3e-15
    FT conserva dim=1: grade= 2,41e-12  marcadores= 2,41e-12  erro=3,7e-16

A força espalhada na malha escalonada IGUALA a dos marcadores a ~1e-15 — o cerne
do acoplamento está correto. (O integral dos marcadores ~1e-13 confirma o
equilíbrio de Laplace no círculo.)

**Passo 5d — o teste de Laplace, FEITO e PASSOU.**

**Oráculo:** gota em repouso → salto de pressão Δp = σ/R através da interface, e
correntes parasitas (max|u|) pequenas e estáveis. Setup estático: entrada zerada
(`boundary_velocity` id=0 → 0), gota R=0,25 em (2,0), σ=1 → Δp esperado = 4,0.
O `ns-example-2d.c` amostra a pressão no centro (dentro) e em (2, 0,9) (fora).

**Resultado (200 passos, dt=0,001, np=1):**

    LAPLACE  p_in=4.001181  p_out=0.000568  Dp=4.000613  sigma/R=4.000000  erro_rel=0.0002

- **Δp = 4,0006 vs σ/R = 4,0000 — erro relativo 0,02%.** Salto de Laplace exato.
- **Correntes parasitas: Vmax ≈ 8,2e-4, estáveis** (8,2039e-4 em passos
  consecutivos, não crescem). A gota fica em repouso, como Laplace exige.

O acoplamento completo produz a física correta — não é só andaime. A força
σκ⃗ espalhada na malha escalonada equilibra o gradiente de pressão e gera o salto.

## Passo 6 — advecção acoplada (os marcadores andam com o fluido)

É a **mudança 1** das três do projeto, e o que separa o front-tracking da
fronteira imersa rígida: onde o corpo rígido interpola u para calcular a força
que impõe u=0, aqui interpola u para **mover** o marcador. Mesma interpolação,
destino oposto.

**Implementação:** `_campo_da_malha` no adaptador interpola a velocidade da malha
escalonada nos marcadores — mesma fórmula de `fi_interpola`
(`u(X) = Σ u(faceta)·peso`, os pesos de Roma somando 1), reusando
`fi_suporte_facetas`. É passada como callback ao `ft_advecta` já verificado no B1.
Ligada por `FT_ADVECTA=1`; **o padrão continua frente fixa**, para não alterar o
resultado B2 já commitado.

**Ordem no passo:** o gancho roda pré-preditor, quando `ns->dpu` é u^n — a
velocidade final **já projetada** (discretamente livre de divergência) do passo
anterior. Então: move a frente com u^n → cirurgia → espalha a força nas posições
novas, que é o que o preditor vê. Campo congelado em u^n dentro do passo (não
existe u em t+dt/2), então o RK2 do `ft_advecta` vira avaliação de ponto médio no
**espaço** e o esquema é de primeira ordem no tempo — coerente com o resto do
acoplamento explícito.

**Oráculos (gota estática, 200 passos, `FT_ADVECTA=1`):**

    passo    n     area        dA/A      deriva    circ      unidade   max|u_marc|
        0  128  0.19627070  0.00e+00  0.00e+00  0.999799  1.23e-14   0.000e+00
       50  128  0.19626975  4.83e-06  2.49e-11  0.999799  1.23e-14   2.410e-04
      100  128  0.19626878  9.75e-06  7.99e-11  0.999799  1.24e-14   2.410e-04
      150  128  0.19626782  1.46e-05  1.50e-10  0.999799  1.24e-14   2.410e-04

- **Partição da unidade 1,2e-14** — os pesos somam 1 em precisão de máquina. É o
  oráculo forte: prova que um campo uniforme seria interpolado exatamente, e pega
  suporte incompleto (que daria velocidade pequena demais, em silêncio). Foi esta
  checagem que pegou a força "exatamente pela metade" no corpo rígido.
- **Área conservada a 1,5e-5 em 150 passos** — advecção por campo livre de
  divergência preserva a área fechada.
- **Deriva do centroide 1,5e-10** e **circularidade constante** (0,999799 é o
  valor de discretização do polígono de 128 lados, não deriva) — a gota fica
  parada e circular, como o equilíbrio de Laplace exige.
- **Laplace segue valendo:** Δp=4,000549 vs 4,0 → **0,01%**.

**A checagem que evitou um falso verde:** uma gota que fica parada é exatamente o
que se veria se a interpolação devolvesse **zero em silêncio** — área perfeita,
deriva nula, circularidade constante. Por isso o `max|u_marc|`: ele dá 2,41e-4
(não-nulo, e **exatamente 0,000e+00 no passo 0**, quando a condição inicial é
u=0, então o número acompanha a realidade do campo). O campo euleriano tem
Vmax=8,05e-4; os marcadores veem 2,41e-4 porque o kernel suaviza os picos do
campo parasita no suporte de 1,5 células.

**Sem regressão:** com o flag desligado o Laplace reproduz o commit `ffc8718`
bit a bit (Δp=4.000613, erro 0,0002) e nenhum diagnóstico vaza.

**O que este teste NÃO exercita, honestamente:** a gota não deforma, então `n`
fica 128 o tempo todo — a **cirurgia não dispara** no regime acoplado, e o caveat
da curvatura (marcador recém-inserido sobre a corda tem κ=0) segue não exercitado
aqui. Os dois entram em cena no B3 (oscilação) e no B4 (bolha subindo).

## Passo 7 — comparação com o VOF (gota estática)

**Por que é legítima:** o `example2d_VOF` já É o teste de Laplace, sem adaptação
— gota circular R=1/6 em (0,5;0,5), duas fases newtonianas, tensão superficial
ligada, u=0 inicial. Alinhei tudo o que importa: mesma malha ([0,1]², 61×61,
h=1/61), mesma gota, **10,2 células por raio nos dois**, ρ=μ=1, σ_efetivo=1,
Re=1, dt=1e-3, 101 passos, np=1, mesmo solver, mesmas sondas.

O σ efetivo do VOF não era óbvio: a força entra como `IF/(We·ρ)` com
`We = Ca·Re = 1`, logo σ=1 — o mesmo do FT. Verificado no código, não suposto.

**Um confundidor removido:** o caso VOF impõe na parede superior uma velocidade
que *cresce no tempo* (`8(1+tanh(8t-4))x²(1-x)²`). Com ela ligada a velocidade
medida não é corrente parasita — é parasita **mais escoamento forçado**. Zerada
atrás de `VOF_ESTATICO=1`, preservando o padrão do exemplo.

**Resultado:**

| métrica | front-tracking | VOF | razão | melhor |
|---|---|---|---|---|
| salto de pressão (erro rel.) | **7,67e-04** | 4,73e-03 | 6,2× | **FT** |
| correntes parasitas (max\|u\|) | 5,06e-04 | **5,23e-05** | 9,7× | **VOF** |
| deriva de massa no tempo | 1,29e-05 | **2,86e-15** | 4,5e9× | **VOF** |
| forma na inicialização | 1,61e-03 | **1,27e-05** | 127× | **VOF** |

Ambas as correntes parasitas são estáveis e decaem; nenhuma cresce.

**Dois resultados confirmam a literatura, um a contraria:**

1. **FT ganha no salto de pressão (6×)** — esperado: a curvatura do FT é exata
   no círculo (7e-13, medido no passo 3); a do VOF vem de derivadas de um campo
   de fração suavizado.
2. **VOF conserva massa em zero de máquina** — esperado: advecta um escalar
   conservado.
3. **VOF ganha nas correntes parasitas (10×)** — **contraria** a expectativa.

**A explicação da nº 3, e é testável:** o problema não é a curvatura (onde o FT
é melhor), é o **balanço discreto**. Para a gota em equilíbrio o contínuo exige
∇p = σκn·δ, e a velocidade parasita mede quanto os *operadores discretos* dos
dois lados falham em satisfazer isso. No VOF/CSF a força é σκ∇f, e com κ quase
constante no círculo isso é ≈ ∇(σκf) — quase exatamente um gradiente discreto,
cancelável pelo gradiente de pressão (propriedade *balanced-force*, Francois et
al. 2006). No meu FT a força é espalhada pelo kernel de Roma em facetas
escalonadas, e esse operador **não é** o gradiente discreto de escalar nenhum
que a projeção produza. Coerente com o salto sair melhor: a *integral* da força
está certa (conservação 1e-15), a *distribuição* não se escreve como gradiente.

*Como testar:* espalhar um potencial e tomar dele o **mesmo** gradiente discreto
que a projeção usa. Se as parasitas caírem sem piorar o salto, confirma.

**Um erro de medição meu, corrigido no caminho:** a primeira comparação dava
empate em massa (1,3e-5 dos dois lados). Era artefato de **medir coisas
diferentes** — no VOF eu comparava com a área analítica (inclui inicialização),
no FT com a área inicial (deriva pura). Separadas, a 1,3e-5 do VOF é erro de
*inicialização* e a deriva verdadeira é zero de máquina.

**Limitações:** uma única resolução (10,2 células/raio, sem estudo de
convergência em h), um único tempo (t=0,101), serial, propriedades iguais, e
gota estática — que não deforma, logo não exercita a cirurgia.

**Relatório LaTeX para Overleaf:** `doc/relatorio-front-tracking/` (main.tex +
secoes/ + dados/ + referencias.bib). As 11 referências foram conferidas no
Crossref, incluindo a linhagem do sistema: Freeflow (Castelo et al., 2000),
GENSMAC (Tomé & McKee, 1994), o front-tracking do próprio grupo (de Sousa et
al., JCP 198, 2004) e o HiGFlow (Sousa et al., JCP 396, 2019). **Não compilado**
— não há pdflatex na máquina; o que foi verificado está no LEIAME.

## Passo 8 — a hipótese do balanço, TESTADA E CONFIRMADA

**A hipótese:** as correntes parasitas não vêm da curvatura (onde o FT é
melhor), vêm do **balanço discreto**. A força precisa *pertencer à imagem do
gradiente discreto que a projeção inverte* — não basta cada lado estar certo.

**A reformulação** (`FT_BALANCEADO=1`): força remontada como **F = σκ∇H**, com
- **H** = fração de área da gota por célula, da *geometria* da frente (recorte
  de polígono Sutherland–Hodgman + fórmula do laço) — `ft_area_na_caixa`;
- **∇** = o **mesmo** operador que `higflow_final_velocity` usa para a pressão
  (`compute_center_p_left/right` + `compute_dpdx_at_point`). Essa é a exigência
  que faz "balanceado"; outro gradiente, ainda que consistente, não cancelaria.

Com κ constante (o círculo), isso é exatamente ∇(σκH) → a pressão cancela termo
a termo.

**Indicadora verificada antes de usar:** somar H·área_célula sobre uma grade que
cobre a frente devolve a área do polígono — **3,7e-15** (4 frentes × 4 grades).
Célula cheia 4,0e-14, célula vazia exatamente 0.

**Resultado — mesmo caso, só trocando a formulação:**

| métrica | FT espalhado | FT balanceado | fator |
|---|---|---|---|
| erro rel. de Δp | 7,67e-04 | **1,7e-07** | 4500× |
| correntes parasitas (final) | 5,06e-04 | **<5e-11** | >10⁴× |
| deriva de área | 1,29e-05 | **1,9e-15** | 7e9× |
| deriva do centroide | 3,4e-07 | **1,9e-14** | 2e7× |

**Confirmada nos dois sentidos que previa:** as parasitas caem >4 ordens **e** o
salto de pressão *melhora* (não piora). Se a curvatura fosse o problema, as duas
coisas não andariam juntas.

**A assinatura temporal é mais clara que os números finais:** o espalhado sobe e
**estaciona** em 5e-4 (forçamento espúrio persistente que a viscosidade
equilibra mas não elimina); o balanceado **decai** monotonicamente — 2,4e-7,
1,6e-9, 2e-10 — até sumir sob a precisão de impressão. Não há o que equilibrar.

**Falso verde descartado:** velocidade nula é também o que eu veria se a força
não estivesse sendo aplicada. O que descarta é o próprio salto: **Δp = 6,000001**
contra o exato 6. Sem força não haveria salto nenhum. As duas medidas juntas —
pressão certa *e* velocidade nula — são a assinatura do equilíbrio; nenhuma
sozinha bastaria.

**Consequência para a comparação** — o FT balanceado passa a superar o VOF nas
três métricas:

| métrica | FT balanceado | VOF | vantagem |
|---|---|---|---|
| erro rel. de Δp | **1,7e-07** | 4,73e-03 | 2,8e4× |
| correntes parasitas | **<5e-11** | 5,23e-05 | >10⁶× |
| deriva de massa | 1,9e-15 | 2,9e-15 | equivalentes |

A massa deixou de ser desvantagem, mas **por uma razão que não se deve
superestimar**: a área não deriva porque o escoamento *é* estático, e marcador
parado não perde área. Em escoamento de verdade a deriva volta, e a conservação
por construção do VOF continua sendo vantagem *estrutural*.

**A limitação mais importante:** o balanço é *exato* só com κ constante, pois aí
σκ∇H = ∇(σκH). Com curvatura variável sobra um resíduo ~H∇(σκ) que **este teste
não mede** — e é justamente ele que governaria uma gota deformada. O próximo
teste tem de ser gota elíptica.

**Sem regressão:** o caminho padrão (espalhado) reproduz `Dp=5.995400` bit a bit.

## Passo 9 — gota elíptica: o teste que achou um defeito

**O caso:** elipse a=0,2357 b=0,1179 (razão 2:1, **mesma área** do círculo R=1/6),
que relaxa para o círculo de mesma área → R_eq=1/6 → Δp final esperado = **6,0**.
1200 passos (t=1,2, ~4 tempos de relaxação). É também o **primeiro caso em que a
frente deforma**, logo o primeiro que exercita a **cirurgia** com o solver no
circuito (n: 108 → 102-104).

**O balanceado FALHOU.** Circularidade sobe normal até o passo 900 (0,9887) e
então **reverte** — 0,976 → 0,962 — com a velocidade crescendo. Δp = 10,99
contra 6,0: **83% de erro**. A falha aparece exatamente quando a gota fica quase
circular, isto é, quando a força física que a movia já decaiu e o artefato ficou
sozinho.

**A causa era a simplificação que eu já tinha registrado como limitação:** κ na
faceta vinha do **marcador mais próximo**. No círculo é exato (κ constante). Na
elipse é um campo **descontínuo** — salta quando a atribuição troca de marcador.
As descontinuidades movem a frente, o que muda a atribuição, o que muda a
força: realimentação.

**Conserto dirigido confirmou o diagnóstico:** troquei κ por **média ponderada**
dos marcadores vizinhos (kernel de Roma, largura h). Nada mais. A instabilidade
**sumiu**, e a trajetória passou a coincidir com a do espalhado até a 4ª casa
(0,998477 vs 0,998491). **Sem regressão no círculo**: Δp=6,000001, u<5e-11 —
κ constante, e a média ponderada de constante devolve a constante.

**Resultados (relaxação da elipse):**

| grandeza | espalhado | bal. κ vizinho | bal. κ suave | VOF |
|---|---|---|---|---|
| circularidade final | 0,998491 | **0,961697** (revertendo) | 0,998477 | — |
| Δp final (exato 6,0) | 5,9239 | **10,99** | 5,9089 | 5,9285 |
| erro rel. do salto | 1,27% | **83%** | 1,52% | 1,19% |
| deriva de massa | 1,09e-04 | 3,20e-05 | 2,25e-04 | **2,47e-13** |

**Três leituras, duas delas corrigindo impressões anteriores:**

1. **Na elipse o balanço NÃO traz vantagem** — 1,27% / 1,52% / 1,19%, os três
   dentro de 0,3% um do outro. O ganho espetacular do passo 8 era **específico
   do caso estático com κ constante**, onde o balanço exato é alcançável.
2. **A vantagem estrutural do VOF em massa aparece agora.** Eu havia advertido
   que a conservação perfeita do balanceado vinha de o escoamento *ser* estático
   e que "em escoamento de verdade a deriva volta". Voltou: 1e-4/2e-4 contra
   **2,5e-13** do VOF. Este é o primeiro teste com movimento real.
3. **Os três concordam na física** — 5,909 / 5,924 / 5,929 para o exato 6,0. O
   déficit comum de 1,2% não é de método: é relaxação incompleta (circ ainda
   0,9985, velocidade residual ~9e-3 contribuindo pressão dinâmica). Três
   discretizações independentes errando *junto* e pelo mesmo tanto é o que se
   espera de erro físico compartilhado.

**A lição de método:** o resultado do círculo era verdadeiro e continua válido —
faltava saber *sobre o que* ele falava. Um caso em que uma quantidade é
constante não distingue uma implementação que a trata bem de outra que a trata
de qualquer jeito, e aqui a diferença era entre funcionar e ir à instabilidade.

## Passo 10 — B3: oscilação de gota contra a frequência de Lamb

**Oráculo:** ω² = n(n²−1)σ/[R³(ρ_i+ρ_o)] → para n=2, σ=1, ρ=1, R=1/6:
**T = 0,24683**. É oráculo de *frequência*, insensível à calibração da medida.

**Dois parâmetros mudaram, por razão medida:**
- Re=1 é **fortemente superamortecido** (β≈144 contra ω=25) — por isso a elipse
  relaxou monotonicamente. Usei **Re=200**.
- A restrição capilar dá dt < 8,4e-4: **o dt=1e-3 que eu vinha usando estava
  acima do limite**. Passei a 5e-4.

**A medida:** D = (a−b)/(a+b) pelos **momentos de área** — no FT em forma fechada
do polígono, no VOF da integral da fração. *Mesma* grandeza dos dois lados.
Verificada contra 4 elipses: erro 1,3e-5.

| grandeza | front-tracking | VOF | teoria |
|---|---|---|---|
| período medido | 0,26791 | 0,26210 | **0,24683** |
| erro vs Lamb | 8,54% | 6,19% | — |
| amortecimento | 1,613 | 1,377 | ~0,90 |
| oscilações limpas | **7** | **7** | — |
| deriva de massa | **3,9%** | **0** | 0 |

Os dois dão **7 oscilações amortecidas limpas**, período estável, e concordam
entre si a **2,2%**. O excesso comum de ~7% sobre Lamb é confinamento (paredes a
2R) — a correção viscosa da frequência responde por <0,5%. Dois métodos
independentes errando para o *mesmo lado* e quase pelo mesmo tanto é assinatura
de efeito físico compartilhado.

**A conservação de massa sob movimento sustentado é a diferença decisiva.** A
progressão conta a história inteira: gota estática, o FT parecia conservar
perfeitamente (1,9e-15) — mas só porque nada se movia; elipse relaxando, a
deriva apareceu (1,1e-4); **oscilação sustentada, 3,9%**. O VOF manteve
`A=0.08726527` em *todos* os dígitos nos três regimes. Isso vem da advecção dos
marcadores, e nenhum ajuste da força remove.

**O balanceado falhou catastroficamente** e a corrida foi abortada: marcadores
64 → **136.707**, área −27,5%, circularidade 0, partição da unidade degradada de
6e-15 para **0,833**. Com o amortecimento menor, a instabilidade que a suavização
de κ havia *contido* (não curado) reapareceu e despedaçou a frente.

**Conclusão sobre o balanço, revista pela terceira vez:** espetacular na gota
estática (4 ordens), sem vantagem na elipse, **inutilizável** sob oscilação pouco
amortecida. **A formulação padrão continua sendo a espalhada.**

**Um erro de montagem que o próprio oráculo pegou:** a primeira corrida do VOF
não oscilou. O σ efetivo dele é 1/(Ca·Re) — então **subir Re de 1 para 200
dividiu o σ do VOF por 200**, deixando-o superamortecido, enquanto o FT usa σ
explícito. Corrigi com Ca=1/Re=0,005. A relação tinha sido *verificada no código*
quando Re valia 1, e ainda assim escapou ao mudar o parâmetro: **uma equivalência
conferida uma vez não permanece conferida quando se mexe em algo de que ela
depende.**

**Figuras:** o relatório LaTeX ganhou 3 figuras geométricas (formas reais das
frentes, de `dados/formas/*.dat`): Rider-Kothe círculo→filamento→círculo,
relaxação da elipse, e as 4 fases de um período da oscilação.

## Resumo da corrida autônoma

Cinco passos, todos verificados por oráculo e commitados:
1. B1 (advecção + cirurgia) — Rider-Kothe reversível, área a 6e-5.
2. Portão do B1 (código de saída para CI).
3. Curvatura — círculo exato, elipse 2ª ordem.
4. Força de tensão σκn — conservação 1e-16 (standalone).
5. Acoplamento com o solver — conservação no domínio real 1e-15; **Laplace
   Δp=σ/R a 0,02%, correntes parasitas 8e-4**.

**B2 (gota estática, lei de Laplace) alcançado.** Aberto para quando o usuário
voltar: (a) advecção acoplada + correntes parasitas ao longo do tempo (o campo
interpolado nos marcadores, `FT_ADVECTA`); (b) paralelo (reduce de franja no
espalhamento); (c) B3 (oscilação de gota, frequência de Lamb). Decisões de
formulação (#2/#3, salto de propriedade) seguem adiadas — propriedades iguais
bastaram para o Laplace.
