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

**Falta (passo 5d):** o teste de Laplace de fato — BCs estáticas (sem entrada,
a copiar do canal seria fluxo forçado), medir Δp interno/externo vs σ/R e as
correntes parasitas max|u|. É a parte que depende de setup de contorno cuja
correção só o próprio oráculo (u≈0 e Δp certo) confirma.
