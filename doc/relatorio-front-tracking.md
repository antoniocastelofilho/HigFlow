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
