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
