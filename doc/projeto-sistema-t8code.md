# Um sistema sobre o t8code: escalonamento dele, multiblocos do HiGFlow

Levantamento de **11/10/2026**, a pedido do usuário. A pergunta era se é
possível construir um sistema novo sobre o t8code que tenha as características
de escalonamento dele e seja multiblocos como o HiGFlow.

**É possível, e boa parte do caminho já está andada neste repositório.** Este
documento separa o que existe do que falta, porque a resposta muda conforme
isso. Tudo abaixo foi apurado no código; onde é opinião, está dito.

---

## 1. O que já existe

### 1.1 Multiblocos — já casa

`higflow/examples-common/malha-t8.cxx` produz malha pelo t8code a partir da
mesma informação AMR que o solver leria de qualquer jeito. Levantamento que o
próprio cabeçalho registra: **oito exemplos têm um bloco uniforme, o
`example2d_Newt_contraction` tem dois** — e a fonte já o cobre, uma árvore por
bloco.

O `cmesh` do t8code é nativamente multiárvore. Essa é a peça que casa melhor:
o multiblocos do HiGFlow e a multiárvore do t8code são o mesmo conceito.

### 1.2 Escalonamento — já prototipado e verificado

`higtree/tests/test-particao-t8code.c` tem **a partição final sendo a do
t8code**:

> Até aqui o t8code produzia a malha e o `lbal` a reparticionava. Aqui não há
> `lb_calc_partition`: cada rank fica com o que a curva de preenchimento lhe
> deu, e o `partition_graph` é montado a partir das caixas.

E o oráculo é o certo. Não conta células nem volume — confere o **conteúdo da
franja** depois do `dp_sync`, porque a construção do grafo *falha em silêncio*
se o pareamento entre quem envia e quem recebe sair trocado: o domínio se monta,
o sync roda, e a célula de franja recebe o valor de outra posição.

### 1.3 Refino local — existe, com guarda geométrica

`t8_produz_por_rank_brick_refinado`, com `test-refino-caixa-t8code` afirmando:

- **razão no máximo 2:1** entre células vizinhas;
- **nenhuma árvore de uma célula** — uma árvore de uma célula não tem interior,
  toda faceta dela é de fronteira, e a franja declarada do HiGFlow é de cinco
  células. O cabeçalho do teste registra que esse defeito aparecia por dois
  caminhos distintos e *"só apareceu como escoamento divergindo numa simulação
  de meia hora"*.

### 1.4 A suíte

Cinco arreios em `higtree/tests/`: `particao`, `partition`, `production`,
`rank-production`, `refino-caixa`. O contrato de mesh de 18 cláusulas está
verificado nos **dois** backends, e o `ci/run_higtree_tests.py` lê as cláusulas
do próprio cabeçalho e exige caso passando para cada uma — cláusula sem teste
sai como falha, não como omissão.

Instruções de uso em `higtree/tests/t8code/README.md`.

---

## 2. O limite de hoje

O `malha-t8.cxx` declara em letras:

> **LIMITE, e ele decide quais exemplos podem usar isto: malha UNIFORME, um
> nível e um patch por bloco.** [...] Com refino a fonte RECUSA em vez de
> aproximar — uma malha diferente daria números diferentes, e o silêncio aí
> seria pior.

Ou seja: **hoje o t8code entra como produtor de malha, e a árvore hierárquica do
HiGFlow continua dona da estrutura.** Um sistema que tenha o escalonamento do
t8code precisa do inverso — a floresta sendo a estrutura, com o HiGFlow lendo
dela.

Essa inversão é o trabalho principal, e não é mecânica: muda quem é dono da
partição, de onde vem a franja, e quem decide o refino.

---

## 3. A armadilha de projeto a não herdar

Esta é a lição mais cara de 11/10/2026, e vale diretamente para o desenho novo.

O domínio da fração volumétrica do VOF (`sdmult`) é, **por construção, uma malha
uniforme no nível mais fino**, separada da árvore adaptativa (`sdp`). Está no
comentário de `higflow_create_amr_info_mult`:

> Use the given mesh information to create amr info about the multiphase domain.
> Corresponds to a uniform higtree with the cell sizes of the last given level.

As sete funções de advecção 3D percorriam `sdp` e usavam o índice do instantâneo
**direto** como índice de uma propriedade que vive em `sdmult`. Medido nos campos
salvos:

| propriedade | células |
|---|---|
| `dpfracvol` (em `sdmult`) | 375 000 = 50×50×150 |
| `dpu2` (na árvore `sdp`) | ~81 000 |

**Em malha uniforme os dois domínios coincidem célula a célula, e o índice
errado é acidentalmente o certo.** Por isso todo o 3D testado até então passou.
Com um nível de refino, 80 mil iterações escreviam num vetor de 375 mil em
posições sem relação geométrica com a célula de origem — e a fração parou de
advectar, em silêncio, enquanto a quantidade de movimento continuava correta.

O caminho 2D fazia certo, e o contraste é total: as funções de mesmo nome em
`hig-flow-step-multiphase.c` abrem com `psd_get_local_domain(ns->ed.mult.psdmult)`.
E não era só o domínio — o estêncil divergia por completo (2D lê `dpfracvol` com
`ns->ed.mult.stn` 49 vezes e com `ns->ed.stn` 0; o 3D fazia o inverso, 48 contra 0).

**Decisão de projeto para o sistema novo**, e ela é do usuário: um domínio só, ou
dois com tradução de índice **explícita e verificada**. Repetir "dois domínios
paralelos, um uniforme e um adaptado" sem tradução verificada herda a classe
inteira de defeito. Paga-se barato no início e caro depois.

---

## 4. O que recomendo verificar antes de comprometer o desenho

Em ordem de retorno:

1. **Medir escalabilidade do t8code no regime 3D.** O 2D rendeu só 1,6× em np=4
   porque 9377 células divididas por quatro dão 2300 por processo e a
   comunicação come o ganho. O 3D tem 375 mil células — quarenta vezes mais
   trabalho por processo, que é o regime onde a divisão vence. Sem essa medida,
   o escalonamento é promessa.

2. **Decidir a dona da estrutura** antes de escrever código: floresta do t8code
   com adaptador de leitura, ou árvore do HiGFlow com o t8code produzindo. As
   duas funcionam; só a primeira dá o escalonamento.

3. **Fechar os dois bloqueios do paralelo no FT 3D** — o sensor de esteira foi
   consertado em 11/10/2026 (`MPI_Allreduce`/`MPI_MAX`), mas a regeneração da
   fração a partir da frente ainda não reduz entre ranks, e o diagnóstico do
   arranque já acusa 17% de discrepância em np=2.

4. **Levar a disciplina de oráculo junto.** O que separou causa de sintoma em
   todos os casos deste dia foi experimento controlado com **uma** variável —
   mesmo diretório, mesma malha, mesmo binário, mudando só a que se quer
   isolar. E oráculo com **contraste**: "conserva área" não diz nada sem a média
   laplaciana ao lado, encolhendo 65 vezes mais.

---

## 5. Onde está o resto

- Estado da integração, contrato e suíte: `higtree/tests/t8code/README.md`
- Tutorial: `doc/tutorial-mtree-t8code.pdf`
- O defeito do índice e seu conserto: commit `4f8228e`
- A redução do sensor de esteira: commit `fa1eef8`
