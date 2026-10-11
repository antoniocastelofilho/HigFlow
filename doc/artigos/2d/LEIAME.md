# Os dois artigos: bidimensional e tridimensional

Mesma tese nos dois — os dois métodos de interface construídos sobre a **mesma**
maquinaria, de modo que a diferença medida seja atribuível à representação da
interface e não ao arranjo — e cada um no que a sua dimensão permite medir:

- **2D** (`2d/`): há um benchmark publicado com valores tabelados (Hysing
  et al.), então a comparação é contra número externo.
- **3D** (`3d/`): não há equivalente publicado, e a comparação é entre os
  métodos e contra oráculos de forma fechada.

Cada pasta é **autocontida**: `main.tex`, `preambulo.tex`, `secoes/`, `dados/`,
`referencias.bib`. Nenhuma imagem binária — todos os gráficos saem do `pgfplots`
lendo `dados/*.dat`.

## No Overleaf

Suba a pasta (ou o `.zip` correspondente) e marque `main.tex` como documento
principal. Compilador: **pdfLaTeX**. Classe: `elsarticle`, que o Overleaf já tem.

## Localmente

    latexmk -pdf main

Exige `texlive-publishers` (é de lá que vem o `elsarticle`).

## As três macros de honestidade

Herdadas da convenção desta árvore, e a distinção entre as duas primeiras
importa:

- `\missing{...}` — o que falta **medir**. Vermelho.
- `\caveat{...}` — o que **foi** medido, mas com uma ressalva que muda a
  leitura. Laranja. Confundir com `\missing` esconderia um resultado-com-reserva
  no meio das pendências.
- `\status{...}` — abre cada seção dizendo que respaldo ela tem. Definir
  `\nostatus` antes do `\begin{document}` faz todos sumirem, para a submissão.

E a convenção que sustenta as três: **os números vivem em `dados/*.dat`** e são
lidos pelo `pgfplotstable`, nunca digitados no texto. Regerar uma medida troca o
`.dat` e o texto acompanha, sem chance de a tabela e a prosa discordarem.

## Armadilhas já pagas, registradas para não se pagar de novo

**O BibTeX não trata `%` como comentário no nível de topo.** Qualquer arroba
fora de uma entrada é lida como início de entrada, mesmo numa linha que começa
com `%` — e o sintoma é bibliografia vazia com o PDF compilando sem reclamar.
Para desativar uma entrada, apague-a.

**`\nocite{*}` está nos dois `main.tex` e tem de sair na submissão.** Ele faz a
bibliografia listar todas as entradas do `.bib`, inclusive as ainda não citadas.
Sem ele, uma seção sem nenhum `\cite` sai com bibliografia vazia — foi o que
aconteceu com o 2D quando o removi: `I found no \citation commands`, e o PDF
compilou sem erro.

**Campo matemático dentro de um `.dat` precisa de `$...$`.** O `pgfplotstable`
compõe o conteúdo literalmente, e um `_` solto dá "Missing $ inserted".

**A `elsarticle` tem mancha mais estreita que a `article` com margem 2,7 cm.**
Três tabelas que cabiam passaram a transbordar na conversão; estão envolvidas em
`\adjustbox{max width=\linewidth}`. Uma tabela nova provavelmente precisará do
mesmo.

**O `microtype` exige fonte escalável.** Sem o `lmodern` no preâmbulo a
compilação morre com "auto expansion is only possible with scalable fonts".

## De onde saem as séries temporais

`higflow/examples-common/metricas-superficie.py` lê os VTK de superfície dos
**dois** métodos — o do marching cubes no VOF e o da triangulação no FT — de
modo que a comparação não herde diferença de instrumento. Em 2D o equivalente é
o `hysing-metricas.c`, ligado e chamado nos dois exemplos, pela mesma fórmula.

## Na submissão

- definir `\nostatus`
- remover o `\nocite{*}`
- preencher autores, data e resumo

## Estado, em uma linha cada

**3D — pronto para escrita final.** A matriz viscoelástica uniforme fechou (12
corridas), o par da cúspide na resolução de referência fechou, os doze oráculos
passaram, o custo em malha uniforme está medido. Falta o resumo, a introdução,
o caso 2 e o par adaptativo em produção.

Em 11/10/2026 entrou a explicação da deriva de massa do front-tracking, que
estava listada como pendência: ela só *parecia* crescer com o refino porque as
duas corridas de produção mantiveram `FT3_NSUB=4` fixo enquanto a malha
refinava, de modo que `ds/h` dobrou entre elas (as frentes finais têm 2823 e
2829 pontos — a superfície não refinou). Com `ds/h` constante a deriva cai em
ordem observada 1,87, e a cirurgia a **reduz** em 42% em vez de causá-la. Os
dados vieram de `relatorio-front-tracking/dados/tridim-{duasresolucoes,cirurgia}.dat`.

**2D — estrutura pronta, metade das medidas por fazer.**

| seção | estado |
|---|---|
| Formulação | portar do 3D, restrita ao plano |
| Método | escrita quanto ao implementado |
| Verificação | oráculos 2D por inventariar, com guarda de vacuidade |
| Newtoniano | **fechado** — os dois a t=3, 12.000 passos, cinco alvos contra o publicado, e o achado do viés compartilhado |
| Viscoelástico | **não medido**: VOF tem 2 corridas de 101 passos; FT tem deck e nenhuma corrida |
| Custo | não medido em 2D |

Os dados 2D **já estavam medidos** — moravam em `doc/relatorio-front-tracking/`,
não no artigo. Foram trazidos, com a vírgula decimal convertida em ponto e os
rótulos de linha traduzidos. O que falta ao 2D agora é só a metade
viscoelástica, que precisa de um par de corridas que não existe.
