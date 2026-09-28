# Relatório: front-tracking bifásico no HiGFlow

Projeto LaTeX. No Overleaf: envie a pasta inteira; o arquivo principal é
`main.tex`.

```
main.tex            preâmbulo e ordem das seções
secoes/*.tex        uma seção por arquivo
referencias.bib     11 entradas, todas conferidas no Crossref
dados/*.dat         AS MEDIDAS
figuras/            (vazio por enquanto)
```

Compila com `pdflatex + bibtex + pdflatex + pdflatex`, ou `latexmk -pdf main`.

## Os números vivem nos `.dat`, não no texto

Mesma convenção do `relatorio-fronteira-imersa`: as tabelas e o gráfico são
gerados por `pgfplotstable`/`pgfplots` lendo `dados/*.dat`. Regenerar uma medida
é trocar o arquivo; o relatório acompanha. Evita o erro clássico de a tabela e o
texto discordarem depois de uma correção.

| arquivo | o que traz |
|---|---|
| `b1-rider-kothe.dat` | B1: conservação e reversão, com e sem cirurgia |
| `curvatura.dat` | erro de curvatura contra círculo e elipse analíticos |
| `tensao-conservacao.dat` | espalhamento: conservação e equilíbrio |
| `adveccao-acoplada.dat` | advecção acoplada: área, deriva, circularidade, partição da unidade, `max\|u_marc\|` |
| `comparacao-laplace.dat` | **a comparação FT × VOF**, métricas relativas |
| `comparacao-bruto.dat` | a mesma comparação em valores brutos |
| `parasitas-serie.dat` | série temporal das correntes parasitas nos dois métodos |
| `hysing-referencia.dat` | os valores publicados por Hysing et al., três códigos e o consenso |
| `hysing-parcial.dat` | Hysing: $y_c$ e $V_c$ nos dois métodos, **trecho parcial** $t\in[0;0{,}35]$ |
| `hysing-forma.dat` | Hysing: concordância de forma entre os dois métodos |
| `paralelo.dat` | o que $n_p=2$ mede nos dois métodos (seção de discussão) |
| `formas/hy_ft_*.dat`, `formas/hy_vof_*.dat` | a interface nos dois métodos nos mesmos instantes |

## Estado da verificação deste projeto LaTeX

**Não foi compilado.** Não há `pdflatex` nem `bibtex` na máquina onde o relatório
foi escrito, então o PDF não foi gerado nem conferido. O que *foi* verificado,
por conferência estática automatizada:

- todos os `\input` existem;
- todos os `\dat{...}` apontam para arquivos que existem;
- toda coluna referenciada em `columns/<nome>/` existe no cabeçalho do `.dat`
  correspondente, **e** toda coluna do `.dat` tem `column name` (nenhuma sairia
  com o nome cru);
- toda `\cite` tem entrada no `.bib`; toda `\ref`/`\eqref` tem `\label`;
- chaves balanceadas em cada arquivo;
- número de colunas consistente em cada `.dat`.

Sublinhados foram removidos dos `.dat` de propósito: em coluna `string type` eles
quebram a compilação (`Missing $ inserted`).

Ao abrir no Overleaf, espere ajustes de espaçamento de tabela — a largura das
colunas não pôde ser conferida aqui.

## As referências

Todas as 11 entradas do `referencias.bib` foram conferidas contra o registro
canônico do editor via API do Crossref (`api.crossref.org/works/<DOI>`), não
contra memória nem citação de terceiros. Título, autores na ordem, volume,
páginas e ano vieram de lá.

Duas observações que ficaram registradas no próprio `.bib`:

- O Crossref transcreve o primeiro autor do artigo do Freeflow como
  "A. Castello F."; o correto é **Castelo**, e é assim que consta.
- O tutorial do sistema FreeFlow data o artigo de 1999; o registro do editor diz
  2000 (Comput. Vis. Sci. 2(4), março de 2000). Prevalece o publicado.

## A seção do Hysing está incompleta de propósito

`secoes/10-hysing.tex` monta o *benchmark*, traz os valores publicados, registra
o erro de definição da circularidade e declara de antemão como ler cada
desfecho — mas **não** fecha contra os valores publicados. As corridas até
$t=3$ (12 000 passos, ~6 h cada) não terminaram a tempo do texto. O que está
medido é `hysing-parcial.dat`, o trecho $t\in[0;0{,}35]$, e a subseção que o
apresenta diz isso na primeira linha.

Para fechar: substituir a subseção "Resultado parcial" pelos três números de
cada método ($c_{\min}$ e seu instante, $V_{c,\max}$ e seu instante,
$y_c(t=3)$) contra `hysing-referencia.dat`.

As figuras de forma são regeneráveis:

```bash
higflow/scripts/gera-formas-hysing.sh 0.0 1.0 2.0 3.0
```

Ele copia o despejo lagrangeano do front-tracking e extrai o contorno
`FracVol=0,5` do VTK do VOF (via `higflow/scripts/contorno-vof.py`) para os
mesmos instantes, escrevendo `dados/formas/hy_ft_<k>.dat` e `hy_vof_<k>.dat`. A
figura do relatório lê esses nomes, então trocar os instantes não exige editar o
`.tex` — só ajustar os rótulos do eixo e a janela do painel ampliado.

O extrator do contorno foi conferido contra o círculo analítico em $t=0$: desvio
radial máximo 1,34e-4, erro de área 3,25e-4. Ele mediania os valores de célula
nos vértices antes do *marching squares*, o que **suaviza** — serve para
desenhar a forma, não para medir perímetro, e é por isso que a circularidade do
lado VOF continua não sendo reportada.

A linhagem do código está citada explicitamente na introdução: Freeflow →
GENSMAC → o front-tracking do próprio grupo (de Sousa et al., JCP 198, 2004) →
HiGFlow (Sousa et al., JCP 396, 2019).
