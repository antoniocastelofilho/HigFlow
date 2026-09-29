# Relatório: front-tracking bifásico no HiGFlow

Projeto LaTeX, autossuficiente. **No Overleaf:** envie a pasta inteira (ou o
`relatorio-front-tracking.tar.gz`); o arquivo principal é `main.tex`, o compilador
é pdfLaTeX e a bibliografia sai por BibTeX — todos os padrões do Overleaf.

Nenhum caminho é absoluto e nada é referenciado fora da pasta. São 13 `.tex`,
55 `.dat`, 1 `.bib`, e o PDF tem 37 páginas.

Compila com `pdflatex + bibtex + pdflatex + pdflatex`, ou `latexmk -pdf main`.
Verificado a partir do zero (`latexmk -C` e recompilar): zero erros, zero
estouros de margem, zero referências ou citações não resolvidas.

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
| `paralelo-consertos.dat` | o que foi consertado no paralelo, e o que se revelou não ser defeito |
| `paralelo-posse.dat` | as três regras de posse no espalhamento, medidas |
| `hysing-vcmax.dat` | o primeiro alvo publicado fechado: $V_{c,\max}$ nos dois métodos |
| `hysing-divergencia.dat` | a divergência entre os métodos crescendo com a deformação |
| `hysing-massa-amr.dat` | a conservação do VOF em malha fixa contra malha adaptativa |
| `hysing-fechamento.dat` | **os cinco alvos fechados**, com o desvio em unidades da dispersão dos códigos |
| `hysing-final.dat` | o que as corridas custaram: passos, retomadas, deriva de massa, marcadores |
| `formas/hy_ft_*.dat`, `formas/hy_vof_*.dat` | a interface nos dois métodos nos mesmos instantes |

## Estado da verificação deste projeto LaTeX

**Compilado.** `latexmk -pdf main` roda limpo: 32 páginas A4, **zero erros, zero
estouros de margem, zero referências ou citações não resolvidas**. O PDF está em
`main.pdf`.

Além da compilação, uma conferência estática automatizada verifica:

- todos os `\input` existem, e nenhuma seção está órfã;
- todos os `\dat{...}` apontam para arquivos que existem;
- toda coluna referenciada em `columns/<nome>/` existe no cabeçalho do `.dat`,
  **e** toda coluna do `.dat` tem `column name`;
- toda `\cite` tem entrada no `.bib`; toda `\ref`/`\eqref` tem `\label`;
- chaves balanceadas em cada arquivo;
- número de colunas consistente em cada `.dat`, **contando com reconhecimento de
  `{agrupamento}`** — sem isso, uma célula com espaços dentro de chaves conta
  errado.

### Duas armadilhas que só a compilação pegou

**Tabela vazia sem erro.** O `pgfplotstable` pode descartar todas as linhas de
dados e ainda assim compilar limpo — cabeçalho impresso, corpo em branco, código
de saída zero. Aconteceu com `hysing-fechamento.dat`, que é a tabela principal do
relatório. A causa era célula composta sem `{chaves}`. **Verificar código de saída
não basta:** é preciso conferir que o conteúdo apareceu.

**Janela de eixo cortando a curva.** O painel ampliado da figura da bolha tinha
`ymax` abaixo do topo das formas, e mostrava fragmentos. Nenhuma verificação
estática vê isso — só renderizar a página e olhar.

Para repetir as duas checagens: compile, extraia o texto (`pdftotext main.pdf -`)
e procure um número distintivo da última linha de cada tabela; e renderize as
páginas de figura (`pdftoppm -f N -l N -r 110 -png main.pdf saida`).

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

## A seção do Hysing está FECHADA

Os três alvos publicados foram medidos nos dois métodos (cinco resultados, porque
a circularidade só existe no front-tracking). Ver `sec:fechamento` para a tabela e
`sec:vies` para o achado principal: os cinco desvios são negativos, o que aponta
viés compartilhado pelos dois métodos e põe a convergência em $h$ como a pergunta
aberta.

A subseção `sec:leituras` ("O que decidiria") foi **mantida como estava**, com as
três leituras declaradas antes de os números existirem. Ela não foi reescrita de
propósito: alterá-la apagaria a evidência de que a interpretação foi fixada de
antemão. O mesmo vale para a subseção do resultado parcial, cuja conclusão a
`sec:divergencia` mostra não se estender ao regime deformado.

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
