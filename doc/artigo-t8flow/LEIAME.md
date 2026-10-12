# Artigo: t8Flow — estratégias e verificação do Poisson por MMS

Compila com

    pdflatex main && bibtex main && pdflatex main && pdflatex main

ou `latexmk -pdf main`. No Overleaf: enviar a pasta inteira; o principal é
`main.tex`.

## Os números vivem em dados/

Mesma convenção dos outros artigos desta árvore: **nenhum número é digitado no
texto**. As tabelas são lidas de `dados/*.dat` por `pgfplotstable`, e os `.dat`
são gerados rodando o código:

    cd doc/artigo-t8flow && ./gera-dados.sh

O script precisa do exemplo compilado (`cd t8flow && make 2d`). Ele imprime o que
gerou e quantos níveis cada arquivo tem. Regenerar uma medida troca o `.dat` e o
texto acompanha, sem a chance de a tabela e a prosa discordarem.

Para mudar os níveis de refino: `NIV=3,4,5,6 ./gera-dados.sh` (mais caro; a
ordem 4 do MLS no nível 6 leva vários minutos).

## O desenho é uma matriz esquema × malha

|                      | uniforme      | adaptada (2:1) |
|----------------------|---------------|----------------|
| fluxo 2ª ordem, 2 pt | grau 1        | grau 1         |
| fluxo 2ª ordem, 3 pt | grau 2        | grau 1         |
| compacto 4ª ordem    | grau 3        | **recusa**     |
| MLS + Gauss, ordem 3 | grau 3        | grau 3         |
| MLS + Gauss, ordem 4 | grau 4        | grau 3 (rebaixa) |

As dez células dessa matriz — incluindo a recusa, que é um resultado — são o que
o artigo traz, e cada uma está medida. A estrutura das seções segue do desenho em
vez de separar validação de resultados.

## Os blocos \status{} e \falta{} / \reserva{}

Cada seção abre com um `\status{}` dizendo que respaldo ela tem. Definir
`\nostatus` antes de `\begin{document}` para suprimi-los na versão de submissão.

- `\falta{...}` marca o que falta **medir**;
- `\reserva{...}` marca o que foi medido mas com ressalva que **muda a leitura**.

Confundir os dois esconderia o segundo tipo entre pendências, quando ele é um
resultado com reserva.

## O que o artigo declara que não traz

Seção 7.2, e não diluído entre os resultados: execução paralela, três dimensões,
Robin dentro do operador, e o fechamento de parede de ordem 5.
