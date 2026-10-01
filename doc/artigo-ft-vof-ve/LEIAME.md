# Artigo: FT contra VOF em 3D, newtoniano e viscoelastico, malha uniforme e adaptativa

Esqueleto do artigo. Compila com

    pdflatex main && bibtex main && pdflatex main && pdflatex main

ou `latexmk -pdf main`. No Overleaf: enviar a pasta inteira; o arquivo
principal e' `main.tex`.

## O desenho e' um 2x2x2

Dois metodos (front-tracking, VOF) x duas reologias (newtoniana, viscoelastica)
x duas malhas (uniforme, adaptativa), tudo em tres dimensoes. As oito celulas
dessa matriz sao as oito corridas de producao que o artigo tem de trazer, e a
estrutura das secoes segue do desenho em vez de separar validacao de resultados.

## Os numeros vivem em dados/

Mesma convencao dos relatorios desta arvore: nenhum numero e' digitado no
texto. Regenerar uma medida troca o `.dat` e o texto acompanha, sem a chance de
a tabela e a prosa discordarem. A serie temporal sai de
`higflow/examples-common/metricas-superficie.py`, que le os VTK de superficie
dos DOIS metodos -- o do marching cubes no VOF e o da triangulacao no FT -- de
modo que a comparacao nao herde diferenca de instrumento.

Campos matematicos dentro de um `.dat` precisam de `$...$`: o pgfplotstable os
compoe literalmente, e um `_` solto para o LaTeX com "Missing $ inserted".

## Dois marcadores, visiveis de proposito

- `\status{...}` abre cada secao dizendo o que ela tem de lastro medido. Defina
  `\nostatus` antes do `\begin{document}` para suprimi-los.
- `\missing{...}` marca uma lacuna pontual dentro do texto.

Eles existem para que ninguem -- inclusive nos -- leia um paragrafo de
planejamento como resultado.

## Na submissao

- definir `\nostatus`
- remover o `\nocite{*}`, que esta' la' so' para a bibliografia do esqueleto
  listar as entradas disponiveis
- preencher autores, data e resumo

## Armadilha do BibTeX, registrada porque ja' custou duas vezes

O BibTeX NAO trata `%` como comentario no nivel de topo. Qualquer arroba fora de
uma entrada e' lida como inicio de entrada, mesmo numa linha que comeca com
`%` -- e o sintoma e' bibliografia vazia com o PDF compilando sem reclamar.
Para desativar uma entrada, apague-a.
