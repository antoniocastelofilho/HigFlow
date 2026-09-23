# Relatório: fronteira imersa no HiGFlow

Projeto LaTeX compilável. No Overleaf: envie a pasta inteira; o arquivo
principal é `main.tex`.

```
main.tex            preâmbulo e ordem das seções
secoes/*.tex        uma seção por arquivo
referencias.bib     11 entradas, todas verificadas em fonte primária
dados/*.dat         AS MEDIDAS
figuras/            (vazio por enquanto)
```

## Os números vivem nos `.dat`, não no texto

As tabelas são geradas por `pgfplotstable` lendo `dados/*.dat`. Regenerar uma
medida é trocar o arquivo; o relatório acompanha. Isso evita o erro clássico de
a tabela e o texto discordarem depois de uma correção.

Consequência prática: **não edite números dentro do `.tex`.** Se um valor
estiver errado, ele está errado no `.dat`.

## O que ainda falta

`dados/refino-cd.dat` tem a linha da terceira resolução (80 células/D)
**comentada**, aguardando a corrida. Quando ela fechar:

1. descomentar a linha e preencher os três valores;
2. atualizar, na Seção 6, a razão entre os erros e a inclinação de convergência;
3. a Seção 7 tem um item "em curso" que passa a conclusão.

## Figuras

Três, todas geradas por `pgfplots` lendo os mesmos `.dat` das tabelas --- não há
imagem pré-renderizada para ficar desatualizada:

1. resíduo × `dt`, defasado contra pós-preditor (a rampa e o patamar);
2. resíduo × `dt` nas duas malhas (duas retas paralelas: refinar piora);
3. `Cd` × resolução, com a linha da referência (a troca de sinal).

## NÃO FOI COMPILADO

Não há LaTeX nesta máquina — este projeto **nunca passou por um compilador**.
A sintaxe foi escrita com cuidado e o risco mais provável foi removido (valores
como `1e104` estouram o parser numérico do TeX; aquelas colunas são `string
type`), mas isso não substitui compilar.

Compile no Overleaf antes de confiar. Se algo quebrar, os suspeitos, em ordem:
`pgfplotstable` lendo os `.dat`, as colunas `S` do `siunitx`, e os acentos nas
legendas de coluna.
