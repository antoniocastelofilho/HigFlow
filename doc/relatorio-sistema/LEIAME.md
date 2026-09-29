# Relatório de sistema: inventário do HiGTree/HiGFlow

Projeto LaTeX autossuficiente. **No Overleaf:** envie a pasta inteira; o arquivo
principal é `main.tex`, compilador pdfLaTeX, bibliografia por BibTeX.

Mesma convenção dos outros relatórios desta pasta: os números vivem em
`dados/*.dat` e são lidos por `pgfplotstable`, não digitados no texto.

## O que este relatório é

Um **inventário**: o que existe no código, em quantas linhas, qual o papel de
cada módulo. Os números vêm de contagem direta da árvore (`wc -l`, `ls`,
`git log`), não de estimativa.

**Não é uma avaliação.** A presença de um módulo não diz que ele esteja correto,
verificado ou em uso. A seção final lista explicitamente o que o inventário não
diz.

## Sobre a bibliografia, e o que a busca alcançou

As entradas foram conferidas no registro do editor via API do Crossref
(`api.crossref.org/works/<DOI>`), não contra memória.

**A pós-graduação está incompleta, e por um motivo conhecido.** A busca localizou
uma única dissertação da USP com o sistema no título (Godoi, ICMC/São Carlos,
2022). O Crossref indexa metadados — título, autores, resumo —, não texto
integral: uma tese que *use* o HiGFlow sem nomeá-lo no título não aparece.

Dado que o sistema tem módulos de milhares de linhas para tixotropia, bandas de
cisalhamento, suspensões, eletro-osmose e viscoelasticidade integral — cada um do
tipo que rende uma dissertação —, é quase certo que existam outras. Elas não
estão listadas porque **não foram encontradas**, não porque não existam.

Para completar: busca de texto integral no repositório de teses da USP
(`teses.usp.br`), que está fora do alcance das ferramentas usadas aqui.

## Como regenerar os números

Os dados vêm de comandos diretos sobre a árvore:

```bash
wc -l higtree/src/*.c | sort -rn      # módulos da HiGTree
wc -l higflow/src/*.c | sort -rn      # módulos do HiGFlow
ls -d higflow/example*                # os exemplos
git log --format='%an' | sort | uniq -c | sort -rn   # contribuidores
```

Os enums de `flowtype`, modelos e discretizações saem de
`higflow/src/hig-flow-kernel.h`.
