# As premissas de desempenho do projeto CPU/GPU, contra a medida

Projeto LaTeX autossuficiente. **No Overleaf:** envie a pasta inteira; principal
`main.tex`, pdfLaTeX, BibTeX.

Este documento foi separado do `relatorio-sistema` de propósito: aquele é um
**inventário do que existe**, este **discute um projeto futuro e faz
recomendações**. São gêneros diferentes, e misturá-los faria o inventário parecer
uma tomada de posição.

## O que ele não é

Não é uma avaliação do projeto. O núcleo compartilhado entre HiGFlow e MFSim-NG é
a parte principal dele e independe inteiramente do que está aqui — unificar duas
implementações divergentes tem valor próprio, de manutenção e de sustentabilidade,
que nenhuma medida de desempenho toca.

O que este documento oferece é a contribuição de uma medida para a **ordem** das
extensões.

## As medidas, e de onde vêm

Duas campanhas no HiGFlow, sobre o caso de bolha ascendente bifásica (~9.400
células): o custo do remalhamento cronometrado separado do resto em np = 1, 2, 4,
8; e a contagem de iterações do KSP nos mesmos pontos.

**A pequenez do caso é a principal limitação**, e está dita no documento. O
projeto mira malhas grandes, onde o equilíbrio muda.
