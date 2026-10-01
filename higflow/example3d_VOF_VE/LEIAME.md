# Caso mínimo viscoelástico 3D

## A pergunta

O caminho multifásico viscoelástico **liga** em 3D (ver `tests/vof-3d/`), mas
nunca foi executado. Este caso pergunta se ele **roda**.

## A montagem, e por que o oráculo é o equilíbrio

Caixa fechada em repouso, esfera de fase 1 no centro, fase 0 (ambiente)
viscoelástica com De = 1 e beta = 0,5. Sem gravidade, sem entrada, **sem tensão
superficial**: o estado inicial é um equilíbrio exato.

Partindo do repouso com o tensor de conformação na identidade — o equilíbrio do
Oldroyd-B — tudo tem de ficar ali. Qualquer velocidade ou desvio da identidade
mede a maquinaria, não a física.

É teste **fraco de física e forte de encanamento**: pega NaN, laço `[DIM][DIM]`
que supõe duas dimensões, componente *z* não inicializada, índice trocado — o que
pode ter sobrado num caminho nunca compilado em 3D.

## Resultado

```
VE3D EQUILIBRIO  |u|max=0.000e+00  pior |A-I|=0.000e+00  NaN(u)=0  NaN(A)=0
VE3D PASSOU: o equilibrio se manteve
```

Zero **exato**, em 101 passos, 4 096 células.

## Uma armadilha de configuração que custou uma rodada

A primeira execução **falhou** com `|u|max = 7,9e−3`. Não era defeito do código:
**a tensão superficial é ligada por padrão** quando `surf_tension` está ausente
do YAML (`hig-flow-io.c:8240`), e o próprio log avisava — *"Surface Tension is
added to the simulation (Using default value)"*. Com CSF numa interface curva,
repouso não é equilíbrio discreto, e correntes espúrias de ~1e−3 são o que esse
método produz.

O oráculo estava errado, não a maquinaria. A chave está agora explícita no YAML.

## O que isto NÃO prova

Que a física viscoelástica esteja correta. Preservar o equilíbrio trivial é
condição necessária, não suficiente. Um caso com **deformação real** — onde o
tensor de conformação sai da identidade — é o próximo degrau, e aí o oráculo
precisa ser outro.
