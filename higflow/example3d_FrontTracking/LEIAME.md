# example3d_FrontTracking — fase B2 3D, gota estática (lei de Laplace)

Caixa fechada em repouso, uma esfera de front-tracking parada no centro, e a
tensão superficial impondo o salto de pressão. O oráculo é

    Δp = 2 σ / R

**Dois sobre R, não um.** Em 2D o salto de Laplace é σ/R; em 3D a curvatura
média da esfera é 2/R. O B2 bidimensional deste repositório fechou a 0,02%
contra σ/R, e copiar aquele número para cá validaria o número errado.

## Como rodar

```bash
python3 gera-malha.py 20 1.0     # malha 20^3 no cubo unitário
./roda-b2.sh 20 3                # N=20, icosfera com 3 subdivisões
```

O resultado sai numa linha só, no fim do log:

```
=+=+=+= LAPLACE 3D  p_in=...  p_out=...  Dp=...  2*sigma/R=...  erro_rel=...  |u|max=... =+=+=+=
```

`|u|max` é a **corrente espúria**: num equilíbrio exato a velocidade seria zero,
e o que sobra mede o erro do acoplamento.

## O par (N, nsub) não é livre

Front-tracking exige espaçamento de marcador **da ordem de h**. Refinar só a
malha deixa a frente grossa demais e o resultado **piora**: medido, N=30 com
`nsub=3` (ds/h = 1,13) deu erro maior que N=20 com o mesmo `nsub` (ds/h = 0,75).

Os pares que mantêm ds/h ≈ 0,75, com R = 0,25:

| N  | h      | nsub | ds     | ds/h |
|----|--------|------|--------|------|
| 20 | 0,0500 | 3    | 0,0377 | 0,75 |
| 40 | 0,0250 | 4    | 0,0188 | 0,75 |

## Três coisas que travaram isto, e que travam de novo

**1. As seis faces em Neumann deixam a pressão singular.** O sistema fica
determinado só a menos de constante, e o Krylov não converge para nada útil:
Δp saía ~0 com a força chegando certa à malha. Uma face (`bc5`) tem pressão
Dirichlet, que fixa o nível. O exemplo 2D que fecha o B2 faz igual.

**2. `explicit_euler` é instável aqui.** Com h = 0,05 e Re = 1, o número de
difusão é 0,2 por direção — **1,2 em três dimensões**, acima de 1. A velocidade
crescia até 23. Com `semi_implicit_euler` ela estabiliza em 1,8e−3.

**3. Os `.o` de `../src` são compartilhados entre TODOS os exemplos, e carregam
a dimensão.** Construir este exemplo recompila `../src/*.o` com `-DDIM=3`. Os
exemplos 2D passam então a linkar objetos 3D **em silêncio**, porque a regra
`%.o: %.c libhig$(DIM)d.a` só reconstrói se a biblioteca for mais nova que o
objeto — e depois de um build 3D ela não é.

Ao voltar para um exemplo 2D, force a reconstrução:

```bash
make -C .. clean && make -C .. DIM=2
```

Isto não é defeito deste exemplo; é propriedade da árvore, e vale para qualquer
alternância entre `example2d_*` e `example3d_*`.

## O que este exemplo NÃO faz

**Monofásico.** ρ e μ uniformes, de propósito: isola o termo de tensão do salto
de propriedade, que é o que esta fase existe para verificar. O bifásico precisa
da fração de volume de uma célula dentro de uma superfície triangulada, que é
trabalho da fase seguinte.

**A superfície fica parada** (`FT3_ADVECTA=0`). Advectar com a velocidade da
malha sem um oráculo de movimento seria exatamente o que os degraus anteriores
existem para evitar; o adaptador aborta se for pedido.

**Serial.** Em `np=1` não há franja e o espalhamento é completo. A regra de posse
que o 2D adotou (acumular só em faceta própria, com a frente replicada) se
transporta sem mudança quando o paralelo entrar.
