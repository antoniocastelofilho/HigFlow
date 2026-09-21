# Schäfer–Turek 2D-1, com o cilindro por fronteira imersa

Escoamento em canal sobre cilindro, `Re = 20`, estacionário. É o benchmark com
valores publicados de `Cd`, `Cl` e `Δp` — e é por isso que ele existe aqui: o
oráculo é **número contra número**, não a figura da esteira, que parece certa
numa faixa larga de implementações erradas.

## A escala, e por que não é a do artigo

O original é dimensional: canal 2,2 × 0,41, cilindro `D = 0,1` centrado em
(0,2 ; 0,2), `ν = 1e-3`, entrada `u(y) = 4·Um·y(H−y)/H²` com `Um = 0,3`. Daí
`ū = 2Um/3 = 0,2` e `Re = ū·D/ν = 20`.

O HiGFlow aplica `1/Re` no termo viscoso, então **dar a geometria dimensional e
`Re = 20` produziria um Reynolds efetivo de 0,4**. Aqui tudo está
adimensionalizado por `D` e pela velocidade média:

```
canal      22 × 4,1          cilindro   D = 1, centro (2 ; 2)
entrada    u(y) = 6y(4,1−y)/4,1²        média 1, máximo 1,5
malha      440 × 82, h = 0,05           20 células no diâmetro
Re         20
```

`Cd` e `Cl` são adimensionais e não mudam com essa escala. `Δp` fica
adimensionalizado por `ρū²`.

**O centro em `y = 2` num canal `[0 ; 4,1]` está fora do eixo de propósito** — 2
do fundo, 2,1 do topo. A assimetria é do benchmark, e é ela que fixa o `Cl`
publicado. Centralizar invalida a comparação.

## O custo, medido

```
0,59 s por passo em np=3 (16 núcleos na máquina)
um atravessamento (t = 22) com dt = 1e-3:  22.000 passos, 3,6 h
o mesmo com dt = 1e-2:                                     0,4 h
```

O estado estacionário pede alguns atravessamentos. **O maior ganho disponível é
o passo de tempo**: o CFL convectivo permite até `~0,028`, então `dt = 1e-2` é
dez vezes menos trabalho. O que limita não é o CFL, é a fronteira imersa — a
força defasada impõe o não escorregamento a `O(dt)` e o ganho da realimentação
cresce com `dt`. Achar o maior `dt` estável **com o corpo ativo** é o primeiro
experimento, e vale por dez de qualquer outro.

## O que ainda não existe

**O cálculo do arrasto.** `fi_forca_total` já dá a força total sobre o corpo; o
`Cd` é menos a reação, adimensionalizado por `½ρū²D`. Conferir o sinal contra um
caso de sinal óbvio antes de confiar no número.

E os valores de referência do benchmark **não estão escritos aqui de propósito**:
devem ser lidos da fonte quando forem usados, não da memória de quem escreveu
este arquivo.
