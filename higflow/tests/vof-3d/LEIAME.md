# VOF 3D viscoelástico — o que está verificado, e o que não

## A pergunta

Uma versão anterior do código (junho/2026) trazia um exemplo de bolha ascendente
3D chamando `higflow_solver_step_multiphase_viscoelastic`. A dúvida era se a
versão atual da árvore, que **reestruturou** esse caminho, ficou incompleta para
três dimensões.

## O teste

`./testa-link-3d.sh` compila tudo com `-DDIM=3` e liga um binário que referencia
`higflow_solver_step_multiphase_viscoelastic`. Precisa de `PETSC_DIR` apontando
para um PETSc compilado.

**Ele compila os próprios objetos**, em `obj3d/`, em vez de reusar os de
`higflow/src`. Isso não é zelo: os `.o` compartilhados **carregam a dimensão**, e
reusá-los tornaria o teste refém de qual exemplo foi construído por último — uma
passagem poderia significar apenas que alguém acabou de construir um exemplo 3D,
e uma falha, que acabou de construir um 2D. Verificado nos dois estados da
árvore.

**Resultado: liga com zero símbolos não resolvidos.**

### O que isso prova, e o que não prova

| prova | não prova |
|---|---|
| nenhuma função falta em 3D | que o resultado esteja certo |
| o ramo `#if DIM == 3` resolve | que rode sem quebrar |

Falta **caso e verificação**, não código.

## O que a comparação com a versão de junho mostrou

**As fontes da árvore são a linhagem nova, não a velha.**

- `hig-flow-vof-plic-3D.c` recebeu quatro commits em 18/09/2026, um deles
  corrigindo que *"os ramos degenerados do PLIC 3D estavam presos ao eixo x"*. A
  escolha de ramo por cúbicas foi trocada por bissecção sobre forma fechada, e
  `73504ad` removeu as seis funções órfãs — que são as `solver_equation_*` que a
  versão de junho ainda tem.
- Os outros quatro arquivos VOF 3D (advecção, HF, mehta, curvatura) têm o
  **mesmo conjunto de funções** nas duas versões.
- `hig-flow-step-multiphase-viscoelastic.c` foi **reestruturado**: a versão de
  junho duplica a maquinaria por fase (`_phase0`/`_phase1`); a atual interpola os
  parâmetros pela fração e chama o solver constitutivo **compartilhado** — o
  mesmo que viscoelástico, elastoviscoplástico, viscosidade variável, *shear
  banding* e *shear thickening* usam. Só 4 de 43 funções coincidem.
- O arquivo atual **não tem guarda de dimensão nenhuma** e não chama função
  `_2D`: é DIM-genérico. O passo multifásico base tem ramo 3D de verdade
  (`mehta`, `HF_padrao_3D`, `distance_3D`).

**Copiar `src/` da versão de junho reverteria isso.** O `hig-flow-io.c` iria de
10 972 para 2 413 linhas, levando junto o instantâneo e o t8code.

## Uma advertência sobre o caso de junho

No arquivo de junho, a viscosidade polimérica está **fixada em zero no código**,
com os valores reais comentados ao lado:

```c
//real visc_p_0   = 9.53251434971;
real visc_p_0   = 0.0;
real visc_p_1   = 0.0;
//real visc_p_1   = 100.0;
```

Com isso `S0` e `S1` saem nulos e não há contribuição viscoelástica: **o caso 3D
daquela versão rodou newtoniano-newtoniano**. Não existe, de nenhum dos lados,
uma execução 3D que de fato exercite viscoelasticidade.

A versão atual não fixa nada: lê `De`, `beta` e as viscosidades de
`ns->ed.mult.ve.par0/par1` e interpola pela fração.

## O que vem depois

Portar o caso `examples3D_rising_NN` — bolha ascendente 3D — que é o que a árvore
não tem. O porte tem três partes: entradas do formato antigo para YAML, a API de
oito ponteiros (`higflow_set_external_functions`) para o objeto de problema
(`higflow_set_problem`), e o Makefile 3D. E, ao portar, **usar viscosidade
polimérica não nula** — senão o caso repete o newtoniano de junho com outro nome.

## Nota de construção

`hig-flow-vof-9-cells.c` e `hig-flow-vof-elvira.c` **não estão** na lista
`MODULES` do Makefile de `higflow`: são compilados por exemplo. O script os
compila à mão.
