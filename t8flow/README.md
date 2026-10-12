# t8Flow

Sistema sobre o **t8code**, com o escalonamento dele (floresta como estrutura,
partição pela curva de preenchimento, nenhum processo segurando a malha inteira)
e **multiblocos** como o HiGFlow — o `cmesh` do t8code é nativamente multiárvore,
e é o mesmo conceito.

O projeto está em [`doc/projeto-t8flow.md`](../doc/projeto-t8flow.md): as cinco
decisões tomadas, o que é medido e o que é opinião. Leia-o antes de mexer aqui.

## As decisões, em uma linha cada

1. **Núcleo herdado** — solver e MLS vêm da HiGTree, porque já não sabem o que é
   uma malha (zero referências a `hig_cell`/`sim_domain` em todos os backends).
2. **Partição e domínio do t8code** — a floresta é a estrutura; sai o `lbal`.
3. **A suíte nasce com o código** — e cobra um contrato.
4. **Um domínio só** — sem o segundo domínio uniforme no nível mais fino, que é
   a malha cheia por construção e não encolhe com o AMR.
5. **C++ na norma do t8code** (`gnu++20`), sem `extern "C"` e sem `-fpermissive`.
6. **Vive em worktree**, branch `t8flow`.
7. **Código separado para 2D e 3D**, começando pelo 2D — o t8code tem malhas 2D e
   3D nativas (`T8_ECLASS_QUAD` e `T8_ECLASS_HEX` são classes distintas), então o
   código novo é escrito duas vezes em vez de uma vez com `#if DIM`. O núcleo
   herdado continua por macro, porque `Point` **é** `real[DIM]`.

## Construir e rodar

```bash
python3 ci/run_t8flow_tests.py
```

O driver constrói o lado 2D, roda os testes e cobra o contrato. Os prefixos de
dependência têm padrão apontando para a árvore principal e podem ser trocados por
`--t8code`, `--higtree-src` e `--petsc`.

À mão, a partir deste diretório:

```bash
make 2d && mpirun -n 1 ./tests/2d/test-mls2d
```

## CPU ou GPU: a escolha é de quem compila

O lote de estênceis (`src/2d/lote2d.hxx`) tem dois back-ends, e a flag `GPU`
decide:

```bash
make 2d            # detecta OpenCL; se houver, a GPU fica DISPONÍVEL, padrão CPU
make 2d GPU=sim    # exige OpenCL e faz a GPU ser o PADRÃO do lote
make 2d GPU=nao    # CPU só, mesmo que haja OpenCL instalado
```

São três estados e não dois: "detectei e deixei disponível" é diferente de "o
usuário pediu", e tratar os dois igual tira a escolha de quem compila. Com
`GPU=sim` sem OpenCL, o build **para** com a causa — quem pediu placa quer saber
que não há, não receber CPU em silêncio.

**A versão de CPU é de primeira classe**, não plano B: funciona sem OpenCL
instalado, sem OpenCL compilado e sem placa, e é ela que define o resultado
correto. Quem não tem GPU tem a suíte inteira verde — os casos de placa se
declaram `SALTADO` na saída, e o teste imprime `gpu-exercitada sim|nao` para o
verde nunca significar "a GPU foi testada" quando não foi.

**A configuração entra na dependência do make.** Trocar `GPU=sim` por `GPU=nao`
reconstrói sozinho: a configuração é um selo em `obj/.config`, conferido no parse.
Sem isso, o binário do build anterior ficava no lugar anunciando a configuração
errada — aconteceu aqui, e é a mesma classe de "build sem o que se quer medir"
que já custou duas corridas medidas como t8code que eram malha uniforme.

**Para a Radeon, o rusticl precisa de `RUSTICL_ENABLE=radeonsi`** — sem a
variável a plataforma aparece com zero dispositivos, que se lê como driver
quebrado em vez de desabilitado.

Estrutura, pela decisão 7:

    src/comum/   sem dimensão: arranque, contrato
    src/2d/      o discretizador (mls2d), o solver (solver2d), a malha
                 multiblocos com refino local e numeracao global (malha2d)
    src/3d/      quando existir
    tests/2d/

`make 3d` **diz** que o lado 3D ainda não existe, em vez de falhar obscuro.

**As dependências não estão neste worktree.** O t8code e o PETSc vivem em
`bibliotecas/`, que o git ignora; por isso os caminhos são parâmetros absolutos,
e não caminhos relativos. Caminho relativo repassado a um `make -C` resolve a
partir do diretório errado e o sintoma é "não construído" — falha silenciosa.

## A ordem de arranque não é livre

Descoberta com um SIGSEGV dentro do PETSc, sem uma linha que apontasse para o
nosso código:

1. `higtree_initialize(&argc, &argv)` — faz o `MPI_Init` **e guarda argc/argv**;
2. `sc_init` + `t8_init`, com o MPI já inicializado;
3. o PETSc entra preguiçoso, pelo próprio solver.

O motivo: `_slv_create_petsc` chama `_try_initialize_petsc` (`utils.c:34`), que
desreferencia um `static char **_argv` preenchido **só** por
`higtree_initialize`. E a flag `is_petsc_initialized` é `static`: inicializar o
PETSc por fora não a desarma. Não chame `sc_MPI_Init` (dois `MPI_Init` são erro)
nem `MPI_Finalize`/`PetscFinalize` (o `atexit` do `higtree_initialize` fecha os
dois).

Todo binário vai por `mpirun`, inclusive em np=1: rodado direto, trava no
`MPI_Init` sem imprimir nada e aparece como timeout, não como falha.

## O contrato

`src/t8flow-contract.h`. O driver **lê as cláusulas do arquivo** e exige, para
cada uma, um caso que saiu `PASS`. Cláusula sem `COBERTA-POR`, com caso ausente
ou com caso reprovado conta como falha da suíte — não como omissão. Conferido no
sentido inverso: apontando uma cláusula para um caso inexistente, o driver sai
com `CASO AUSENTE` e código 1.

Uma cláusula entra no contrato **quando tem caso que a afirme**, nunca antes. As
previstas e ainda sem código vivem no documento de projeto; se entrassem aqui, a
suíte nasceria vermelha e o vermelho deixaria de significar "quebrou".

## O que o verde não prova

Medido por sabotagem controlada em 11/10/2026, e esta seção é a mais útil do
arquivo:

- **deslocando todos os centroides do suporte meia célula, os três casos de F2
  continuam verdes.** A reprodução polinomial verifica a *discretização*, não a
  malha. Quem verificaria a malha é "o suporte cai na grade e é local" — será a
  cláusula F3, e ainda não existe;
- o ponto de consulta **no centro de simetria** da nuvem fez o ajuste de ordem 2
  reproduzir um campo cúbico com erro 1,1e-16 (os termos ímpares cancelam contra
  pesos simétricos). O caso de contraste reprovou dizendo "o oráculo não
  discrimina" — que é o seu trabalho. Os pontos de consulta são assimétricos de
  propósito;
- nenhuma física roda aqui, np vai só até 1, não há partição, nem franja, nem
  domínio. Tudo isso é o trabalho adiante;
- **o refino local não tem a guarda contra árvore de uma célula** — uma árvore de
  uma célula não tem interior, toda faceta dela é de fronteira, e a franja
  declarada do HiGFlow é de cinco. A guarda entra quando a franja entrar;
- **"centroide dentro da caixa" não pega espelhamento** — medido: trocar a ordem
  dos vértices do QUAD reflete o mapeamento dentro da *mesma* caixa, e todo
  centroide continua dentro dela. Quem acusa é a ordem de percurso (a folha 1
  deslocada em X da folha 0). O comentário do teste afirmava o contrário até a
  sabotagem ser rodada;
- **diferença exatamente zero entre dois back-ends não é boa notícia** — medido:
  com `GPU=sim` o padrão é a placa, e o vetor de referência que o teste chamava de
  "cpu" foi calculado na GPU; a comparação deu 0,000e+00 porque comparava GPU com
  GPU. Hoje o caso força a CPU e **recusa** diferença zero;
- **a contagem de junções não prova que houve junção** — medido: com o
  `t8_cmesh_set_join` removido, o contador continua certo e só o conteúdo (face
  entre blocos tem vizinho) acusa;
- o motor do solver **instalado** não é consultável — não existe
  `slv_get_engine`. N5 cobra a tabela de nomes, não a instalação.

Controles que confirmam o contrário — cada guarda foi removida e o teste
reprovou:

    ordem do ajuste 2 -> 1            N1 reprova nas duas metades
    guarda de suporte minimo          N3 reprova nomeando os numeros
    guarda de linha faltando          o processo morre no abort do herdado, sem
                                        linha de resumo

A última é a que justifica a fachada: ela converte um abort três camadas abaixo em
exceção que nomeia a linha que falta.

Na numeração global, os controles são:

    id_global sem o offset da árvore   G2 reprova: 4 ids repetidos, e o processo
                                         MORRE na guarda de linha repetida do
                                         solver — sem linha de resumo
    id_global deslocado de 1           G2 reprova: 1 id fora de [0, 20)
    id_do_ponto com intervalo fechado  G4 reprova: o desempate passa a depender
                                         da ordem de visita

E o limite da G5, que é tentador ler como mais do que é: ela pega **divergência
entre as duas vias**, não erro consistente na numeração. Com todos os ids
deslocados de 1 ela continua verde — escrita e leitura erram igual. Quem pega erro
consistente é a bijeção (G2) e a aritmética de Morton (G3).

No refino local, os controles são:

    t8_forest_set_balance removido   R3 reprova com salto de 2 níveis — e SÓ ela
    adapt não recursivo              R1 reprova: nível máximo 3 em vez de 4

E dois testes meus estavam errados antes de o código estar: "toda folha de nível
máximo tem o centro dentro da caixa" é falso (o critério julga a célula *antes* de
dividir, e os filhos herdam a decisão do pai), e um `unref` a mais deixou a malha
**vazia** — sem erro, contagens em zero. Quem pegou foi a conservação de área.
