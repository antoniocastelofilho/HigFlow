# Contêiner do HiGFlow para o cluster Euler (ICMC)

| arquivo | o que é |
|---|---|
| `higflow-euler.def` | a definição Apptainer — a receita |
| `preparar-fonte.sh` | junta o que `%files` precisa, por `git archive` de um commit |
| `inventario.sh` | o que está dentro da imagem (é o `%runscript` padrão) |
| `verificar-euler.sh` | o roteiro de verificação, com os portões e as guardas de vacuidade |
| `Dockerfile`, `Dockerfile.ci`, `Dockerfile.petsc` | **DEFASADOS** — Ubuntu 22.04 com PETSc 3.14, sem t8code. Mantidos por referência; o que vale neles são os ganchos `/.singularity.d/env/` |

## Uso

```
conteiner/preparar-fonte.sh                 # exige árvore limpa e commitada
apptainer build --fakeroot higflow-euler.sif conteiner/higflow-euler.def
conteiner/verificar-euler.sh higflow-euler.sif          # aqui
conteiner/verificar-euler.sh higflow-euler.sif --euler  # no nó de submissão
```

Sem `--fakeroot` nem root, construa por Docker e converta:
`apptainer build higflow-euler.sif docker-daemon://higflow:euler`.

## Quatro pontos que dependem do Euler

Marcados `[E1]`..`[E4]` na definição. **Confirmar antes de construir:**

1. **`[E1]` ISA dos nós** — `TARGET_ARCH`. Nunca `-march=native`: a imagem é
   construída aqui e roda lá; instrução inexistente dá *illegal instruction* em
   tempo de execução, e o sintoma não aponta para a causa.
2. **`[E2]` PMIx/Slurm** — o lançamento híbrido (`srun --mpi=pmix`) exige PMIx
   compatível com o do OpenMPI 5 da imagem.
3. **`[E3]` Interconexão** — InfiniBand/Omni-Path precisa de verbs/UCX casando
   com os drivers do host, senão cai para TCP.
4. **`[E4]` Montagens** — onde ficam entrada e saída (`$SCRATCH`).

## Duas diferenças deliberadas do ambiente de desenvolvimento

1. **PETSc otimizado.** O ambiente local usa `--with-debugging=1`; num cluster
   isso custa 2 a 5×. A imagem usa `--with-debugging=0 -O3`, e o portão 1 da
   verificação **confere** — a troca esquecida é silenciosa.
2. **Um único libfyaml.** Na máquina há dois (apt e `/usr/local`); com os dois,
   qual resolve depende da ordem do `ld.so`.

O inventário detecta as duas condições. Rodado contra a árvore local, ele
**reprova** — de propósito: é o teste positivo do detector.

## Três armadilhas que custaram tempo e estão fechadas no código

- **`mpirun --version` não serve de sonda** neste empilhamento: no OpenMPI
  5.0.10 do Ubuntu 26.04 ele não acha `help-schizo-ompi.txt`, imprime erro de
  ajuda em vez da versão e **sai com 0**. Use `ompi_info --version`.
- **`git archive` omite o não rastreado em silêncio.** Na primeira execução do
  `preparar-fonte.sh` os próprios arquivos desta pasta estavam não rastreados e
  o tarball saiu sem eles. A guarda agora recusa por **fonte** não rastreada e
  apenas conta o resto (malhas geradas), para não virar ruído.
- **No build in-place do PETSc os cabeçalhos não moram juntos**: `petscconf.h`
  em `$PETSC_DIR/include`, `petscversion.h` na fonte (`$PETSC_DIR/../include`).
