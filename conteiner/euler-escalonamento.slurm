#!/bin/bash
#==============================================================================
# euler-escalonamento.slurm -- submissao do estudo de escalonamento.
#
#   sbatch --partition=<A|B> conteiner/euler-escalonamento.slurm forte
#
# OS DOIS TIPOS DE NO' DO EULER exigem decisoes diferentes, e as duas estao
# marcadas abaixo.  Preencher ANTES de submeter: os valores aqui sao os que eu
# NAO pude confirmar daqui, e um palpite errado numa diretiva do Slurm nao
# falha -- a fila aceita, o job roda num no' que nao era o pretendido, e o
# numero sai sem dono.
#
#   no' A   256 nucleos, sem acelerador
#   no' B    64 nucleos + H100
#
# [E5] O NOME DAS PARTICOES nao e' conhecido daqui.  `sinfo -s` no no' de
#      submissao lista.  NAO submeta com o padrao abaixo sem conferir.
#SBATCH --job-name=t8flow-escal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64          # [E5] ajustar ao no': 256 no A, 64 no B
#SBATCH --time=04:00:00               # [E5] teto da fila
#SBATCH --output=escal-%j.out
#SBATCH --exclusive                   # NAO NEGOCIAVEL numa medida de tempo:
                                      # vizinho no mesmo no' mede o vizinho.

set -eu

# [E2] LANCAMENTO.  Dentro de Apptainer com Slurm, o caminho e' `srun
#      --mpi=pmix`; o PMIx do host tem de casar com o do OpenMPI 5 da imagem.
#      Se nao casar, o sintoma tipico e' cada rank achar que e' o rank 0 --
#      a corrida "funciona" e mede N vezes o problema serial.  A GUARDA e'
#      imprimir o tamanho do comunicador antes de qualquer medida.
IMAGEM="${IMAGEM:-$HOME/higflow-euler.sif}"
export LANCA="srun --mpi=pmix -n"

# [E4] MONTAGENS.  Entrada e saida fora de $HOME; $SCRATCH costuma ser o lugar.
#      Um VTK por processo em centenas de processos castiga sistema de arquivos
#      paralelo -- o protocolo nao grava VTK durante a medida, de proposito.
TRAB="${SCRATCH:-$PWD}/escal-$SLURM_JOB_ID"
mkdir -p "$TRAB"
cd "$TRAB"

echo "== ambiente =="
echo "no's:      $SLURM_JOB_NUM_NODES"
echo "tarefas:   $SLURM_NTASKS"
echo "particao:  ${SLURM_JOB_PARTITION:-?}"

# A GUARDA DO [E2], antes de qualquer tempo: se o comunicador tiver tamanho 1
# com N tarefas pedidas, o lancamento nao formou um job MPI e tudo o que vier
# depois e' N copias do caso serial -- que parecem rodar, e dao speedup 1,0
# com a forma de uma curva.
apptainer exec "$IMAGEM" $LANCA "$SLURM_NTASKS" \
  bash -c 'echo "rank $PMIX_RANK de ${SLURM_NTASKS}"' | sort | head -3
echo

MODO="${1:-forte}"
apptainer exec "$IMAGEM" \
  /opt/higflow/t8flow/conteiner/escalonamento-euler.sh "$MODO" \
  --np "${NPS:-1,2,4,8,16,32,64}" --reps "${REPS:-3}" --saida "$TRAB"

echo
echo "resultados em $TRAB"
