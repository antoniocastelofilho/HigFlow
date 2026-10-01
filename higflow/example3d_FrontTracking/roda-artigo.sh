#!/bin/bash
# BOLHA ASCENDENTE 3D PELO FRONT-TRACKING -- mesmo caso do artigo de referencia
# (Silva, Fernandes, Organista, Souza & Castelo, Mathematics 2023, 11(18), 3900,
# Tabela 3), para comparar FT contra VOF e contra a literatura.
#
#   $1 = celulas por unidade (50 reproduz a malha 50x50x100 do artigo)
#   $2 = numsteps (1501 da' t = 3.0 com dt = 2e-3)
#   $3 = caso, 1 ou 2
#   $4 = nucleo a fixar (opcional); sem ele nao fixa
#
# O QUE ESTA COMPARACAO ISOLA.  FT e VOF compartilham todo o resto -- rho(x),
# mu(x), momento de coeficiente variavel, projecao --, e o gancho do FT entra
# no MESMO ponto do passo multifasico onde o VOF advecta a fracao.  A unica
# coisa diferente e' como a interface e' representada: campo de fracao com
# reconstrucao PLIC contra superficie triangulada explicita.  Por isso a
# diferenca medida e' atribuivel a' representacao, e nao a outra coisa.
#
# DUAS DIFERENCAS DE CONFIGURACAO QUE SAO DO METODO, NAO DESCUIDO:
#   surf_tension: false  -- a forca vem da FRENTE, por integral de linha por
#     triangulo (f = sigma (b-a) x n), nao do CSF do VOF.  Por isso sigma entra
#     por FT3_SIGMA e nao pelo numero capilar.
#   a fracao inicial vem da SUPERFICIE (front_tracking_3d_fracvol_inicial), pela
#     forma fechada do corte plano-caixa aplicada a' triangulacao.  O VOF a
#     inicializa por subdivisao subcelular.  Cada um usa a sua maquinaria, o
#     que e' o certo -- mas os dois tem de nascer com a MESMA massa, e e' isso
#     que os dois numeros impressos no quadro 0 permitem conferir.
#
# nsub=4: 5120 triangulos, aresta 0.0188 = 0.94 h na malha de 50 -- a aresta do
# triangulo casada com a celula, que e' a condicao da grade espacial do modulo.
M=/home/castelo/HiGFlow-System/HigFlow
cd $M && set -a && . ./varsrc && set +a
cd $M/higflow/example3d_FrontTracking

N=${1:-50}; PASSOS=${2:-1501}; CASO=${3:-1}; NUC=${4:-}

python3 gera-malha.py $N 1.0 1.0 2.0
sed -i "s/^  numsteps: .*/  numsteps: ${PASSOS}/" input/gota.load.par.contr.yaml

if [ "$CASO" = "2" ]; then
    export FT3_RHO0=1000.0 FT3_RHO1=1.0 FT3_MU0=10.0 FT3_MU1=0.1 FT3_SIGMA=1.96
else
    export FT3_RHO0=1000.0 FT3_RHO1=100.0 FT3_MU0=10.0 FT3_MU1=1.0 FT3_SIGMA=24.5
fi
export FT3_R=0.25 FT3_CX=0.5 FT3_CY=0.5 FT3_CZ=0.5
export FT3_NSUB=4
export FT3_ADVECTA=1              # a superficie TEM de se mover: e' bolha subindo
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# FIXAR O NUCLEO.  `taskset` sozinho prende o mpirun, nao o solver: o Open MPI
# reatribui a afinidade dos filhos.  Precisa do --bind-to none para a mascara
# do taskset sobreviver ate' o processo que faz a conta.
if [ -n "$NUC" ]; then
    PRE="taskset -c $NUC"; BIND="--bind-to none"
else
    PRE=""; BIND=""
fi

$PRE /usr/bin/time -f "PAREDE %e s" mpirun $BIND -use-hwthread-cpus -n 1 ./ns-example-3d \
  input/gota.load output/g.save VTKS/g.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 \
  > /tmp/ft-c${CASO}-N${N}.log 2>&1
echo "FT caso=$CASO N=$N nucleo=${NUC:-livre} exit=$? em $(date +%H:%M)" >> /tmp/sobe-status
