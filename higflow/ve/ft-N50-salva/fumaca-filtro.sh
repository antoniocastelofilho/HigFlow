#!/bin/bash
# Prova que o filtro esta' LIGADO NA CORRIDA -- coisa diferente de provar que o
# algoritmo esta' certo.  Os oraculos de higflow/tests/front-tracking/teste-tsur
# provam o algoritmo (esfera intocada, ondulacao -80%, volume 0,000%, laplaciano
# perdendo 9,74%, e o elipsoide que expoe a tensao artificial).  Aqui se verifica
# a fiacao: a cadencia chega do ambiente, o bloco e' alcancado no laco, a
# topologia sobrevive, e o ponto de salvamento inclui a superficie.
#
# Roda em N=12 e em nucleo separado, para custar minutos e nao competir.
set -u
M=/home/castelo/HiGFlow-System/HigFlow
T=${FUMACA_DIR:-/tmp/ft-fumaca-filtro}
NUC=${FUMACA_NUC:-3}
PASSOS=${FUMACA_PASSOS:-45}
CADA=${FUMACA_CADA:-20}
LIM=${FUMACA_LIM:-0.05}
falhas=0
ok()   { echo "  [ok]    $*"; }
nao()  { echo "  [FALHA] $*"; falhas=$((falhas+1)); }

echo "=== 1. os oraculos do algoritmo ==="
if ( cd $M/higflow/tests/front-tracking && ./teste-tsur >/tmp/fumaca-oraculos.log 2>&1 ); then
    ok "teste-tsur passou  ($(LC_ALL=C grep -ac '\[ok\]' /tmp/fumaca-oraculos.log) assercoes)"
else
    nao "teste-tsur FALHOU"; tail -5 /tmp/fumaca-oraculos.log | sed 's/^/          /'
fi

echo "=== 2. a fiacao, em N=12 com $PASSOS passos ==="
rm -rf $T && mkdir -p $T/{output,VTKS,DATA,amrs}
cp $M/higflow/ve/ft-N50-salva/gera-malha.py $T/
cp -r $M/higflow/ve/ft-N50-salva/input $T/
cp $M/higflow/example3d_FrontTracking/ns-example-3d $T/
cd $T
python3 gera-malha.py 12 1.0 1.0 3.0 >/dev/null
ln -sfn ../input/gota.load.domain.yaml output/o.save.domain.yaml
ln -sfn ../input/gota.load.bc.yaml     output/o.save.bc.yaml
sed -i "s/^  numsteps: .*/  numsteps: $PASSOS/;s/^  dts: .*/  dts: 0.05/;s/^  dtp: .*/  dtp: 0.05/" input/gota.load.par.contr.yaml
cd $M && set -a && . ./varsrc && set +a
cd $T
export FT3_RHO0=1000.0 FT3_RHO1=100.0 FT3_MU0=10.0 FT3_MU1=1.0 FT3_SIGMA=24.5
export FT3_R=0.25 FT3_CX=0.5 FT3_CY=0.5 FT3_CZ=0.5 FT3_NSUB=4 FT3_ADVECTA=1
export FT3_SALVA=1 FT3_CAMPO_CADA=0
export FT3_TSUR_CADA=$CADA FT3_TSUR_LIM=$LIM FT3_TSUR_RLX=0.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
taskset -c $NUC stdbuf -oL -eL mpirun --bind-to none -use-hwthread-cpus -n 1 \
  ./ns-example-3d input/gota.load output/o.save VTKS/p.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > log 2>&1
R=$?
[ $R -eq 0 ] && ok "a corrida terminou (exit 0)" || nao "a corrida saiu com $R"

# (a) a cadencia e o limiar chegaram do ambiente
if LC_ALL=C grep -aq "BALANCE+TSUR a cada $CADA passo(s)" log && \
   LC_ALL=C grep -aq "limiar=$LIM" log; then
    ok "cadencia e limiar vieram do ambiente: $(LC_ALL=C grep -ao 'BALANCE+TSUR.*' log | head -1)"
else
    nao "o anuncio do filtro nao reflete o ambiente (CADA=$CADA LIM=$LIM)"
    LC_ALL=C grep -a 'Ondulacoes submalha' log | sed 's/^/          /'
fi

# (b) o bloco foi ALCANCADO no laco -- e' isto que distingue ligado de nao ligado
n=$(LC_ALL=C grep -ac '===> TSUR passo' log)
esp=$(( PASSOS / CADA ))
if [ "$n" -ge 1 ] && [ "$n" -ge "$esp" ]; then
    ok "o filtro foi aplicado $n vez(es) em $PASSOS passos (esperado >= $esp)"
    LC_ALL=C grep -a '===> TSUR passo' log | sed 's/^/          /'
else
    nao "o filtro NAO foi aplicado ($n vezes, esperado >= $esp) -- fiacao morta"
fi

# (c) topologia intacta em TODA aplicacao
mal=$(LC_ALL=C grep -a '===> TSUR passo' log | LC_ALL=C grep -avc 'chi=2')
[ "$mal" -eq 0 ] && ok "chi=2 em todas as $n aplicacoes" \
                 || nao "$mal aplicacao(oes) com chi != 2"

# (d) o ponto de salvamento inclui a SUPERFICIE
P=$(awk '/step:/{print $2}' output/o.save.init.yaml 2>/dev/null)
O=$(awk '/^passo/{print $2}' output/o.save.ok 2>/dev/null)
if [ -n "${P:-}" ] && [ "$P" = "${O:-x}" ] && [ -s output/o.save.frente ] && \
   [ ! -f output/o.save.gravando ]; then
    ok "salvamento coerente: passo $P, frente com $(stat -c%s output/o.save.frente) bytes, sem sinal pendente"
else
    nao "salvamento incompleto: init=$P marca=$O frente=$(stat -c%s output/o.save.frente 2>/dev/null || echo ausente)"
fi

# (d2) GUARDA DE VACUIDADE, e e' a assercao que faz este teste poder falhar.
# Acima o filtro foi aplicado e moveu ZERO vertices -- correto, porque em 45
# passos a bolha e' quase esferica e o limiar bloqueia.  Mas "bloqueou certo" e
# "fiacao morta que nunca move nada" produzem a MESMA saida, e um teste que nao
# os distingue nao testa.  Entao: repete com limiar 0, onde o filtro TEM de
# mover, no MESMO binario de producao.  Se mover zero ali, esta' morto.
echo "=== 3. guarda de vacuidade: limiar 0 TEM de mover ==="
rm -rf $T/vac && mkdir -p $T/vac/{output,VTKS,amrs}
cp -r $T/input $T/vac/; cp $T/ns-example-3d $T/vac/; cp -r $T/amrs/. $T/vac/amrs/
cd $T/vac
ln -sfn ../input/gota.load.domain.yaml output/o.save.domain.yaml
ln -sfn ../input/gota.load.bc.yaml     output/o.save.bc.yaml
sed -i "s/^  numsteps: .*/  numsteps: 11/" input/gota.load.par.contr.yaml
FT3_TSUR_CADA=10 FT3_TSUR_LIM=0.0 FT3_SALVA=0 \
taskset -c $NUC stdbuf -oL -eL mpirun --bind-to none -use-hwthread-cpus -n 1 \
  ./ns-example-3d input/gota.load output/o.save VTKS/p.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > log 2>&1
lin=$(LC_ALL=C grep -a '===> TSUR passo' log | head -1)
mov=$(echo "$lin" | sed -n 's/.*BALANCE moveu \([0-9]*\), TSUR \([0-9]*\).*/\1 \2/p')
a=$(echo $mov | cut -d' ' -f1); b=$(echo $mov | cut -d' ' -f2)
if [ -n "${a:-}" ] && [ $(( ${a:-0} + ${b:-0} )) -gt 0 ]; then
    ok "com limiar 0 o filtro moveu ${a:-0}+${b:-0} vertices -- a fiacao esta' VIVA"
    echo "          $lin"
else
    nao "com limiar 0 o filtro moveu ZERO -- fiacao MORTA no binario de producao"
    echo "          ${lin:-nenhuma linha TSUR}"
fi
cd $T

# (e) a serie de quadros existe, que e' o que a animacao consome
q=$(ls VTKS/*frente_* 2>/dev/null | wc -l)
[ "$q" -ge 1 ] && ok "$q quadro(s) de superficie para a animacao" \
               || nao "nenhum quadro de superficie"

echo
[ $falhas -eq 0 ] && { echo "FUMACA DO FILTRO: PASSOU"; exit 0; } \
                  || { echo "FUMACA DO FILTRO: $falhas FALHA(S)"; exit 1; }
