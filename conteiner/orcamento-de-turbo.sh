#!/bin/bash
# CONFIRMAR (ou derrubar) o orcamento de turbo como causa do padrao
# "primeira corrida rapida" que apareceu em QUATRO varreduras seguidas.
#
# A HIPOTESE.  O processador roda em limite de potencia curto enquanto o
# orcamento dura, depois cai para o sustentado.  Clock menor da' potencia
# menor, entao a temperatura CAI JUNTO com o desempenho.
#
# A PREVISAO QUE A DISTINGUE de estrangulamento termico: na hipotese termica,
# quente => lento, e a temperatura SOBE ate' o limite.  Aqui e' o contrario --
# as corridas rapidas sao as QUENTES, e a temperatura cai quando o desempenho
# cai.  Ja' observado nos tempos; aqui se mede a FREQUENCIA, que e' a grandeza
# que a hipotese afirma diretamente.
#
# ENERGIA RAPL NAO ENTRA: /sys/class/powercap/intel-rapl:0/energy_uj e'
# 0400 root:root nesta maquina (mitigacao PLATYPUS).  Sem sudo nao ha' leitura
# de potencia, e inventar uma estimativa a partir da frequencia seria medir a
# hipotese com ela mesma.  Fica a frequencia, que e' leitura direta.
CARGA=/home/castelo/t8flow-corridas/perfil/carga-fixa
TEMP=/sys/class/thermal/thermal_zone2/temp
N=${N:-14}

freq_media() {   # MHz, media sobre as 16 CPUs logicas
  # O \n NAO E' DETALHE: sem ele as amostras do laco saem concatenadas numa
  # linha so' e o awk de fora le um numero gigante.  Aconteceu na primeira
  # tentativa e a coluna inteira de frequencia foi perdida.
  awk '{s+=$1; n++} END {printf "%.0f\n", s/n/1000}' \
      /sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq
}
temp_c() { echo $(( $(cat $TEMP 2>/dev/null || echo 0) / 1000 )); }

echo "# repouso de 180 s para o orcamento se recompor"
echo "# antes do repouso: $(freq_media) MHz, $(temp_c) C"
sleep 180
echo "# depois do repouso: $(freq_media) MHz, $(temp_c) C"
echo
printf "%4s %9s %9s %9s %9s %7s %7s\n" "n" "tempo(s)" "freq med" "freq min" "freq max" "T ini" "T fim"
for i in $(seq 1 $N); do
  t0=$(temp_c)
  # amostra a frequencia durante a carga, nao antes nem depois
  ( for k in $(seq 1 60); do freq_media; sleep 0.5; done > /tmp/f.$$ ) &
  AMOSTRA=$!
  t=$(mpirun -np 8 $CARGA 2>/dev/null | cut -d' ' -f1)
  kill $AMOSTRA 2>/dev/null; wait $AMOSTRA 2>/dev/null
  f=$(awk '{s+=$1;n++} END {if(n)printf "%.0f", s/n; else print 0}' /tmp/f.$$)
  fmax=$(sort -n /tmp/f.$$ | tail -1); fmin=$(sort -n /tmp/f.$$ | head -1)
  rm -f /tmp/f.$$
  printf "%4d %9s %9s %9s %9s %7d %7d\n" "$i" "$t" "$f" "$fmin" "$fmax" "$t0" "$(temp_c)"
done
echo "FIM"
