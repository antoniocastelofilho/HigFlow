#!/bin/bash
# Acha o maior dt ESTAVEL com o corpo ativo, no caso 2D.
#
# CRITERIO: numero FIXO de passos, nao tempo fixo.  Instabilidade de
# realimentacao cresce por PASSO -- cada passo e' uma oportunidade de
# amplificacao -- entao comparar a tempo fixo daria a cada dt um numero
# diferente de oportunidades, e o dt grande pareceria melhor do que e'.
#
# O veredito e' o maior |u| ao fim.  A entrada tem maximo 1,5; qualquer coisa
# acima de poucas unidades e' crescimento, nao fisica.
set -u
cd /home/castelo/HiGFlow-System/HigFlow/higflow/example2d_SchaeferTurek
S="${1:-/tmp}"; PASSOS=${PASSOS:-400}
cp input/example-2d.load.par.contr.yaml "$S/par2.bak"
printf "  %-10s %-8s %-14s %s\n" "dt" "passos" "max|u| final" "veredito"
for DT in ${DTS:-0.001 0.002 0.005 0.01 0.02 0.05 0.1}; do
  sed -e "s/dt: .*#/dt: $DT                       #/" \
      -e "s/numsteps: .*#/numsteps: $PASSOS                         #/" \
      -e "s/dtp: .*#/dtp: 1000.0                         #/" \
      -e "s/dts: .*#/dts: 1000.0                         #/" \
      "$S/par2.bak" > input/example-2d.load.par.contr.yaml
  timeout 1800 mpirun --mca btl_base_warn_component_unused 0 -use-hwthread-cpus -n 3 \
      ./ns-example input/example-2d.load "$S/v.save" "$S/v.print" > "$S/dt-$DT.log" 2>&1
  saida=$?
  vmax=$(grep -E "Vmin" "$S/dt-$DT.log" | tail -2 | grep -oE "Vmax *= *[-0-9.]+" | grep -oE "[-0-9.]+" | sort -g | tail -1)
  res=$(grep "RESIDUO" "$S/dt-$DT.log" | tail -1 | grep -oE "[0-9.]+e[-+][0-9]+" | head -1)
  if [ "$saida" != 0 ]; then v="ABORTOU"
  elif [ -z "$vmax" ]; then v="sem saida"
  else v=$(python3 -c "import sys; x=float('$vmax'); print('estavel' if abs(x)<5 else 'DIVERGIU')")
  fi
  printf "  %-10s %-8s %-14s %s  residuo %s\n" "$DT" "$PASSOS" "${vmax:-?}" "$v" "${res:-?}"
done
cp "$S/par2.bak" input/example-2d.load.par.contr.yaml
