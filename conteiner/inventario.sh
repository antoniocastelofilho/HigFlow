#!/usr/bin/env bash
#==============================================================================
# inventario.sh -- o que esta' DENTRO da imagem.  E' o %runscript padrao.
#
# Serve para conferir a imagem num no' do cluster antes de gastar uma alocacao:
# roda em segundos, nao exige MPI entre nos, e cada linha e' um fato que o
# roteiro de verificacao depois testa de verdade.
#==============================================================================
set -u
: "${PETSC_DIR:=/opt/petsc-3.25.4/x86_64}"
: "${T8CODE:=/opt/t8code}"
: "${HIGFLOW_RAIZ:=/opt/HigFlow}"

echo "=============================================================="
echo " HiGFlow -- inventario da imagem"
echo "=============================================================="
[ -f $HIGFLOW_RAIZ/conteiner/_build/PROCEDENCIA ] && \
  sed 's/^/  /' $HIGFLOW_RAIZ/conteiner/_build/PROCEDENCIA
[ -f /.procedencia ] && sed 's/^/  /' /.procedencia

echo "-- base"
. /etc/os-release 2>/dev/null && echo "   $PRETTY_NAME"
echo "   glibc $(ldd --version | head -1 | awk '{print $NF}')"
echo "   gcc   $(gcc -dumpversion 2>/dev/null)"

echo "-- MPI"
# `mpirun --version` NAO serve de sonda neste empilhamento: no OpenMPI 5.0.10 do
# Ubuntu 26.04 ele nao acha help-schizo-ompi.txt, imprime um erro de ajuda em vez
# da versao -- e SAI COM 0.  Uma checagem baseada nele passa sempre, dizendo nada.
# `ompi_info --version` responde certo.
if command -v ompi_info >/dev/null 2>&1; then
  echo "   $(ompi_info --version 2>/dev/null | head -1)"
else
  V=$(mpirun --version 2>/dev/null | head -1)
  case "$V" in *Sorry*|"") echo "   *** nao consegui a versao do MPI ***" ;; *) echo "   $V" ;; esac
fi
echo "   mpicc: $(mpicc -show 2>/dev/null | cut -c1-80)"

echo "-- PETSc"
# No build IN-PLACE os dois cabecalhos que interessam NAO moram juntos:
#   petscconf.h   -> $PETSC_DIR/include        (gerado pelo configure)
#   petscversion.h-> $PETSC_DIR/../include     (fica na fonte)
# A primeira versao disto procurava a versao em $PETSC_DIR/include e concluia
# "AUSENTE" num PETSc perfeitamente instalado -- alarme falso que mascararia o
# portao seguinte, que e' o que importa.
PCONF="$PETSC_DIR/include/petscconf.h"
PVER=""
for c in "$PETSC_DIR/include/petscversion.h" "$PETSC_DIR/../include/petscversion.h"; do
  [ -f "$c" ] && { PVER="$c"; break; }
done
if [ -f "$PCONF" ]; then
  if [ -n "$PVER" ]; then
    # ANCORADO em ^#define: o petscversion.h usa PETSC_VERSION_MAJOR dentro de
    # macros de comparacao (linhas 23, 32, 33), e sem a ancora o awk colhia
    # "(MAJOR)" dali -- a versao saia "(MAJOR).(MINOR).(SUBMINOR)".
    V=$(awk '/^#define PETSC_VERSION_MAJOR/{a=$3}
             /^#define PETSC_VERSION_MINOR/{b=$3}
             /^#define PETSC_VERSION_SUBMINOR/{c=$3}
             END{print a"."b"."c}' "$PVER")
  else
    V="(petscversion.h nao achado)"
  fi
  echo "   versao $V   prefixo $PETSC_DIR"
  if grep -q 'define PETSC_USE_DEBUG' "$PCONF" 2>/dev/null; then
    echo "   *** BUILD DE DEPURACAO -- 2 a 5x mais lento; NAO usar em producao ***"
  else
    echo "   otimizado (sem PETSC_USE_DEBUG)"
  fi
  grep -ao '\-march=[a-z0-9-]*' "$PCONF" 2>/dev/null | sort -u | sed 's/^/   /'
  [ -f "$PETSC_DIR/lib/libHYPRE"*.so ] 2>/dev/null && echo "   HYPRE presente"
  ls "$PETSC_DIR"/lib/libHYPRE*.so* >/dev/null 2>&1 && \
    echo "   HYPRE $(basename $(ls "$PETSC_DIR"/lib/libHYPRE*.so | head -1))"
else
  echo "   *** AUSENTE: nao achei $PCONF ***"
fi

echo "-- t8code"
if [ -f "$T8CODE/lib/libt8.so" ]; then
  echo "   $(basename $(readlink -f $T8CODE/lib/libt8.so))   prefixo $T8CODE"
else
  echo "   ausente (opcional: usado pelos exemplos 2D e example3d_complex)"
fi

echo "-- libfyaml"
pkg-config --modversion libfyaml 2>/dev/null | sed 's/^/   /' || echo "   *** AUSENTE ***"
n=$(ls /usr/lib/*/libfyaml.so.0 /usr/local/lib/libfyaml.so.0 2>/dev/null | wc -l)
[ "$n" -gt 1 ] && echo "   *** $n libfyaml instalados: qual resolve depende do ld.so ***"

echo "-- executaveis"
for b in "$HIGFLOW_RAIZ"/higflow/example3d_FrontTracking/ns-example-3d \
         "$HIGFLOW_RAIZ"/higflow/example3d_VOF_rising/ns-example-3d \
         "$HIGFLOW_RAIZ"/higflow/example3d_VOF_rising-ve/ns-example-3d; do
  if [ -x "$b" ]; then
    nf=$(ldd "$b" 2>/dev/null | grep -c 'not found')
    printf "   %-58s %s\n" "$(basename $(dirname $b))" \
      "$([ "$nf" -eq 0 ] && echo ok || echo "*** $nf biblioteca(s) nao resolvida(s) ***")"
  fi
done

echo "-- arreios de verificacao"
ls "$HIGFLOW_RAIZ"/higflow/tests/front-tracking/teste-* 2>/dev/null | wc -l | \
  awk '{print "   "$1" arreio(s) compilado(s)"}'

echo "-- montagens esperadas  [E4]"
for d in /scratch /work /dados; do
  printf "   %-10s %s\n" "$d" "$([ -d $d ] && echo existe || echo ausente)"
done
echo "=============================================================="
