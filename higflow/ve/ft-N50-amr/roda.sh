#!/bin/bash
# FRONT-TRACKING ADAPTATIVO, lambda=0.4 beta=0.2, COM O FILTRO.
#
# COMPARAVEL AO N=50 UNIFORME, e isso e' um requisito e nao um detalhe:
#
#   base 25x25x75 (h=0.04) + 1 NIVEL de refino  ->  mais fino h=0.02
#
# que e' exatamente o h da corrida uniforme 50x50x150.  Base 50 com 2 niveis
# daria h=0.01 e compararia resolucoes diferentes; base 13 com 2 niveis daria
# h=1/52, proximo mas nao igual.  Um nivel sobre base 25 e' o unico arranjo que
# iguala o mais fino.
#
# E TODO O RESTO e' igual por construcao: os tres yaml (par.contr, bc, init) sao
# COPIA da corrida uniforme, nao versoes paralelas.  Isso importa porque o
# exemplo AMR vinha configurado como NEWTONIANO (flowtype0: newtonian, De=1.0,
# beta=0.5) -- rodar assim compararia um caso newtoniano adaptativo com um
# viscoelastico uniforme, ou seja, nada.  Conferido campo por campo.
#
# RETOMAVEL.  Este roteiro roda UM TRECHO: se o diretorio ja' tem um ponto de
# salvamento valido, retoma dele; senao comeca do zero.  Chamado de novo, segue
# de onde parou.  Quem decide quando chamar e' o `lanca.sh`.
#
# O que faz um ponto ser VALIDO: o `.ok` existe e o passo dele e' o mesmo do
# init.yaml.  O `.ok` e' gravado por rename depois de todos os `._dp*`, logo
# divergencia entre os dois significa gravacao interrompida no meio -- conjunto
# metade novo e metade velho, que nao serve.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$AQUI/registro" 2>/dev/null
M=/home/castelo/HiGFlow-System/HigFlow
TOTAL=2001                       # passos da corrida inteira
cd "$M" && set -a && . ./varsrc && set +a
cd "$AQUI"

LOG=$AQUI/registro/corrida.log
# A impressao digital do binario no registro.  O diretorio da corrida tem uma
# COPIA do executavel, e ela ja' ficou atras de duas correcoes do ponto de
# salvamento -- a corrida teria rodado sem elas e o registro nao diria nada.
# Registrar faz a divergencia aparecer no mesmo lugar onde se olha o progresso.
echo "$(date +%H:%M) binario $(md5sum ns-example-3d | cut -c1-12) (fonte: $(md5sum "$M/higflow/example3d_FrontTracking-amr/ns-example-3d" 2>/dev/null | cut -c1-12))" >> $AQUI/registro/status.log
PASSO=0
if [ -f output/o.save.init.yaml ] && [ -f output/o.save.ok ]; then
    PI=$(awk '/step:/{print $2}' output/o.save.init.yaml)
    PO=$(awk '/^passo/{print $2}' output/o.save.ok)
    # um `.gravando` sobrando e' informacao certa: aquela gravacao nao
    # terminou, e o conjunto corrente nao serve mesmo que os dois numeros
    # batam (o `.ok` pode ser o da geracao anterior).
    if [ -f output/o.save.gravando ]; then
        echo "$(date +%H:%M) gravacao INTERROMPIDA (sinal .gravando presente)" >> $AQUI/registro/status.log
        PI=-1
    fi
    if [ "$PI" = "$PO" ] && [ "${PI:-0}" -gt 0 ]; then
        PASSO=$PI
    else
        echo "$(date +%H:%M) ponto de salvamento INCONSISTENTE (init=$PI ok=$PO)" >> $AQUI/registro/status.log
        # a copia de seguranca, se houver, e' o conjunto anterior COMPLETO
        if [ -f seguranca/o.save.ok ]; then
            PS=$(awk '/^passo/{print $2}' seguranca/o.save.ok)
            echo "$(date +%H:%M) voltando a' copia de seguranca, passo $PS" >> $AQUI/registro/status.log
            cp -f seguranca/o.save.* output/ 2>/dev/null
            rm -f output/o.save.gravando
            PASSO=$PS
        else
            echo "$(date +%H:%M) sem copia de seguranca -- RECOMECANDO do zero" >> $AQUI/registro/status.log
            rm -f output/o.save.gravando
            PASSO=0
        fi
    fi
fi

export FT3_RHO0=1000.0 FT3_RHO1=100.0 FT3_MU0=10.0 FT3_MU1=1.0 FT3_SIGMA=24.5
export FT3_R=0.25 FT3_CX=0.5 FT3_CY=0.5 FT3_CZ=0.5 FT3_NSUB=4 FT3_ADVECTA=1
export FT3_SALVA=1
# O FILTRO DE ONDULACOES, que e' a razao desta corrida existir.  Na metade SEM
# filtro a ondulacao vai de 0 a 21% das arestas com maximo 0,997 -- aresta quase
# perpendicular a' superficie local --, e toda ela na metade inferior da bolha.
#
# O LIMIAR NAO E' OPCIONAL.  Sem ele, 75 aplicacoes num elipsoide liso mudam a
# esfericidade em 0,021 -- NOVE VEZES a diferenca FT-VOF de 0,0023 que o
# trabalho mede.  Reduzir area a volume fixo e' o que a tensao superficial faz, e
# conservar volume nao protege disso.  Em 0,05 o elipsoide liso fica INTOCADO
# (zero vertices movidos por passe) e a ondulacao ainda cai 82%.
export FT3_TSUR_CADA=20
export FT3_TSUR_LIM=0.05
export FT3_TSUR_RLX=0.0
# ---- O REFINO ADAPTATIVO ----
# NIVEIS=1 sobre base 25 da' mais fino 50, igual ao uniforme (ver o cabecalho).
# CELMIN=5 e' o piso de celulas finas de cada lado da interface; medido, a banda
# que sai e' de 7 a 8 celulas, com folga sobre o exigido.
# CADA=25 remalha a cada 25 passos.  O laco dinamico foi verificado: banda
# cumprida, zero posicoes sem valor, suporte_nivel_trocado=0.
export FT3_AMR_NIVEIS=1
export FT3_AMR_CELMIN=5
# BANDA FISICA, ancorada para a serie poder convergir.  Com celmin FIXO cada
# resolucao refina largura fisica DIFERENTE -- a mais fina refina metade do
# espaco --, e foi isso que impediu o estudo de convergencia em 2D: nenhuma das
# cinco series ficou monotona, e a nao-monotonicidade foi de 12 a 40 vezes o piso
# de amostragem.  0,10 = 5 celulas finas NESTA resolucao, logo aqui nada muda; o
# ganho e' que uma corrida com outro h refina a MESMA largura.  O piso de 5 vence
# se a banda pedida for fina demais, e o exemplo avisa quando vence.
#
# A BANDA NAO RESOLVE A ESTEIRA: 0,10 em torno de uma bolha de R=0,25 cobre a
# casca de r=0,15 a 0,35, e a esteira se estende bem abaixo.  Alargar a banda
# ate' alcanca-la refinaria uma casca espessa inteira -- caro e no lugar errado.
# A esteira pede SENSOR (criterio-esteira-3d), nao banda.
export FT3_AMR_BANDA=0.10
export FT3_AMR_CADA=25
export FT3_AMR_NX=25 FT3_AMR_NY=25 FT3_AMR_NZ=75
export FT3_AMR_LX=1.0 FT3_AMR_LY=1.0 FT3_AMR_LZ=3.0
# O CAMINHO TEM DE SER O QUE O YAML LE'.  Escrever o criterio num arquivo que
# ninguem rele' faz a corrida rodar UNIFORME em silencio -- o log anuncia o
# refino, o .amr e' gerado, e o solver monta a malha base.  Era o que este
# roteiro fazia quando eu o derivei do uniforme: criterio.amr contra
# amrs/gota_d/dominio-0.amr do domain.yaml.
export FT3_AMR_CAMINHO=amrs/gota_d/dominio-0.amr
# CADENCIA DO CAMPO.  Em N=50 um quadro de campo sai com ~300 MB e leva ~20 min
# -- mais que um passo.  A animacao da bolha NAO precisa dele: a interface ja'
# sai com p, vel, S e T interpolados nos vertices.  O campo serve as' linhas de
# corrente, e para isso tres amostras bastam: 80 quadros de dtp=0.025 e' t=2.0.
export FT3_CAMPO_CADA=80
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

if [ "$PASSO" -gt 0 ]; then
    REST=$(( TOTAL - PASSO ))
    if [ "$REST" -le 0 ]; then
        echo "$(date +%H:%M) ja' completou $PASSO de $TOTAL -- nada a fazer" >> $AQUI/registro/status.log
        exit 0
    fi
    sed -i "s/^  numsteps: .*/  numsteps: $REST/" output/o.save.par.contr.yaml
    echo "$(date +%H:%M) RETOMANDO do passo $PASSO, faltam $REST" >> $AQUI/registro/status.log
    CARGA=output/o.save
else
    echo "$(date +%H:%M) comecando do zero, $TOTAL passos" >> $AQUI/registro/status.log
    CARGA=input/gota.load
fi

# nucleo fisico 1 (o 0 fica para o sistema); os pares SMT deste i9 sao (k,k+8),
# logo 1 e 9 sao o MESMO nucleo -- usar numeros contiguos nao da' nucleos
# distintos, e e' por isso que ha' so' um numero aqui.
taskset -c 1 stdbuf -oL -eL /usr/bin/time -f "PAREDE %e s" mpirun --bind-to none -use-hwthread-cpus -n 1 \
  ./ns-example-3d "$CARGA" output/o.save VTKS/p.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 \
  >> $LOG 2>&1
R=$?
echo "$(date +%H:%M) trecho terminou exit=$R, passo salvo=$(awk '/step:/{print $2}' output/o.save.init.yaml 2>/dev/null)" >> $AQUI/registro/status.log
exit $R
