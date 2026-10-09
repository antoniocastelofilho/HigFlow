#!/bin/bash
# Acompanha o ORCAMENTO DE TENSAO SUPERFICIAL ARTIFICIAL do filtro, em producao.
#
# O teste do elipsoide PROJETOU o custo (0,021 de esfericidade em 75 aplicacoes
# sem limiar, contra 0,00000 com limiar).  Isto MEDE o custo na corrida de
# verdade, somando o dA/A que cada aplicacao imprime.  A relacao e' exata a
# volume conservado: Psi = (36 pi V^2)^(1/3)/A, logo dPsi/Psi = -dA/A.
#
# A projecao a partir das primeiras aplicacoes seria otimista: o numero de
# vertices acima do limiar CRESCE conforme a ondulacao se desenvolve.  Por isso
# isto acompanha em vez de extrapolar uma vez.
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$AQUI/registro"
L=$AQUI/registro/corrida.log
REG=$AQUI/registro/orcamento.log
python3 - "$L" "$REG" <<'PY'
import re, sys, time, os
L, REG = sys.argv[1], sys.argv[2]
vistos = 0
while True:
    try: txt = open(L, errors='replace').read()
    except OSError: txt = ""
    lin = re.findall(r'===> TSUR passo (\d+): BALANCE moveu (\d+), TSUR (\d+), '
                     r'chi=(-?\d+), volume ([-+0-9.e]+), AREA ([-+0-9.e]+)', txt)
    if len(lin) > vistos:
        dA = sum(float(x[5]) for x in lin); dV = sum(float(x[4]) for x in lin)
        mov = [x for x in lin if int(x[1]) + int(x[2]) > 0]
        chi_ok = all(x[3] == '2' for x in lin)
        psi = -dA * 0.96
        with open(REG, "a") as f:
            f.write("%s  aplic %3d (com movimento %3d)  soma dA/A %+.3e  dV/V %+.3e  "
                    "-> dPsi %+.2e  (%.1f%% de 0,0023)  chi %s\n"
                    % (time.strftime("%m-%d %H:%M"), len(lin), len(mov), dA, dV, psi,
                       100*abs(psi)/0.0023, "ok" if chi_ok else "VARIOU"))
        vistos = len(lin)
        if not chi_ok:
            with open(REG,"a") as f: f.write("  *** chi != 2 numa aplicacao -- topologia rompida ***\n")
            sys.exit(2)
        # O CRITERIO DE ABORTO MUDOU, porque o primeiro estava errado.
        #
        # Ele disparava quando a soma de -dA/A passava de 1/4 do sinal FT-VOF
        # (0,0023), tratando toda reducao de area como tensao artificial.  Ela
        # nao e': remover zigue-zague remove a area ESPURIA que o zigue-zague
        # acrescentou, e isso e' correcao, nao contaminacao.  As duas misturadas
        # num numero so' nao se separam.
        #
        # O discriminador medido contra o gemeo SEM filtro em N=50, no mesmo t:
        # a area do filtrado ficou 6,6e-04 menor e a esfericidade 4,8e-04 maior,
        # mas a ALTURA da bolha (zmax-zmin) ficou +2,7e-04 -- mais alta, nao mais
        # baixa.  Tensao artificial faria mais baixa e mais redonda.  Logo a area
        # sai da escala submalha e a forma resolvida nao se move.
        #
        # O que SERIA tensao artificial e' o filtro deixar de ser seletivo e
        # passar a mover o grosso da superficie.  E' isso que se vigia agora:
        # a fracao de vertices movidos numa aplicacao.  A soma de area continua
        # no registro como informacao, sem abortar por ela.
        NV = 2562                      # icosfera nsub=4
        pior = max((int(x[1])+int(x[2])) for x in lin)
        if pior > 0.25*NV:
            with open(REG,"a") as f:
                f.write("  *** O FILTRO DEIXOU DE SER SELETIVO: %d de %d vertices "
                        "numa aplicacao ***\n" % (pior, NV))
                f.write("  *** acima de 1/4 da superficie nao e' cirurgia, e' suavizacao\n"
                        "      global -- o limiar precisa subir\n")
            sys.exit(3)
    if not os.system("pgrep -f 'lanca\\.sh$' >/dev/null 2>&1") == 0:
        with open(REG,"a") as f: f.write("%s  o lancador saiu; encerrando o acompanhamento\n" % time.strftime("%m-%d %H:%M"))
        sys.exit(0)
    time.sleep(900)
PY
