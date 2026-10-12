#!/usr/bin/env python3
"""Conduz as medidas de custo da franja e sumariza.

POR QUE UM CONDUTOR EM PYTHON, e nao um laco no shell: a maior parte do que torna
uma medida confiavel esta' FORA do cronometro -- saber quantos nucleos a maquina
tem, se havia processo alheio consumindo CPU, e se o binario e' o da arvore
atual.  Sem isso, o numero e' sem procedencia, e neste projeto um build sem o que
se queria medir ja' produziu duas corridas "estaveis" que eram outra coisa.

A SAIDA E' UMA TABELA por nivel, com a fracao do tempo de solver -- que e' a unica
forma em que o custo da franja e' interpretavel.  Tempo absoluto de franja nao e'
resultado: 50 ms pode ser desprezivel ou dominante.
"""

import argparse
import os
import re
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(AQUI, "franja")


def estado_da_maquina():
    nucleos = os.cpu_count() or 1
    try:
        with open("/proc/loadavg") as f:
            carga = f.read().split()[:3]
    except OSError:
        carga = ["?"]
    return nucleos, " ".join(carga)


def alheios(limiar=20.0):
    """Processos consumindo CPU que nao sao nossos -- a medida nao vale com eles."""
    try:
        r = subprocess.run(["ps", "-eo", "pcpu,comm", "--sort=-pcpu"],
                           capture_output=True, text=True, timeout=10)
    except Exception:
        return []
    fora = []
    for linha in r.stdout.splitlines()[1:8]:
        p = linha.split(None, 1)
        if len(p) != 2:
            continue
        try:
            uso = float(p[0])
        except ValueError:
            continue
        nome = p[1].strip()
        if uso >= limiar and "franja" not in nome and "mpirun" not in nome:
            fora.append("%5.1f%%  %s" % (uso, nome))
    return fora


def roda(np, nivel, camadas, fracao, reps, cpu_set=None, sem_poisson=False):
    cmd = ["mpirun", "--oversubscribe", "-np", str(np)]
    if cpu_set:
        cmd += ["--cpu-set", cpu_set, "--bind-to", "core"]
    cmd += [BIN, "--nivel=%d" % nivel, "--camadas=%d" % camadas,
            "--fracao=%.2f" % fracao, "--reps=%d" % reps]
    if sem_poisson:
        cmd.append("--sem-poisson")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        print("  TIMEOUT em np=%d nivel=%d" % (np, nivel))
        return {}
    if r.returncode != 0:
        print("  NAO TERMINOU em np=%d nivel=%d; ultimas linhas:" % (np, nivel))
        print("\n".join("    " + l for l in (r.stdout + r.stderr).splitlines()[-8:]))
        return {}
    saida = {}
    for linha in r.stdout.splitlines():
        m = re.match(r"medida (\S+) (.*)", linha)
        if not m:
            if linha.startswith("AVISO"):
                print("  " + linha)
            continue
        etapa = m.group(1)
        campos = dict(kv.split("=", 1) for kv in m.group(2).split() if "=" in kv)
        chave = (etapa, int(campos.get("camadas", 0)))
        saida[chave] = campos
    return saida


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nps", default="1,2,3,4")
    p.add_argument("--niveis", default="5,6,7")
    p.add_argument("--camadas", type=int, default=4)
    p.add_argument("--fracao", type=float, default=0.0)
    p.add_argument("--reps", type=int, default=4)
    p.add_argument("--cpu-set", default=None,
                   help="lista de nucleos para o mpirun, p.ex. 4-15; use quando "
                        "houver outra corrida presa a nucleos especificos")
    a = p.parse_args()

    if not os.path.exists(BIN):
        print("binario ausente: %s  (construa com `make bench`)" % BIN)
        return 1

    nucleos, carga = estado_da_maquina()
    print("maquina: %d nucleos, carga %s" % (nucleos, carga))
    fora = alheios()
    if fora:
        print("ATENCAO -- processos alheios consumindo CPU durante a medida:")
        for x in fora:
            print("  " + x)
        print("  (o tempo do rank mais lento inclui a disputa com eles)")
    print("binario: %s  (%d bytes)" % (BIN, os.path.getsize(BIN)))
    print()

    nps = [int(x) for x in a.nps.split(",")]
    niveis = [int(x) for x in a.niveis.split(",")]

    for nivel in niveis:
        print("=" * 78)
        print("NIVEL %d  (fracao de refino aleatorio: %.2f)" % (nivel, a.fracao))
        print("=" * 78)
        cab = ("  np  camadas  celulas  franja     bytes/cel   franja(s)  "
               "fluxo2(s)  mls(s)  %fluxo  %mls")
        print(cab)
        for np in nps:
            d = roda(np, nivel, a.camadas, a.fracao, a.reps, a.cpu_set)
            if not d:
                continue
            t_f2 = float(d.get(("poisson-fluxo2", 0), {}).get("seg", "nan"))
            t_mls = float(d.get(("poisson-mlsgauss", 0), {}).get("seg", "nan"))
            for k in range(1, a.camadas + 1):
                c = d.get(("franja", k))
                if not c:
                    continue
                cel = int(c["globais"])
                by = int(c["bytes"])
                t = float(c["seg"])
                pf = 100.0 * t / t_f2 if t_f2 == t_f2 and t_f2 > 0 else float("nan")
                pm = 100.0 * t / t_mls if t_mls == t_mls and t_mls > 0 else float("nan")
                print("  %2d  %7d  %7d  %6s  %10.1f  %10.6f  %9.4f  %6.4f  %5.1f  %5.1f"
                      % (np, k, cel, c["franjamax"], by / max(1, cel), t,
                         t_f2, t_mls, pf, pm))
            rep = d.get(("reparticao", 3))
            if rep:
                print("      reparticao (3 camadas): %.6f s" % float(rep["seg"]))
        print()

    print("LEITURA: `camadas=1` tem de sair em tempo DESPREZIVEL -- `exige_camadas(1)`")
    print("nao faz rodada nenhuma.  Se sair com tempo, a medida nao esta' medindo as")
    print("rodadas de troca, e sim outra coisa.")
    print()
    print("`bytes/cel` e' um TETO: conta o que foi posto no canal, e o t8code")
    print("transmite menos -- so' as posicoes de franja que interessam a cada vizinho.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
