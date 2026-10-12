#!/usr/bin/env python3
"""Roda a matriz de escalabilidade e tabula, com a configuracao da maquina na saida.

A configuracao vai na saida porque a maquina tem vizinhos: medir sem declarar
quantos nucleos estao ocupados por outra coisa e' como medir sem dizer a
configuracao -- e processos alheios a 99% de CPU ja' fizeram eu medir 283 s no que
custa 5 s neste projeto.

Cada ponto e' o MENOR de `--repete` execucoes: com vizinhos na maquina, o minimo
e' menos sensivel a eles que a media.
"""

import argparse
import os
import re
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(AQUI, "escalabilidade")


def estado_da_maquina():
    """Nucleos, carga e processos alheios pesados -- para a medida ser interpretavel."""
    nucleos = os.cpu_count() or 1
    with open("/proc/loadavg") as f:
        carga = f.read().split()[0]
    alheios = []
    try:
        saida = subprocess.run(["ps", "-eo", "pcpu=,comm=,etime="], capture_output=True,
                               text=True, timeout=30).stdout
        for linha in saida.splitlines():
            partes = linha.split(None, 2)
            if len(partes) < 3:
                continue
            try:
                cpu = float(partes[0])
            except ValueError:
                continue
            if cpu > 50.0 and "escalabilidade" not in partes[1]:
                alheios.append("%s (%.0f%% CPU, %s)" % (partes[1], cpu, partes[2].strip()))
    except Exception:
        pass
    return nucleos, carga, alheios


def uma(np, dim, nivel, repete, cpu_set=None):
    """Roda um ponto; devolve {etapa: (globais, locmin, locmax, seg)} com o menor tempo.

    `cpu_set` restringe os ranks a um conjunto de CPUs.  Nao e' refinamento: nesta
    maquina ha' corridas alheias PRESAS a cpus especificas, e medir sem excluir
    aquelas cpus produziu ganho SUPERLINEAR de 18x -- numero que so' podia ser
    ruido.  A maquina tem 8 cores fisicos com hyperthreading, nao 16 cores.
    """
    melhor = {}
    for _ in range(repete):
        cmd = ["mpirun"]
        if cpu_set:
            cmd += ["--cpu-set", cpu_set, "--bind-to", "core"]
        cmd += ["-n", str(np), BIN, "--dim=%d" % dim, "--nivel=%d" % nivel]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        except subprocess.TimeoutExpired:
            print("  TIMEOUT em np=%d dim=%d nivel=%d" % (np, dim, nivel))
            return {}
        if "fim" not in r.stdout:
            print("  NAO TERMINOU em np=%d dim=%d nivel=%d; ultimas linhas:" % (np, dim, nivel))
            print("\n".join("    " + l for l in (r.stdout + r.stderr).splitlines()[-8:]))
            return {}
        for linha in r.stdout.splitlines():
            m = re.match(r"medida (\S+) dim=\d+ nivel=\d+ np=\d+ globais=(\d+) "
                         r"locmin=(\d+) locmax=(\d+) seg=([\d.]+)", linha)
            if m:
                et = m.group(1)
                v = (int(m.group(2)), int(m.group(3)), int(m.group(4)), float(m.group(5)))
                if et not in melhor or v[3] < melhor[et][3]:
                    melhor[et] = v
    return melhor


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dim", type=int, default=3, choices=(2, 3))
    p.add_argument("--niveis", default="5,6")
    p.add_argument("--nps", default="1,2,4,8")
    p.add_argument("--repete", type=int, default=2)
    p.add_argument("--cpu-set", default=None,
                   help="cpus a usar, p.ex. 0,3,4,5 -- excluindo as de corridas alheias")
    a = p.parse_args()

    if not os.path.exists(BIN):
        print("binario ausente: %s  (construa com `make bench`)" % BIN)
        return 2

    nucleos, carga, alheios = estado_da_maquina()
    print("maquina: %d nucleos, carga %s" % (nucleos, carga))
    if alheios:
        print("ATENCAO -- processos alheios consumindo CPU durante a medida:")
        for x in alheios:
            print("  " + x)
        print("  os tempos abaixo estao contaminados por eles; e' por isso que a")
        print("  medida reporta o MENOR de %d execucoes." % a.repete)
    else:
        print("nenhum processo alheio pesado")
    if a.cpu_set:
        print("medindo RESTRITO as cpus %s (as alheias ficam de fora)" % a.cpu_set)
    print()

    niveis = [int(x) for x in a.niveis.split(",")]
    nps = [int(x) for x in a.nps.split(",")]

    for nivel in niveis:
        base = {}
        print("=== dim=%d nivel=%d ===" % (a.dim, nivel))
        print("%-18s %4s %10s %16s %10s %9s %9s %7s" %
              ("etapa", "np", "globais", "locais min/max", "seg", "ganho", "efic", "deseq"))
        linhas = []
        for np in nps:
            r = uma(np, a.dim, nivel, a.repete, a.cpu_set)
            if not r:
                continue
            for et in ("construcao", "percurso", "refino", "reparticao", "percurso-adaptado"):
                if et not in r:
                    continue
                g, lmin, lmax, seg = r[et]
                if np == 1:
                    base[et] = seg
                ganho = (base[et] / seg) if et in base and seg > 0 else float("nan")
                efic = ganho / np if ganho == ganho else float("nan")
                linhas.append((et, np, g, lmin, lmax, seg, ganho, efic))
        # agrupa por etapa, para o ganho ser lido na vertical
        for et in ("construcao", "percurso", "refino", "reparticao", "percurso-adaptado"):
            for (e, np, g, lmin, lmax, seg, ganho, efic) in linhas:
                if e != et:
                    continue
                deseq = (lmax / lmin) if lmin > 0 else float("nan")
                print("%-18s %4d %10d %7d /%7d %10.4f %8.2fx %8.0f%% %6.2fx" %
                      (e, np, g, lmin, lmax, seg, ganho, 100 * efic, deseq))
            print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
