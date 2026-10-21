#!/usr/bin/env python3
"""Gera famílias de malha do example3d_complex para teste de escalonamento.

DOIS BOTÕES INDEPENDENTES, e confundi-los é o erro que este programa existe
para evitar:

  --m      CÉLULAS POR 0,1 DE COMPRIMENTO.  Muda o número de células e nada
           mais.  m = 10 reproduz a malha original (h = 0,01).  Todas as
           coordenadas dos 33 blocos são múltiplos de 0,1, então qualquer m
           inteiro dá contagem inteira em todo bloco.

  --blocos QUANTOS BLOCOS.  Muda a granularidade da partição e NÃO muda a
           malha: um bloco cortado em dois tem exatamente as mesmas células.

Por que os dois são necessários, e por que o segundo é o que estava faltando:
a higtree particiona ÁRVORES INTEIRAS (`lb_add_input_tree` recebe um
`hig_cell *`), então o maior bloco é um piso para o rank mais carregado.  No
caso original o maior bloco tem 24 000 das 162 000 células -- 14,8% --, e isso
limita o ganho a 6,75x por mais núcleos que se dê.  O teto é atingido em np = 8
e nada acima disso muda.  Refinar sozinho NÃO levanta esse teto: multiplica
todos os blocos igualmente e o maior continua com 14,8%.

AS CONDIÇÕES DE CONTORNO NÃO SÃO CORTADAS, e não precisam ser: a busca de BC em
`domain.c` é geométrica -- varre `d->dirichlet_bcs` e procura células numa caixa
em torno do ponto --, sem associação bloco <-> contorno.  Um retalho de BC pode
cobrir dois blocos.  Por isso os 169 retalhos mantêm número, ordem e id, e o
`switch (id)` compilado em ns-complex-3d.c continua válido sem reescrita.
Refinar, esse sim, mexe nos retalhos: a malha do contorno tem de acompanhar a
do domínio.

Uso:
  gera-malha-complexa.py --origem <dir> --destino <dir> [--m 10] [--blocos N]
  gera-malha-complexa.py --origem <dir> --tabela          (só relata, não grava)
"""

import argparse
import os
import re
import sys

EPSMACH = 1.0e-12            # a mesma de higtree/src/utils.h
UNIDADE = 0.1                # todas as coordenadas do caso são múltiplos disto


# --------------------------------------------------------------- leitura

def le_amr(caminho):
    """Lê um .amr uniforme (numlevels = 1, numpatches = 1).

    O formato, de higio_read_amr_info():  l[0] h[0] l[1] h[1] l[2] h[2]
    numlevels / delta[3] numpatches / initcell[3] patchsize[3].  A ordem é
    lo e hi INTERCALADOS por dimensão, não os três lo seguidos dos três hi --
    ler na ordem errada dá uma caixa que parece plausível.
    """
    v = [float(x) for x in open(caminho).read().split()]
    lo = [v[0], v[2], v[4]]
    hi = [v[1], v[3], v[5]]
    niveis = int(v[6])
    if niveis != 1:
        raise SystemExit("%s: numlevels = %d; este gerador só trata malha "
                         "uniforme" % (caminho, niveis))
    delta = v[7:10]
    retalhos = int(v[10])
    if retalhos != 1:
        raise SystemExit("%s: numpatches = %d; só trato 1" % (caminho, retalhos))
    inicial = [int(x) for x in v[11:14]]
    tam = [int(x) for x in v[14:17]]
    return lo, hi, delta, inicial, tam


def paths_do_yaml(caminho):
    """Os `path:` de um yaml de domínio ou de contorno, na ordem do arquivo."""
    return re.findall(r'^\s*path:\s*(\S+)', open(caminho).read(), re.M)


# --------------------------------------------------------------- escrita

def escreve_amr(caminho, lo, hi, tam):
    """Grava um bloco uniforme, e CONFERE a invariante que o leitor exige.

    higio_read_from_amr_info() tem um assert FLT_EQ(delta, (hi-lo)/patchsize)
    com EPSMACH = 1e-12 absoluto.  Se ele falhar, o executável aborta sem dizer
    qual arquivo -- então a conferência é aqui, com o nome do arquivo na mão.
    """
    linhas = []
    linhas.append(" ".join("%.12f %.12f" % (lo[k], hi[k]) for k in range(3)))
    linhas.append("1")
    delta = []
    for k in range(3):
        d = 0.0 if tam[k] == 0 or hi[k] == lo[k] else (hi[k] - lo[k]) / tam[k]
        delta.append(d)
    linhas.append("%.12f %.12f %.12f 1" % tuple(delta))
    linhas.append("1 1 1 %d %d %d" % tuple(tam))
    texto = "\n".join(linhas) + "\n"

    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with open(caminho, "w") as f:
        f.write(texto)

    # Reler do disco: o que o solver vê é o texto, não as variáveis daqui.
    rlo, rhi, rdelta, _, rtam = le_amr(caminho)
    for k in range(3):
        if rtam[k] != tam[k]:
            raise SystemExit("%s: patchsize voltou errado" % caminho)
        esperado = (rhi[k] - rlo[k]) / rtam[k]
        if abs(rdelta[k] - esperado) >= EPSMACH:
            raise SystemExit(
                "%s: dim %d, delta %.17g contra (hi-lo)/n %.17g, diferença "
                "%.3g >= EPSMACH -- o assert de higio_read_from_amr_info vai "
                "abortar" % (caminho, k, rdelta[k], esperado,
                             abs(rdelta[k] - esperado)))


def escreve_yaml_dominio(caminho, caminhos):
    with open(caminho, "w") as f:
        f.write("domain:\n")
        f.write("   number_domains: %d                      "
                "#*** The Total Number of the Domains ***\n" % len(caminhos))
        for i, p in enumerate(caminhos):
            f.write("   domain%d:\n" % i)
            f.write("      id: %d                              "
                    "#*** The Domain Identifier           ***\n" % i)
            f.write("      path: %s       "
                    "#*** The Path to the Domain          ***\n" % p)


def escreve_lista_dominio(caminho, caminhos):
    with open(caminho, "w") as f:
        f.write("%d\n" % len(caminhos))
        for p in caminhos:
            f.write("%s\n" % p)


def reescreve_paths(origem, destino, novos):
    """Copia um yaml trocando só os `path:`, linha a linha.

    Textual de propósito.  O yaml de contorno carrega a física de cada id
    (type, value_type) e a ordem das entradas É o id que o `switch` do exemplo
    usa; reemitir por um serializador arriscaria reordenar.  Aqui o resto do
    arquivo sai byte a byte igual.
    """
    linhas = open(origem).readlines()
    k = 0
    for i, linha in enumerate(linhas):
        m = re.match(r'^(\s*path:\s*)(\S+)(.*)$', linha)
        if m:
            linhas[i] = "%s%s%s\n" % (m.group(1), novos[k], m.group(3).rstrip())
            k += 1
    if k != len(novos):
        raise SystemExit("%s: %d `path:` no arquivo contra %d esperados"
                         % (origem, k, len(novos)))
    with open(destino, "w") as f:
        f.writelines(linhas)


# --------------------------------------------------------------- partição

def lpt(pesos, np_):
    """Maior rank sob o escalonamento guloso LPT -- otimista de propósito.

    Zoltan não é LPT; este número é um LIMITE SUPERIOR do ganho possível, não
    uma previsão do que a corrida vai dar.  Serve para separar "a malha não
    permite" de "o código não aproveita".
    """
    caixas = [0] * np_
    for p in sorted(pesos, reverse=True):
        i = min(range(np_), key=lambda k: caixas[k])
        caixas[i] += p
    return max(caixas)


def relata(pesos, nps=(1, 2, 4, 8, 16, 32, 64, 128, 256)):
    tot = sum(pesos)
    maior = max(pesos)
    out = []
    out.append("   blocos %d, celulas %d" % (len(pesos), tot))
    out.append("   maior bloco %d (%.1f%% da malha) -> teto absoluto %.1fx"
               % (maior, 100.0 * maior / tot, tot / maior))
    out.append("   %6s %12s %12s %10s" % ("np", "maior rank", "ideal", "teto"))
    for n in nps:
        m = lpt(pesos, n)
        out.append("   %6d %12d %12d %9.1fx" % (n, m, tot // n, tot / m))
    return "\n".join(out)


def corta(blocos, alvo):
    """Divide o maior bloco ao meio até haver `alvo` blocos.

    Em índices INTEIROS de célula, então o corte cai sempre numa face de
    célula e a união é preservada exatamente -- não há arredondamento para
    conferir depois.  Cortar sempre o maior é o que minimiza o maior bloco.
    """
    blocos = [list(b) for b in blocos]
    while len(blocos) < alvo:
        i = max(range(len(blocos)), key=lambda k: celulas(blocos[k]))
        b = blocos[i]
        eixo = max(range(3), key=lambda k: b[1][k] - b[0][k])
        if b[1][eixo] - b[0][eixo] < 2:
            # Nada mais a cortar: todo bloco restante tem 1 célula no eixo
            # maior.  Acontece se o alvo for maior que o número de células.
            if all(bb[1][e] - bb[0][e] < 2 for bb in blocos for e in range(3)):
                break
            # Este não dá, mas outro pode: tira-o de cena marcando-o.
            cands = [k for k in range(len(blocos))
                     if any(blocos[k][1][e] - blocos[k][0][e] >= 2
                            for e in range(3))]
            if not cands:
                break
            i = max(cands, key=lambda k: celulas(blocos[k]))
            b = blocos[i]
            eixo = max(range(3), key=lambda k: b[1][k] - b[0][k])
        meio = (b[0][eixo] + b[1][eixo]) // 2
        a_lo, a_hi = list(b[0]), list(b[1])
        c_lo, c_hi = list(b[0]), list(b[1])
        a_hi[eixo] = meio
        c_lo[eixo] = meio
        blocos[i] = [a_lo, a_hi]
        blocos.append([c_lo, c_hi])
    return blocos


def celulas(b):
    return (b[1][0] - b[0][0]) * (b[1][1] - b[0][1]) * (b[1][2] - b[0][2])


# --------------------------------------------------------------- principal

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--origem", required=True,
                    help="diretório do example3d_complex")
    ap.add_argument("--destino", help="onde gravar a malha gerada")
    ap.add_argument("--prefixo", default="example-3d.load",
                    help="prefixo dos arquivos de entrada (nameload)")
    ap.add_argument("--m", type=int, default=10,
                    help="células por 0,1 de comprimento (10 = malha original)")
    ap.add_argument("--blocos", type=int, default=0,
                    help="alvo de blocos; 0 mantém os originais")
    ap.add_argument("--para-np", type=int, default=0,
                    help="corta até o maior bloco caber em células/(folga*np); "
                         "é o alvo em intenção, --blocos é o alvo em número")
    ap.add_argument("--folga", type=float, default=2.0,
                    help="quantas vezes menor que a fatia ideal o maior bloco "
                         "deve ficar (padrão 2); 1 só garante o teto absoluto, "
                         "e o empacotamento precisa de margem acima disso")
    ap.add_argument("--raiz-relativa", default="",
                    help="prefixo a pôr nos `path:` dos yamls, relativo ao "
                         "diretório de onde o executável roda")
    ap.add_argument("--tabela", action="store_true",
                    help="só relata a granularidade da origem e sai")
    a = ap.parse_args()

    org = a.origem
    dom_yaml = os.path.join(org, "input", a.prefixo + ".domain.yaml")
    bc_yaml = os.path.join(org, "input", a.prefixo + ".bc.yaml")
    dom_paths = paths_do_yaml(dom_yaml)
    bc_paths = paths_do_yaml(bc_yaml)

    # --- lê os blocos de domínio e passa para índices inteiros de célula
    orig = []
    for p in dom_paths:
        lo, hi, delta, ini, tam = le_amr(os.path.join(org, p))
        if ini != [1, 1, 1]:
            raise SystemExit("%s: initcell %s inesperado" % (p, ini))
        for k in range(3):
            n_unid = (hi[k] - lo[k]) / UNIDADE
            if abs(n_unid - round(n_unid)) > 1e-9:
                raise SystemExit(
                    "%s: a extensão %g na dim %d não é múltiplo de %g; o botão "
                    "--m supõe que é" % (p, hi[k] - lo[k], k, UNIDADE))
        orig.append((lo, hi, tam))

    h = UNIDADE / a.m
    # índice global de célula: i = round(x / h)
    def idx(x):
        i = x / h
        if abs(i - round(i)) > 1e-6:
            raise SystemExit("coordenada %g não cai na grade de h = %g" % (x, h))
        return int(round(i))

    blocos = []
    for lo, hi, tam in orig:
        blocos.append([[idx(lo[k]) for k in range(3)],
                       [idx(hi[k]) for k in range(3)]])

    total_orig = sum(celulas(b) for b in blocos)

    if a.tabela:
        print("ORIGEM  m = %d" % a.m)
        print(relata([celulas(b) for b in blocos]))
        return

    if not a.destino:
        raise SystemExit("--destino é obrigatório sem --tabela")

    # --- corta
    if a.para_np:
        teto = total_orig / (a.folga * a.para_np)
        alvo = len(blocos)
        while max(celulas(b) for b in corta(blocos, alvo)) > teto:
            alvo *= 2
            if alvo > total_orig:
                raise SystemExit("não há como levar o maior bloco abaixo de "
                                 "%.0f células" % teto)
        blocos = corta(blocos, alvo)
    elif a.blocos:
        blocos = corta(blocos, a.blocos)

    # --- as três conferências que valem a pena
    #  1. nenhuma célula criada nem perdida
    total = sum(celulas(b) for b in blocos)
    if total != total_orig:
        raise SystemExit("o corte mudou a contagem: %d -> %d" % (total_orig, total))
    #  2. nenhum par de blocos se sobrepõe
    for i in range(len(blocos)):
        for j in range(i + 1, len(blocos)):
            if all(blocos[i][0][k] < blocos[j][1][k] and
                   blocos[j][0][k] < blocos[i][1][k] for k in range(3)):
                raise SystemExit("blocos %d e %d se sobrepõem" % (i, j))
    #  3. nenhum bloco vazio
    for i, b in enumerate(blocos):
        if celulas(b) <= 0:
            raise SystemExit("bloco %d ficou vazio" % i)
    # (1) + (2) sobre um conjunto que já cobria o domínio implica que a união
    # continua a mesma: mesma contagem, sem sobreposição, dentro das caixas
    # originais por construção do corte.

    # --- grava
    dst = a.destino
    raiz = a.raiz_relativa
    novos_dom = []
    for i, b in enumerate(blocos):
        rel = "amrs/d/malha-d-%d.amr" % i
        lo = [b[0][k] * h for k in range(3)]
        hi = [b[1][k] * h for k in range(3)]
        tam = [b[1][k] - b[0][k] for k in range(3)]
        escreve_amr(os.path.join(dst, rel), lo, hi, tam)
        novos_dom.append(os.path.join(raiz, rel) if raiz else rel)

    novos_bc = []
    for i, p in enumerate(bc_paths):
        lo, hi, delta, ini, tam = le_amr(os.path.join(org, p))
        novo_tam = []
        for k in range(3):
            if hi[k] == lo[k]:
                novo_tam.append(1)          # a dimensão degenerada do retalho
            else:
                n_unid = (hi[k] - lo[k]) / UNIDADE
                novo_tam.append(int(round(n_unid * a.m)))
        rel = "amrs/bc/malha-bc-%d.amr" % i
        escreve_amr(os.path.join(dst, rel), lo, hi, novo_tam)
        novos_bc.append(os.path.join(raiz, rel) if raiz else rel)

    os.makedirs(os.path.join(dst, "input"), exist_ok=True)
    escreve_yaml_dominio(os.path.join(dst, "input", a.prefixo + ".domain.yaml"),
                         novos_dom)
    escreve_lista_dominio(os.path.join(dst, "input", a.prefixo + ".domain"),
                          novos_dom)
    reescreve_paths(bc_yaml, os.path.join(dst, "input", a.prefixo + ".bc.yaml"),
                    novos_bc)
    # o resto da entrada vai inalterado
    for sufixo in (".par.contr.yaml", ".init.yaml"):
        o = os.path.join(org, "input", a.prefixo + sufixo)
        if os.path.exists(o):
            with open(o) as f:
                conteudo = f.read()
            with open(os.path.join(dst, "input", a.prefixo + sufixo), "w") as f:
                f.write(conteudo)

    pesos = [celulas(b) for b in blocos]
    print("GERADA em %s   m = %d   h = %g" % (dst, a.m, h))
    print(relata(pesos))
    print("   retalhos de contorno: %d (número, ordem e id inalterados)"
          % len(novos_bc))


if __name__ == "__main__":
    main()
