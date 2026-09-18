#!/usr/bin/env python3
"""Converte um exemplo da porta higflow_set_external_functions para HigFlowProblem.

Move as oito funcoes livres para metodos de uma classe, sem tocar nos corpos, e
troca a chamada de registro.  Verifica ao fim que os oito corpos sao identicos
aos de antes -- e' essa checagem que da' confianca no gesto repetido doze vezes.
"""
import re, sys, subprocess

NOMES = [("get_pressure","pressure"), ("get_velocity","velocity"),
         ("get_source_term","source_term"), ("get_facet_source_term","facet_source_term"),
         ("get_boundary_pressure","boundary_pressure"),
         ("get_boundary_velocity","boundary_velocity"),
         ("get_boundary_source_term","boundary_source_term"),
         ("get_boundary_facet_source_term","boundary_facet_source_term")]

def corpo_em(s, i):
    d=0; j=i
    while j < len(s):
        if s[j]=='{': d+=1
        elif s[j]=='}':
            d-=1
            if d==0: return s[i:j+1], j
        j+=1
    raise SystemExit("chave nao fecha")

def norm(x):
    return re.sub(r'\s+',' ',x).strip() if x else None

def corpo_por_padrao(s, pat):
    m = re.search(pat, s, re.M)
    if not m: return None
    c,_ = corpo_em(s, m.end()-1)
    return norm(c)

def converte(caminho, nome_classe, caminho_registro=None):
    """caminho: onde estao as oito funcoes.  caminho_registro: onde esta a chamada
    higflow_set_external_functions, se for outro arquivo (alguns exemplos separam
    as funcoes do usuario num .c incluido pelo principal)."""
    original = open(caminho).read()
    s = original
    ext = {}
    for velho, novo in NOMES:
        m = re.search(r'^real '+velho+r'\(([^)]*)\)\s*\{', s, re.M)
        if not m: raise SystemExit(f"nao achei {velho} em {caminho}")
        corpo, fim = corpo_em(s, m.end()-1)
        ini = m.start()
        k = s.rfind('\n', 0, ini-1)
        com = ''
        if k >= 0 and s[k+1:ini].lstrip().startswith('//'):
            com = s[k+1:ini].strip(); ini = k+1
        ext[velho] = (novo, m.group(1), corpo, com)
        s = s[:ini] + '\x00M_'+velho+'\x00' + s[fim+1:]

    met = []
    for velho, novo in NOMES:
        nv, args, corpo, com = ext[velho]
        if com: met.append("    " + com)
        linhas = corpo.split('\n')
        met.append(f"    real {nv}({args}) {linhas[0]}")      # a chave fica
        met += [("    "+l) if l.strip() else "" for l in linhas[1:]]

    classe = ("// ---------------------------------------------------------------------------\n"
              "// O problema deste exemplo, como um tipo em vez de oito funcoes soltas.\n"
              "// Os corpos sao os mesmos; so' mudaram de lugar e perderam o prefixo get_.\n"
              "// ---------------------------------------------------------------------------\n"
              f"class {nome_classe} : public HigFlowProblem {{\n"
              "public:\n" + '\n'.join(met) + "\n};\n\n"
              f"static {nome_classe} problema;\n")
    s = s.replace('\x00M_get_pressure\x00', classe, 1)
    for velho,_ in NOMES[1:]:
        s = s.replace('\x00M_'+velho+'\x00', '', 1)
    s = re.sub(r'\n{3,}', '\n\n', s)

    alvo = caminho_registro or caminho
    treg = s if alvo == caminho else open(alvo).read()
    r = re.search(r'\n[ \t]*higflow_set_external_functions\(ns,.*?\);', treg, re.S)
    if not r: raise SystemExit(f"nao achei a chamada de registro em {alvo}")
    treg = (treg[:r.start()] +
            "\n    // Registro por objeto: a interface substitui os oito ponteiros.\n"
            "    higflow_set_problem(ns, &problema);" + treg[r.end():])
    if alvo == caminho:
        s = treg
    else:
        open(alvo,'w').write(treg)
    open(caminho,'w').write(s)

    # verificacao: os oito corpos tem que ser identicos aos de antes
    ok = 0
    for velho, novo in NOMES:
        a = corpo_por_padrao(original, r'^real '+velho+r'\([^)]*\)\s*\{')
        b = corpo_por_padrao(s,        r'^\s*real '+novo+r'\([^)]*\)\s*\{')
        if a == b: ok += 1
        else: print(f"   DIVERGE: {velho} -> {novo}")
    print(f"  {caminho.split('/')[-2]}: {ok} de 8 corpos identicos")
    return ok == 8

if __name__ == '__main__':
    reg = sys.argv[3] if len(sys.argv) > 3 else None
    sys.exit(0 if converte(sys.argv[1], sys.argv[2], reg) else 1)
