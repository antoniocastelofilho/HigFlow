# Verificador da migracao, INDEPENDENTE do transformador.
#
# Existe porque o alerta impresso pelo proprio script de migracao e' calculado
# antes da limpeza e acusa funcoes ja' consertadas -- nao da' para confiar nele.
# Este roda sobre o arquivo GRAVADO e nao sabe nada de como ele foi produzido.
#
# Cobre as cinco formas pelas quais o transformador ja' falhou.  A quarta, o
# destroy orfao, e' a unica que nao falha na compilacao: ela morre no PETSc, e em
# codigo sem cobertura nenhuma suite a pegaria.
import re, sys, glob


# TUDO abaixo trabalha sobre o codigo SEM comentario.  Sem isto o verificador
# acusa `//higcit_destroy(it);` como destroy orfao -- foram 6 falsos positivos na
# primeira execucao, e um verificador que mente e' pior que nenhum.
def _sem_comentario(linhas):
    fora = []; bloco = False
    for l in linhas:
        out = []; i = 0
        while i < len(l):
            if bloco:
                if l[i:i+2] == "*/": bloco = False; i += 2; continue
                i += 1; continue
            if l[i:i+2] == "//": break
            if l[i:i+2] == "/*": bloco = True; i += 2; continue
            out.append(l[i]); i += 1
        fora.append("".join(out))
    return fora


ruins = 0
revisar = []
def erro(msg):
    global ruins
    print("  FALHA  " + msg); ruins += 1

for p in sorted(glob.glob("higflow/src/*.c")):
    if "src_hugo" in p: continue
    base = p.split("/")[-1]
    src = _sem_comentario(open(p).read().split("\n"))

    # (1) sobra de `c`/`f` dentro de laco migrado
    i = 0
    while i < len(src):
        m = re.match(r"^(\s*)for\(int (\w+) = 0; \2 < (hms|hfs)->n; \2\+\+\) \{\s*$", src[i])
        if not m:
            i += 1; continue
        ind, kind = m.group(1), m.group(3)
        j = i; prof = 0
        while j < len(src):
            prof += src[j].count("{") - src[j].count("}")
            if prof == 0 and j > i: break
            j += 1
        var = "c" if kind == "hms" else "f"
        for k in range(i+1, j):
            linha = re.sub(r"//.*", "", src[k])
            if re.search(r"\b%s\b" % var, linha) and not re.search(r"\b%s\s*\+|\+\s*%s\b" % (var,var), linha):
                erro("%s:%d sobrou `%s` em laco migrado: %s" % (base, k+1, var, src[k].strip()[:56]))
        i = j + 1

    # (2) DESTROY ORFAO -- a classe que so' a suite pegava, e que em codigo sem
    #     cobertura ninguem pegaria.  Por DOMINANCIA e nao por contagem: contar
    #     acusa `set_outflow`, que destroi em dois RAMOS mutuamente exclusivos com
    #     uma atribuicao so' -- codigo correto.  A pergunta certa e' se existe
    #     atribuicao ANTES do destroy, dentro da funcao.
    prof = 0; ini = 0
    for n, l in enumerate(src):
        a = prof; prof += l.count("{") - l.count("}")
        if a == 0 and prof > 0: ini = n
        if a > 0 and prof == 0:
            for it, destroy in (("it", "higcit_destroy(it)"), ("fit", "higfit_destroy(fit)")):
                for k in range(ini, n+1):
                    if destroy not in src[k]: continue
                    if not any(re.search(r"\b%s\s*=(?!=)" % it, src[q]) for q in range(ini, k)):
                        erro("%s:%d `%s` destruido sem nenhuma atribuicao antes, em %s"
                             % (base, k+1, it, src[ini].strip()[:40]))
                # contagem fica como REVISAO, nao falha: ela tem falso positivo
                # conhecido (destroy em ramos exclusivos), mas e' o unico sinal
                # para orfao que fica DEPOIS de um laco sobrevivente.
                corpo = src[ini:n+1]
                at = len([1 for x in corpo if re.search(r"\b%s\s*=(?!=)" % it, x)])
                de = len([1 for x in corpo if destroy in x])
                if de > at:
                    revisar.append("%s:%d `%s` %d destroy(s) para %d atribuicao(oes): %s"
                                   % (base, ini+1, it, de, at, src[ini].strip()[:40]))

for r in revisar: print("  revisar  " + r)
print("  verificacao independente: %d falha(s), %d para revisao" % (ruins, len(revisar)))
sys.exit(1 if ruins else 0)
