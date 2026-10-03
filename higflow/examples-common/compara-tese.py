#!/usr/bin/env python3
"""Compara a matriz viscoelastica medida com os valores que a tese publica.

POR QUE ESTE SCRIPT EXISTE SEPARADO.  A tese publica CURVAS (Figuras 42 e 43),
nao tabelas; os unicos numeros dela sao os que o texto cita, com "aproximadamente".
Entao a comparacao tem duas naturezas, e misturar as duas seria o erro:

  QUANTITATIVA, e so' onde a tese da' numero: z_c em t=4, w maximo e w
  estacionario, para (lambda=0.1, beta=0.2), (lambda=0.4, beta=0.2) e o
  newtoniano de referencia.  Tolerancia frouxa de proposito -- o alvo tem dois
  digitos e vem precedido de "aproximadamente".

  QUALITATIVA, para as afirmacoes de TENDENCIA: mais polimero sobe mais alto e
  mais rapido; o viscoelastico sobe mais que o newtoniano; e a queda abrupta da
  velocidade depois do maximo em lambda=0.4, beta=0.2.  Essas sao ordenacoes e
  formas de curva, e sao verificaveis sem digitalizar figura nenhuma.

Uso:  compara-tese.py <dir da matriz ve>  [--dt-quadro 0.025]
"""
import sys, os, glob, math, subprocess

def metricas(d, dtq=0.025):
    """Chama o medidor de superficie e devolve a serie (t, z_c)."""
    aqui = os.path.dirname(os.path.abspath(__file__))
    med = os.path.join(aqui, 'metricas-superficie.py')
    csv = '/tmp/_cmp_%s.csv' % os.path.basename(d.rstrip('/'))
    r = subprocess.run([sys.executable, med, os.path.join(d,'VTKS'),
                        '--dt-quadro', str(dtq), '--csv', csv],
                       capture_output=True, text=True)
    if not os.path.exists(csv): return None
    import csv as _csv
    with open(csv) as f: rows = list(_csv.DictReader(f))
    s = [(float(x['t']), float(x['z_c'])) for x in rows if x['t']]
    return s or None

def velocidades(s):
    """w maximo, w estacionario (media do ultimo decimo), e se ha' queda
    abrupta depois do maximo."""
    if len(s) < 5: return None
    w = []
    for q in range(1, len(s)-1):
        w.append(((s[q+1][1]-s[q-1][1])/(s[q+1][0]-s[q-1][0]), s[q][0]))
    wmax, tmax = max(w)
    cauda = [x for x,_ in w[-max(3,len(w)//10):]]
    west = sum(cauda)/len(cauda)
    # queda abrupta: o minimo DEPOIS do maximo, antes da cauda
    dep = [x for x,t in w if t > tmax]
    queda = (wmax - min(dep))/wmax if dep else 0.0
    # e recupera depois?
    recupera = (dep and west > min(dep)*1.02)
    return wmax, tmax, west, queda, bool(recupera)

def main():
    if len(sys.argv) < 2:
        print(__doc__); return 1
    raiz = sys.argv[1]
    dtq = 0.025
    for q,a in enumerate(sys.argv):
        if a == '--dt-quadro' and q+1 < len(sys.argv): dtq = float(sys.argv[q+1])

    casos = sorted(d for d in glob.glob(os.path.join(raiz,'*'))
                   if os.path.isdir(d) and os.path.isdir(os.path.join(d,'VTKS')))
    dados = {}
    for d in casos:
        nome = os.path.basename(d)
        s = metricas(d, dtq)
        if s is None:
            print("  %-20s (sem serie)" % nome); continue
        v = velocidades(s)
        dados[nome] = (s, v)

    print("\n=== medido ===\n")
    print("  %-20s %8s %8s %8s %8s %10s %9s" % (
          "caso","t final","z_c","w max","em t=","w estac.","queda %"))
    for nome in sorted(dados):
        s,v = dados[nome]
        if v is None: continue
        print("  %-20s %8.2f %8.4f %8.4f %8.2f %10.4f %8.1f%%" % (
              nome, s[-1][0], s[-1][1], v[0], v[1], v[2], 100*v[3]))

    # ---- QUANTITATIVO: so' onde a tese da' numero --------------------
    alvos = {('0.1','0.2'): (1.9, 0.38, 0.37),
             ('0.4','0.2'): (2.0, 0.42, 0.40),
             ('newt','newt'): (1.8, 0.36, 0.35)}
    print("\n=== contra os numeros da tese (secao 8.2.5) ===\n")
    print("  %-20s %-22s %-22s %-22s" % ("caso","z_c (tese)","w max (tese)","w estac. (tese)"))
    for nome in sorted(dados):
        s,v = dados[nome]
        if v is None: continue
        if 'newt' in nome: chave=('newt','newt')
        else:
            p=nome.split('-')
            lam = '0.'+p[1][2:]; bet = '0.'+p[2][2:]
            chave=(lam,bet)
        if chave not in alvos: continue
        az,aw,ae = alvos[chave]
        f=lambda m,a: "%.4f vs %.2f  (%+.1f%%)"%(m,a,100*(m-a)/a)
        print("  %-20s %-22s %-22s %-22s"%(nome,f(s[-1][1],az),f(v[0],aw),f(v[2],ae)))

    # ---- QUALITATIVO: ordenacoes e forma de curva ---------------------
    print("\n=== afirmacoes de tendencia da tese ===\n")
    for met in ('ft','vof'):
        for lam in ('01','04'):
            tri=[]
            for bet in ('08','05','02'):
                n='%s-l%s-b%s'%(met,lam,bet)
                if n in dados and dados[n][1]:
                    tri.append((bet, dados[n][0][-1][1], dados[n][1][0], dados[n][1][2]))
            if len(tri)<3: continue
            # beta 0.8 -> 0.5 -> 0.2 e' MAIS polimero; z_c e w devem CRESCER
            zc_ok = tri[0][1] < tri[1][1] < tri[2][1]
            wm_ok = tri[0][2] < tri[1][2] < tri[2][2]
            we_ok = tri[0][3] < tri[1][3] < tri[2][3]
            print("  %s lambda=0.%s  mais polimero sobe mais alto: %-9s  (%.4f, %.4f, %.4f)"%(
                  met.upper(),lam[1],"CONFIRMA" if zc_ok else "REFUTA",*[t[1] for t in tri]))
            print("  %s lambda=0.%s  mais polimero mais rapido:   %-9s  w max (%.4f, %.4f, %.4f)"%(
                  met.upper(),lam[1],"CONFIRMA" if wm_ok else "REFUTA",*[t[2] for t in tri]))
            print("  %s lambda=0.%s  idem no estacionario:        %-9s  (%.4f, %.4f, %.4f)"%(
                  met.upper(),lam[1],"CONFIRMA" if we_ok else "REFUTA",*[t[3] for t in tri]))
            print()
    # viscoelastico contra newtoniano
    for met in ('ft','vof'):
        nn='%s-newt'%met
        if nn in dados and dados[nn][1]:
            zn=dados[nn][0][-1][1]
            for lam in ('01','04'):
                n='%s-l%s-b02'%(met,lam)
                if n in dados and dados[n][1]:
                    zv=dados[n][0][-1][1]
                    print("  %s lambda=0.%s b=0.2 sobe mais que o newtoniano: %-9s (%.4f vs %.4f)"%(
                          met.upper(),lam[1],"CONFIRMA" if zv>zn else "REFUTA",zv,zn))
    # a queda abrupta em lambda=0.4, beta=0.2
    print()
    for met in ('ft','vof'):
        n='%s-l04-b02'%met
        if n in dados and dados[n][1]:
            v=dados[n][1]
            print("  %s lambda=0.4 b=0.2 queda abrupta apos o maximo: %.1f%% e %s"%(
                  met.upper(),100*v[3],"RECUPERA" if v[4] else "NAO recupera"))
    print()
    return 0

if __name__ == '__main__':
    sys.exit(main())
