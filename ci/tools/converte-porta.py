import re, sys, json
def corpo_em(s,i):
    d=0;j=i
    while j<len(s):
        if s[j]=='{': d+=1
        elif s[j]=='}':
            d-=1
            if d==0: return s[i:j+1], j
        j+=1
    raise SystemExit("chave nao fecha")
def norm(x): return re.sub(r'\s+',' ',x).strip() if x else None
def achar(s,pat):
    m=re.search(pat,s,re.M)
    if not m: return None
    c,_=corpo_em(s,m.end()-1); return norm(c)

def converte(caminho, base, pares, rx_chamada, nova_chamada, rotulo):
    orig=open(caminho).read(); s=orig
    ext={}
    for velho,novo in pares:
        m=re.search(r'^real '+velho+r'\(([^)]*)\)\s*\{',s,re.M)
        if not m: raise SystemExit(f"nao achei {velho} em {caminho}")
        corpo,fim=corpo_em(s,m.end()-1)
        ini=m.start(); k=s.rfind('\n',0,ini-1); com=''
        if k>=0 and s[k+1:ini].lstrip().startswith('//'):
            com=s[k+1:ini].strip(); ini=k+1
        ext[velho]=(novo,m.group(1),corpo,com)
        s=s[:ini]+s[fim+1:]
    met=[]
    for velho,novo in pares:
        nv,args,corpo,com=ext[velho]
        if com: met.append("    "+com)
        L=corpo.split('\n')
        met.append(f"    real {nv}({args}) {L[0]}")
        met+=[("    "+l) if l.strip() else "" for l in L[1:]]
    m=re.search(r'class (\w+) : ([^{]+)\{', s)
    if not m: raise SystemExit("nao achei a classe do exemplo")
    bases=m.group(2).strip()
    if base not in bases:
        s=s[:m.start()]+f"class {m.group(1)} : {bases}, public {base} {{"+s[m.end():]
    fimcls=s.index('\n};', m.start())
    s=s[:fimcls]+f"\n\n    // --- {rotulo} ---\n"+'\n'.join(met)+s[fimcls:]
    s=re.sub(r'\n{3,}','\n\n',s)
    r=re.search(rx_chamada, s, re.S)
    if not r: raise SystemExit(f"nao achei a chamada de registro ({rotulo}) em {caminho}")
    s=s[:r.start()]+nova_chamada+s[r.end():]
    open(caminho,'w').write(s)
    ok=sum(1 for a,b in pares if achar(orig,r'^real '+a+r'\([^)]*\)\s*\{')==achar(s,r'^\s*real '+b+r'\([^)]*\)\s*\{'))
    print(f"  {caminho.split('/')[-2]:<24} {rotulo:<26} {ok} de {len(pares)} corpos identicos")
    return ok==len(pares)

if __name__=='__main__':
    cfg=json.loads(sys.argv[1])
    sys.exit(0 if converte(**cfg) else 1)
