import re, sys
PARES=[("get_tensor","tensor"),("get_kernel","kernel"),
       ("get_kernel_inverse","kernel_inverse"),("get_kernel_jacobian","kernel_jacobian")]
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
def conv(caminho, registro=None):
    orig=open(caminho).read(); s=orig
    ext={}
    for velho,novo in PARES:
        m=re.search(r'^real '+velho+r'\(([^)]*)\)\s*\{',s,re.M)
        if not m: raise SystemExit(f"nao achei {velho}")
        corpo,fim=corpo_em(s,m.end()-1)
        ini=m.start(); k=s.rfind('\n',0,ini-1); com=''
        if k>=0 and s[k+1:ini].lstrip().startswith('//'):
            com=s[k+1:ini].strip(); ini=k+1
        ext[velho]=(novo,m.group(1),corpo,com)
        s=s[:ini]+s[fim+1:]
    met=[]
    for velho,novo in PARES:
        nv,args,corpo,com=ext[velho]
        if com: met.append("    "+com)
        L=corpo.split('\n')
        met.append(f"    real {nv}({args}) {L[0]}")
        met+=[("    "+l) if l.strip() else "" for l in L[1:]]
    # herda tambem a interface viscoelastica e injeta os metodos antes do '};'
    m=re.search(r'class (\w+) : public HigFlowProblem \{', s)
    if not m: raise SystemExit("nao achei a classe do exemplo")
    s=s[:m.start()]+f"class {m.group(1)} : public HigFlowProblem, public HigFlowViscoelasticProblem {{"+s[m.end():]
    fimcls=s.index('\n};', m.start())
    s=s[:fimcls]+"\n\n    // --- modelo viscoelastico ---\n"+'\n'.join(met)+s[fimcls:]
    s=re.sub(r'\n{3,}','\n\n',s)
    alvo=registro or caminho
    t=s if alvo==caminho else open(alvo).read()
    r=re.search(r'higflow_create_domain_viscoelastic\(ns,[^;]*?\);',t,re.S)
    if not r: raise SystemExit(f"nao achei create_domain_viscoelastic em {alvo}")
    args0=re.match(r'higflow_create_domain_viscoelastic\(ns,\s*([^,]+),\s*([^,]+),',r.group(0))
    t=t[:r.start()]+f"higflow_create_domain_viscoelastic(ns, {args0.group(1).strip()}, {args0.group(2).strip()}, &problema);"+t[r.end():]
    if alvo==caminho: s=t
    else: open(alvo,'w').write(t)
    open(caminho,'w').write(s)
    ok=sum(1 for a,b in PARES if achar(orig,r'^real '+a+r'\([^)]*\)\s*\{')==achar(s,r'^\s*real '+b+r'\([^)]*\)\s*\{'))
    print(f"  {caminho.split('/')[-2]}: {ok} de 4 corpos identicos")
    return ok==4
if __name__=='__main__':
    sys.exit(0 if conv(sys.argv[1], sys.argv[2] if len(sys.argv)>2 else None) else 1)
