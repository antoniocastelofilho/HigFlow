#!/usr/bin/env python3
"""Monta um HTML AUTOCONTIDO com a animacao: interface e linhas de corrente.

O DESENHO E O ESQUEMA DE DADOS SAO OS DO MONTADOR DA SECAO COMPANHEIRA, e as
duas licoes que ele carrega vieram junto porque as duas foram apanhadas na
pagina RENDERIZADA, nao no fonte:

  `<meta charset>` E' OBRIGATORIO num arquivo autocontido.  Sem ele o navegador
  le' o UTF-8 como Latin-1 e o texto sai corrompido -- e o fonte parece certo;

  SIMBOLO NAO E' TEXTO.  O cabecalho vai em maiusculas por desenho, e
  `text-transform` troca `rho` por `Rho` e `mu` por `Mu` -- simbolos
  DIFERENTES, nao a mesma letra em outra caixa.  A classe `.sym` desliga a
  transformacao, e todo simbolo tem de ir dentro dela.

POR QUE UM ARQUIVO PROPRIO E NAO O REUSO DIRETO.  O montador de la' nao e' so'
desenho: ele AFIRMA coisas sobre os dados no texto fixo -- "seccao vertical
pelo centro da bolha", "nao linhas tridimensionais: a componente fora do plano
e' descartada", e o cabecalho "HiG-Flow".  As tres sao verdade no caso
tridimensional dele e FALSAS aqui, onde o dominio E' o plano e nao ha'
componente fora dele.

Reusa-lo verbatim publicaria tres afirmacoes falsas sobre o proprio resultado,
numa pagina que parece impecavel.  O esquema JSON continua o MESMO, de modo que
os dois lados seguem comparaveis -- o que muda e' so' o que a pagina diz.
"""

import json
import os
import sys

MODELO = """<meta charset="utf-8">
<title>%(tit)s</title>
<style>
:root{--paper:#FBFAF7;--ink:#1A1F26;--muted:#5A6472;--rule:#DFDBD3;
 --surface:#fff;--accent:#1D4E6F;--copper:#A9571F;--linha:#7A8694;
 --sans:"IBM Plex Sans","Segoe UI",system-ui,sans-serif;
 --mono:"IBM Plex Mono",Consolas,monospace;
 --serif:Spectral,Georgia,serif}
@media(prefers-color-scheme:dark){:root:not([data-theme=light]){
 --paper:#14181D;--ink:#E6E4DF;--muted:#9AA3B0;--rule:#2E353E;
 --surface:#1B2027;--accent:#74ADD1;--copper:#D98A4F;--linha:#5B6878}}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
 line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:1040px;margin:0 auto;padding:2.5rem 1.5rem 4rem}
.eyebrow{font-family:var(--mono);font-size:.7rem;letter-spacing:.14em;
 text-transform:uppercase;color:var(--copper);margin:0 0 .6rem}
/* SIMBOLO NAO E' TEXTO -- ver o cabecalho deste arquivo. */
.sym{text-transform:none !important}
h1{font-family:var(--serif);font-weight:700;font-size:clamp(1.6rem,4vw,2.3rem);
 line-height:1.15;margin:0 0 .6rem;letter-spacing:-.015em}
.lede{color:var(--muted);margin:0 0 2rem;max-width:62ch;font-size:.97rem}
.paineis{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));
 gap:1.25rem;margin-bottom:1.25rem}
.painel{background:var(--surface);border:1px solid var(--rule);border-radius:3px;
 padding:.9rem}
.painel h2{font-family:var(--mono);font-size:.74rem;letter-spacing:.1em;
 text-transform:uppercase;color:var(--muted);margin:0 0 .6rem;font-weight:500}
canvas{width:100%%;height:auto;display:block}
.ctl{display:flex;gap:.8rem;align-items:center;flex-wrap:wrap;margin-bottom:1.1rem}
button{font-family:var(--mono);font-size:.8rem;padding:.4rem .9rem;cursor:pointer;
 background:var(--accent);color:#fff;border:0;border-radius:3px}
button.sec{background:transparent;color:var(--accent);
 border:1px solid var(--accent)}
input[type=range]{flex:1;min-width:140px;accent-color:var(--accent)}
.leg{font-family:var(--mono);font-size:.76rem;color:var(--muted)}
.nota{font-size:.87rem;color:var(--muted);max-width:72ch;
 border-left:3px solid var(--copper);padding:.7rem 1rem;background:var(--surface)}
.nota b{color:var(--ink)}
</style>
<div class="wrap">
<p class="eyebrow">t8Flow &middot; <span class="sym">%(sub)s</span></p>
<h1>%(tit)s</h1>
<p class="lede">%(lede)s</p>

<div class="ctl">
  <button id="pp">Pausar</button>
  <button class="sec" id="rw">Início</button>
  <input type="range" id="sl" min="0" max="0" value="0">
  <span class="leg" id="rot"></span>
</div>

<div class="paineis" id="paineis"></div>
<p class="nota" id="nota"></p>
</div>
<script>
const D = %(js)s;
const ROT = %(rot)s;
const DTQ = %(dtq)s, CAIXA = %(caixa)s, NOTA = %(nota)s;
function css(v){return getComputedStyle(document.documentElement)
  .getPropertyValue(v).trim()}

const cont = document.getElementById('paineis');
const cvs = {};
ROT.forEach(r=>{
  const d=document.createElement('div'); d.className='painel';
  d.innerHTML='<h2><span class="sym">'+r+'</span></h2><canvas></canvas>';
  cont.appendChild(d);
  cvs[r]=d.querySelector('canvas');
});

const nq = Math.min(...ROT.map(r=>D[r].quadros.length));
document.getElementById('sl').max = nq-1;

function desenha(r,k){
  const c=cvs[r], q=D[r].quadros[k];
  const dpr=window.devicePixelRatio||1;
  const L=c.clientWidth||320, alt=Math.round(L*CAIXA);
  c.width=L*dpr; c.height=alt*dpr;
  const g=c.getContext('2d'); g.setTransform(dpr,0,0,dpr,0,0);
  g.clearRect(0,0,L,alt);
  const X=x=>x*L, Y=y=>alt-(y/CAIXA)*alt;
  g.strokeStyle=css('--rule'); g.lineWidth=1; g.strokeRect(.5,.5,L-1,alt-1);
  if(q.l){
    g.strokeStyle=css('--linha'); g.lineWidth=.7; g.globalAlpha=.75;
    q.l.forEach(P=>{ g.beginPath();
      P.forEach((p,i)=> i? g.lineTo(X(p[0]),Y(p[1])) : g.moveTo(X(p[0]),Y(p[1])));
      g.stroke(); });
    g.globalAlpha=1;
  }
  g.strokeStyle=css('--accent'); g.lineWidth=2.2;
  q.c.forEach(P=>{ g.beginPath();
    P.forEach((p,i)=> i? g.lineTo(X(p[0]),Y(p[1])) : g.moveTo(X(p[0]),Y(p[1])));
    g.closePath(); g.stroke(); });
}

let k=0, tocando=true;
function pinta(){
  ROT.forEach(r=>desenha(r,k));
  document.getElementById('sl').value=k;
  const q=D[ROT[0]].quadros[k];
  document.getElementById('rot').textContent =
    't = '+(q.q*DTQ).toFixed(3)+'   (quadro '+q.q+')';
  // O QUE FALTA E' DITO NO QUADRO, e nao preenchido: um quadro sem campo
  // desenhado como os outros faria parecer que o escoamento parou.
  const semcampo = ROT.some(r=>!D[r].quadros[k].l);
  document.getElementById('nota').innerHTML =
    semcampo ? '<b>Sem campo de velocidade neste quadro.</b> O contorno é '
             + 'a interface; as linhas não existem aqui porque o campo não foi '
             + 'gravado neste instante.'
             : NOTA;
}
document.getElementById('sl').oninput=e=>{k=+e.target.value; tocando=false;
  document.getElementById('pp').textContent='Tocar'; pinta();};
document.getElementById('pp').onclick=e=>{tocando=!tocando;
  e.target.textContent=tocando?'Pausar':'Tocar';};
document.getElementById('rw').onclick=()=>{k=0; pinta();};
window.addEventListener('resize', pinta);
setInterval(()=>{ if(tocando){ k=(k+1)%%nq; pinta(); } }, 120);
pinta();
</script>
"""


def monta(saida, titulo, sub, pares, dtq, caixa, lede, nota):
    dados = {rot: json.load(open(cam)) for rot, cam in pares}
    html = MODELO % {"tit": titulo, "sub": sub, "lede": lede,
                     "js": json.dumps(dados, separators=(",", ":")),
                     "rot": json.dumps(list(dados.keys())),
                     "dtq": dtq, "caixa": caixa,
                     "nota": json.dumps(nota)}
    open(saida, "w").write(html)
    return os.path.getsize(saida)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("uso: monta-animacao.py <saida.html> <dt por quadro> "
              "<rotulo>=<arq.json> [<rotulo>=<arq.json> ...]")
        sys.exit(2)
    saida, dtq = sys.argv[1], float(sys.argv[2])
    pares = [(a.split("=", 1)[0], a.split("=", 1)[1]) for a in sys.argv[3:]]
    n = monta(saida,
              "Bolha de Hysing no t8Flow",
              "caso 1 &middot; &rho; 1000/100 &middot; &mu; 10/1 &middot; &sigma; 24,5 &middot; Re 35 &middot; Eo 10",
              pares, dtq, 2.0,
              # O TEXTO DA PAGINA VAI ACENTUADO.  A convencao sem acento
              # deste projeto e' dos COMENTARIOS de codigo, onde evita
              # problemas de codificacao em arquivo-fonte; uma pagina que
              # alguem le' nao herda essa restricao -- e o `<meta charset>`
              # acima existe justamente para que ela nao precise herdar.
              "O domínio <b>é</b> o plano: não há secção nem componente fora "
              "dele. O contorno é a interface que o front-tracking transporta; "
              "as linhas são do campo de velocidade <b>deste instante</b>, "
              "integradas por Runge–Kutta de segunda ordem nos dois sentidos a "
              "partir de cada semente.",
              "Linhas de corrente do campo <b>deste quadro</b>, não "
              "trajetórias de partícula: num escoamento que muda no tempo as "
              "duas coisas são diferentes.")
    print("%s  %.2f MB" % (saida, n / 1e6))
