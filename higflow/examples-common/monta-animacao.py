#!/usr/bin/env python3
"""Monta um HTML autocontido com a animacao: contorno da interface e linhas de
corrente, dois metodos em sincronia.  Canvas, sem dependencia externa -- nao ha'
ffmpeg, convert, gnuplot nem matplotlib nesta maquina, e um HTML abre em
qualquer navegador."""
import json, sys, os

def carrega(p):
    d = json.load(open(p)); return d

def monta(saida, titulo, sub, pares, dtq, caixa):
    dados = {}
    for rot, cam in pares:
        dados[rot] = carrega(cam)
    js = json.dumps(dados, separators=(',',':'))
    rotulos = list(dados.keys())
    # <meta charset>: OBRIGATORIO num arquivo AUTOCONTIDO.  Sem ele o
    # navegador le\' o UTF-8 como Latin-1 e o texto sai como "SA~ao a
    # projeA~gA~ao" -- apanhado olhando a pagina renderizada, nao o fonte, que
    # parecia certo.
    html = """<meta charset="utf-8">
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
/* SIMBOLO NAO E' TEXTO.  O cabecalho e o rotulo dos paineis vao em maiusculas
   por desenho, e isso trocava lambda por Lambda, beta por Beta e h por H --
   simbolos DIFERENTES, nao a mesma letra em outra caixa.  Apanhado lendo a
   pagina RENDERIZADA; no fonte os caracteres estavam certos. */
.sym{text-transform:none !important}
h1{font-family:var(--serif);font-weight:700;font-size:clamp(1.6rem,4vw,2.3rem);
 line-height:1.15;margin:0 0 .6rem;letter-spacing:-.015em}
.lede{color:var(--muted);margin:0 0 2rem;max-width:62ch;font-size:.97rem}
.paineis{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));
 gap:1.25rem;margin-bottom:1.25rem}
.painel{background:var(--surface);border:1px solid var(--rule);border-radius:3px;
 padding:.75rem}
.painel h2{font-family:var(--mono);font-size:.72rem;letter-spacing:.1em;
 text-transform:uppercase;color:var(--muted);margin:0 0 .5rem;font-weight:500}
canvas{width:100%%;height:auto;display:block}
.ctl{display:flex;gap:.9rem;align-items:center;flex-wrap:wrap;
 background:var(--surface);border:1px solid var(--rule);border-radius:3px;
 padding:.7rem .9rem;margin-bottom:1.5rem}
button{font-family:var(--sans);font-size:.88rem;font-weight:500;cursor:pointer;
 background:var(--accent);color:var(--paper);border:0;border-radius:3px;
 padding:.4rem .9rem}
button.sec{background:transparent;color:var(--accent);
 border:1px solid var(--accent)}
input[type=range]{flex:1;min-width:140px;accent-color:var(--accent)}
.leg{font-family:var(--mono);font-size:.76rem;color:var(--muted)}
.nota{font-size:.87rem;color:var(--muted);max-width:72ch;
 border-left:3px solid var(--copper);padding:.7rem 1rem;background:var(--surface)}
.nota b{color:var(--ink)}
</style>
<div class="wrap">
<p class="eyebrow">HiG-Flow &middot; %(sub)s</p>
<h1>%(tit)s</h1>
<p class="lede">Secção vertical pelo centro da bolha. O contorno é a interface
que cada método entrega; as linhas são a projeção do escoamento no plano de
corte, integradas no campo de velocidade.</p>

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
const DTQ = %(dtq)s, CAIXA = %(caixa)s;
const cores = {ink:'#1A1F26', linha:'#7A8694', acc:'#1D4E6F'};
function css(v){return getComputedStyle(document.documentElement).getPropertyValue(v).trim()}

const cont = document.getElementById('paineis');
const cvs = {};
ROT.forEach(r=>{
  const d=document.createElement('div'); d.className='painel';
  d.innerHTML='<h2>'+r+'</h2><canvas></canvas>';
  cont.appendChild(d);
  const c=d.querySelector('canvas');
  cvs[r]=c;
});

const nq = Math.min(...ROT.map(r=>D[r].quadros.length));
document.getElementById('sl').max = nq-1;

function desenha(r, k){
  const c=cvs[r], q=D[r].quadros[k];
  const dpr=window.devicePixelRatio||1;
  const L=c.clientWidth||320;
  // a caixa e' 1 x CAIXA (altura); escala para caber
  const alt=Math.round(L*CAIXA);
  c.width=L*dpr; c.height=alt*dpr;
  const g=c.getContext('2d'); g.setTransform(dpr,0,0,dpr,0,0);
  g.clearRect(0,0,L,alt);
  const X=x=>x*L, Z=z=>alt-(z/CAIXA)*alt;
  // moldura
  g.strokeStyle=css('--rule'); g.lineWidth=1;
  g.strokeRect(.5,.5,L-1,alt-1);
  // linhas de corrente
  if(q.l){
    g.strokeStyle=css('--linha'); g.lineWidth=.7; g.globalAlpha=.75;
    q.l.forEach(L2=>{ g.beginPath();
      L2.forEach((p,i)=> i? g.lineTo(X(p[0]),Z(p[1])) : g.moveTo(X(p[0]),Z(p[1])));
      g.stroke(); });
    g.globalAlpha=1;
  }
  // interface
  g.strokeStyle=css('--accent'); g.lineWidth=2.2;
  q.c.forEach(L2=>{ g.beginPath();
    L2.forEach((p,i)=> i? g.lineTo(X(p[0]),Z(p[1])) : g.moveTo(X(p[0]),Z(p[1])));
    g.closePath(); g.stroke(); });
}

let k=0, tocando=true;
function pinta(){
  ROT.forEach(r=>desenha(r,k));
  const q=D[ROT[0]].quadros[k];
  document.getElementById('rot').textContent =
    't = '+(q.q*DTQ).toFixed(3)+'   (quadro '+q.q+')';
  document.getElementById('sl').value=k;
  const semcampo = ROT.some(r=>!D[r].quadros[k].l);
  document.getElementById('nota').innerHTML = semcampo
    ? '<b>Sem linhas de corrente neste quadro.</b> O campo de velocidade no volume '
      +'só foi gravado em alguns instantes desta corrida, e interpolar entre eles '
      +'para preencher a animação mostraria um escoamento que não foi calculado.'
    : '<b>Linhas de corrente do campo deste quadro.</b> São a projeção no plano de '
      +'corte, não linhas tridimensionais: a componente fora do plano é descartada.';
}
document.getElementById('pp').onclick=e=>{tocando=!tocando;
  e.target.textContent=tocando?'Pausar':'Tocar'};
document.getElementById('rw').onclick=()=>{k=0;pinta()};
document.getElementById('sl').oninput=e=>{k=+e.target.value;tocando=false;
  document.getElementById('pp').textContent='Tocar';pinta()};
window.addEventListener('resize',pinta);
pinta();
setInterval(()=>{ if(tocando){ k=(k+1)%%nq; pinta(); } }, 160);
</script>
""" % {"tit":titulo, "sub":sub, "js":js, "rot":json.dumps(rotulos),
       "dtq":dtq, "caixa":caixa}
    open(saida,"w").write(html)
    return os.path.getsize(saida)

if __name__ == "__main__":
    S = sys.argv[1]
    H = sys.argv[2]
    n1 = monta(os.path.join(H,"anim-bolha-newtoniana.html"),
          "Bolha newtoniana: FT contra VOF",
          "caso 1, malha 25&times;25&times;50, t = 0 a 3",
          [("front-tracking", S+"/anim-newt-ft.json"),
           ("VOF", S+"/anim-newt-vof.json")], 0.025, 2.0)
    n2 = monta(os.path.join(H,"anim-bolha-viscoelastica.html"),
          "Bolha viscoelástica: FT contra VOF",
          "Oldroyd-B, &lambda;=0,4 &beta;=0,2 (a cúspide), t = 0 a 4",
          [("front-tracking", S+"/anim-ft-l04-b02.json"),
           ("VOF", S+"/anim-vof-l04-b02.json")], 0.025, 3.0)
    n3 = monta(os.path.join(H,"anim-bolha-ve-contra-newt.html"),
          "Viscoelástico contra newtoniano, mesma geometria",
          "front-tracking, caixa 1&times;1&times;3, t = 0 a 4",
          [("Oldroyd-B &lambda;=0,4 &beta;=0,2", S+"/anim-ft-l04-b02.json"),
           ("newtoniano", S+"/anim-ft-newt.json")], 0.025, 3.0)
    for n,f in ((n1,"anim-bolha-newtoniana.html"),(n2,"anim-bolha-viscoelastica.html"),
                (n3,"anim-bolha-ve-contra-newt.html")):
        print("  %-38s %.1f MB" % (f, n/1e6))
