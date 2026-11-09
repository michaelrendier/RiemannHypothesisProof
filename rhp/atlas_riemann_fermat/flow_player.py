# SPDX-License-Identifier: GPL-3.0-only
"""flow_player.html — watch the extinction flow: drag τ = ln y (or press play) and see each integer die (Telperion) and be born (Laurelin).

Self-contained page (no network).  Panel 1: every n ≤ N at (u = ln n, status).  Panel 2: the survivor fraction against Mertens
(e^{−γ}/τ) with a marker at the current τ.  Colours: blue = still alive (spf > y), grey = dead (spf ≤ y), orange ring = smooth (gpf ≤ y).
Primes are black: they die/are born at the same instant, u = τ.   Code drawn: rhp.wrap.spf_gpf.   Run: python -m rhp.atlas_riemann_fermat.flow_player
"""
import json
import math
from ._common import HERE
from .. import wrap as W

N = 2500

HTML = r"""<!doctype html><meta charset=utf-8><title>Extinction flow</title>
<style>body{font:13px monospace;margin:12px;background:#fff;color:#1b1b1f}canvas{border:1px solid #ccc;margin-right:8px}
@media(prefers-color-scheme:dark){body{background:#15151a;color:#ddd}canvas{border-color:#444}}</style>
<h3>The extinction flow: Telperion (sieve) and Laurelin (un-sieve) as τ = ln y runs</h3>
<p>τ <input id=s type=range min=0 max=1000 value=0 style="width:420px"> <b id=lab></b> <button id=pl>play</button></p>
<canvas id=c1 width=640 height=420></canvas><canvas id=c2 width=420 height=420></canvas>
<p id=ro></p><p>blue = alive (spf &gt; y) · grey = dead (spf ≤ y) · orange ring = smooth (gpf ≤ y) · black = prime · dashed = e<sup>−γ</sup>/τ (Mertens)</p>
<script>
const D=__DATA__, N=D.N, spf=D.spf, gpf=D.gpf, tmax=Math.log(N);
const c1=document.getElementById('c1').getContext('2d'), c2=document.getElementById('c2').getContext('2d');
const s=document.getElementById('s'), lab=document.getElementById('lab'), ro=document.getElementById('ro');
function mer(y){let r=1;for(const p of D.primes){if(p>y)break;r*=1-1/p}return r}
function draw(){
 const tau=+s.value/1000*tmax, y=Math.exp(tau); lab.textContent='τ='+tau.toFixed(2)+'  y='+y.toFixed(1);
 c1.clearRect(0,0,640,420); let alive=1,smooth=1;
 for(let n=2;n<=N;n++){const u=Math.log(n)/tmax*620+10, isP=spf[n]===n;
  const dead=spf[n]<=y, sm=gpf[n]<=y; if(!dead)alive++; if(sm)smooth++;
  const row=Math.log(spf[n])/tmax; const yy=410-row*390;                    // vertical = ln spf (death clock)
  c1.fillStyle=isP?'#000':(dead?'#bbb':'#2a6fb0'); c1.fillRect(u,yy,2,2);
  if(sm&&!isP){c1.strokeStyle='#c9732a';c1.strokeRect(u-1,yy-1,4,4)}}
 c1.strokeStyle='#c0392b';c1.beginPath();c1.moveTo(0,410-tau/tmax*390);c1.lineTo(640,410-tau/tmax*390);c1.stroke();
 c1.fillStyle='#888';c1.fillText('u = ln n  →',560,415);c1.fillText('ln spf(n)',4,12);
 c2.clearRect(0,0,420,420);c2.strokeStyle='#888';c2.beginPath();c2.moveTo(30,10);c2.lineTo(30,400);c2.lineTo(410,400);c2.stroke();
 c2.setLineDash([4,3]);c2.beginPath();for(let k=0;k<=200;k++){const t=0.7+k/200*(tmax-0.7),x=30+t/tmax*380,v=Math.exp(-0.5772)/t;
  const yy=400-Math.min(1,v)*380;k?c2.lineTo(x,yy):c2.moveTo(x,yy)}c2.stroke();c2.setLineDash([]);
 c2.strokeStyle='#1b1b1f';c2.beginPath();for(let k=0;k<=200;k++){const t=k/200*tmax,v=mer(Math.exp(t));const x=30+t/tmax*380,yy=400-v*380;k?c2.lineTo(x,yy):c2.moveTo(x,yy)}c2.stroke();
 c2.fillStyle='#c0392b';c2.beginPath();c2.arc(30+tau/tmax*380,400-alive/N*380,5,0,7);c2.fill();
 c2.fillStyle='#888';c2.fillText('survivor fraction vs τ',40,20);
 ro.textContent='alive '+alive+' / '+N+'  = '+(alive/N).toFixed(4)+'   Mertens Π(1−1/p) = '+mer(y).toFixed(4)+'   smooth (born) '+smooth+' / '+N+'   U = ln N/ln y = '+(tau>0?(tmax/tau).toFixed(2):'∞');}
s.oninput=draw; let tm=null; document.getElementById('pl').onclick=()=>{if(tm){clearInterval(tm);tm=null;return}
 tm=setInterval(()=>{s.value=(+s.value+4)%1001;draw()},40)}; draw();
</script>"""


def build():
    spf, gpf = W.spf_gpf(N)
    data = {"N": N, "spf": spf[:N + 1].tolist(), "gpf": gpf[:N + 1].tolist(), "primes": W.primes_to(N).tolist()}
    out = HERE / "flow_player.html"
    out.write_text(HTML.replace("__DATA__", json.dumps(data)))
    return out


if __name__ == "__main__":
    print(build())
