# Genus 4 (379 graphs from the project's list; enumeration is not the object under test). Evidence only.
import os
import re,time,sys
from fractions import Fraction as Fr
from tf import Closed,TF,Om
from cmp_loc import series,evalser,Qofy
from gapnum import cf_series,ev,dev,rnd,beta
def graphs4():
    out=[]
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_collect_genus4.txt")):
        m=re.match(r'^(\d+)\s+\{([\d, ]*)\}\s+\{(.*?)\}\s+(\d+)\s+\d+\s+\d+',line)
        if not m: continue
        gen=tuple(int(x) for x in m.group(2).replace(' ','').split(',') if x)
        edges=[(int(a)-1,int(b)-1) for a,b in re.findall(r'\{(\d+), *(\d+)\}',m.group(3))]
        out.append((gen,edges,int(m.group(4))))
    return out
gs=graphs4(); print(len(gs),sum(Fr(1,x[2]) for x in gs),flush=True)
def tot(th):
    t=Om(0)
    for gen,ed,aut in gs: t=t+th.graph(gen,ed)*Fr(1,aut)
    return t
e1,e2=3,7; E=Fr(10); wts=[E,-Fr(3),-Fr(7)]
cl=Closed(E,16)
vals=[Fr(625,6048),Fr(52397281,1209600),Fr(-81799253033,282240)]
a,F,gsr=series(30)
w,tc=cf_series(60)
# gap
pc=(-3)**3*beta(6,e1,e2)
print("polar coeff",pc,flush=True)
for d in (Fr(1,100),Fr(1,1000),Fr(1,10000)):
    t0=time.time()
    t=Fr(-1,3)+d; xi=1+27*t**3
    wv=rnd(ev(w,xi)); tcv=rnd(ev(tc,xi)); dw=rnd(dev(w,xi)); Fcf=rnd((xi-1)*dw/wv)
    H=tot(TF(cl,e1,e2,t,Fcf,-wv,'hv',0,wts))
    print(f"delta={d}: G_4-polar={float(H.a-pc/tcv**6):.10f} Im={float(H.b):.1e} [{time.time()-t0:.0f}s]",flush=True)
# large radius, point frame vs localisation
for t in (Fr(1,500),Fr(-1,500)):
    t0=time.time()
    y=t**3; Fv=evalser(F,y); Av=evalser(a,y); Q=Qofy(y,gsr)
    P=tot(TF(cl,e1,e2,t,Fv,Av,'pt',1))
    print(f"t={t}: (Tbar_4-GW_0-Q GW_1)/Q^2={float((P.a-vals[0]-Q*vals[1])/Q**2):.6f}; GW_(4,2)={float(vals[2]):.6f}; Im={float(P.b):.1e} [{time.time()-t0:.0f}s]",flush=True)
