import os
import re,time
from fractions import Fraction as Fr
from tf import Closed,TF,Om
from cmp_loc import series,evalser,Qofy
def graphs4():
    out=[]
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_collect_genus4.txt")):
        m=re.match(r'^(\d+)\s+\{([\d, ]*)\}\s+\{(.*?)\}\s+(\d+)\s+\d+\s+\d+',line)
        if not m: continue
        gen=tuple(int(x) for x in m.group(2).replace(' ','').split(',') if x)
        edges=[(int(a)-1,int(b)-1) for a,b in re.findall(r'\{(\d+), *(\d+)\}',m.group(3))]
        out.append((gen,edges,int(m.group(4))))
    return out
gs=graphs4()
e1,e2=3,7; E=Fr(10)
cl=Closed(E,16)
vals=[Fr(625,6048),Fr(52397281,1209600),Fr(-81799253033,282240)]
a,F,gsr=series(9)
for t in (Fr(1,500),Fr(-1,500)):
    t0=time.time()
    y=t**3; Fv=evalser(F,y); Av=evalser(a,y); Q=Qofy(y,gsr,6)
    th=TF(cl,e1,e2,t,Fv,Av,'pt',1)
    P=Om(0)
    for gen,ed,aut in gs: P=P+th.graph(gen,ed)*Fr(1,aut)
    print(f"t={t}: (Tbar_4-GW_0-Q GW_1)/Q^2={float((P.a-vals[0]-Q*vals[1])/Q**2):.6f}; GW_(4,2)={float(vals[2]):.6f}; Im={float(P.b):.1e} [{time.time()-t0:.0f}s]",flush=True)
