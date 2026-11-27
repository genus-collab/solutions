# Evidence only. Genus-4 frame change using the graph list of the project (hodge_collect_genus4.txt), compared with stored hodge_F40.m.
import os
import re, sys, time
from fractions import Fraction as Fr
from fc import *
from gap_series import F40expr
def graphs4():
    out=[]
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_collect_genus4.txt")):
        m=re.match(r'^(\d+)\s+\{([\d, ]*)\}\s+\{(.*?)\}\s+(\d+)\s+\d+\s+\d+',line)
        if not m: continue
        gen=tuple(int(x) for x in m.group(2).replace(' ','').split(',') if x)
        edges=[(int(a)-1,int(b)-1) for a,b in re.findall(r'\{(\d+), *(\d+)\}',m.group(3))]
        out.append((gen,edges,int(m.group(4))))
    return out
gs=graphs4(); print(len(gs),"graphs; sum 1/Aut =",sum(Fr(1,g[2]) for g in gs),flush=True)
assert len(gs)==379
f4=F40expr()
def tot(th):
    t=W(0,0)
    for gen,edges,aut in gs:
        t=t+th.graph(gen,edges)*Fr(1,aut)
    return t
samples=[(Fr(3),Fr(7),Fr(1,5),Fr(1)),(Fr(1),Fr(1),Fr(-2,7),Fr(-2,3))]
N=18
for e1,e2,t,F in samples[:int(sys.argv[1])]:
    t0=time.time()
    rows=Rows(e1,e2,N); wts=[e1+e2,-e1,-e2]
    P=tot(Theory(rows,t,F,'pt',1)); t1=time.time()
    H=tot(Theory(rows,t,F,'hodge',0,wts))
    ref=eval(f4,{},dict(e1=e1,e2=e2,y=t**3,XX=F,Fr=Fr))
    print(f"g=4 eps=({e1},{e2}) t={t} F={F}: pt-hodge={(P-H).a},{(P-H).b} pt-stored={(P-ref).a},{P.b} value {P.a} [{t1-t0:.0f}s,{time.time()-t1:.0f}s]",flush=True)
