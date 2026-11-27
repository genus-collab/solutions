# Fresh code: all stable graphs of genus 2 and 3 (list read from the stored percolour logs), every colouring,
# compare the colour-and-graph sum of (0.1) with the stored F20, F30 (hodge_F20.m, hodge_F30.m) as series in t_c. Evidence only.
import os
import re,sys,time
from fractions import Fraction as Fr
from itertools import product
from rows import *
from graphs_check import Vtx,graph_value,recon
LP=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "")
CG=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus4_hodge", "")
def read_graphs(fn):
    out=[]
    seen=set()
    for line in open(fn,encoding="utf-8",errors="replace"):
        m=re.search(r"graph (\d+): genera \{([^}]*)\}, edges \{(.*)\}, Aut (\d+)",line)
        if m:
            k=int(m.group(1))
            if k in seen: continue
            seen.add(k)
            gen=[int(x) for x in m.group(2).split(",")]
            eds=re.findall(r"\{(\d+), (\d+)\}",m.group(3))
            edges=[(int(a)-1,int(b)-1) for a,b in eds]
            out.append((k,gen,edges,int(m.group(4))))
    return out
def parse_expr(fn,e1,e2,y,XX,prec):
    s=open(fn,encoding="utf-8").read()
    s=s.replace("\\[Epsilon][1]","e1").replace("\\[Epsilon][2]","e2").replace("\n"," ")
    s=re.sub(r"\^","**",s)
    # implicit multiplication by whitespace does not occur? handle ' ' between factors
    # wrap integers not preceded by ** or letters
    s=re.sub(r"(?<![\*\w])(\d+)",lambda m:"Z(%s)"%m.group(1),s)
    s=re.sub(r"\)\s+(?=[A-Za-z(])",")*",s)
    s=re.sub(r"(\w)\s+(?=[A-Za-z(])",r"\1*",s)
    Z=lambda n:L.const(n,prec)
    env=dict(Z=Z,e1=L.const(fp(e1),prec),e2=L.const(fp(e2),prec),y=y,XX=XX)
    return eval(s,{"__builtins__":{}},env)
def main():
    e1,e2=Fr(sys.argv[1]),Fr(sys.argv[2]); g=int(sys.argv[3])
    th=Theory(e1,e2,N=34,K=13)
    V=Vtx(th)
    y=(th.Del-1)*inv(27)
    # XX = theta log c1, with theta c1 computed from t_c-derivative
    XX=th.theta(th.c1)/th.c1
    print("XX(0)=",recon(XX.coef(0)))
    graphs=read_graphs(LP+"lp2_percolour_g%d_out.txt"%g)
    print("graphs",len(graphs))
    tot=None; t0=time.time()
    bad=[]
    for k,gen,edges,aut in graphs:
        gs=None
        for col in product(range(3),repeat=len(gen)):
            val=graph_value(th,V,gen,edges,col) if edges else V.vertex(gen[0],col[0],())
            nv=val.norm()
            if edges and nv.v<0 and nv.c: bad.append((k,col,nv.v))
            gs=val if gs is None else gs+val
        gs=gs*inv(aut)
        tot=gs if tot is None else tot+gs
    print("time",round(time.time()-t0,1),"edge graphs with a pole for some colouring:",bad)
    F=parse_expr(CG+"hodge_F%d0.m"%g,e1,e2,y,XX,th.prec)
    print("fresh G_g valuation",tot.norm().v,"stored F valuation",F.norm().v)
    d=tot-F
    print("fresh G_g coefs t^-(2g-2)..t^3:",[str(recon(tot.coef(k))) for k in range(2-2*g,4)])
    print("difference with stored:",[d.coef(k) for k in range(2-2*g,6)])
main()
