# Evidence only. BS eq:betagiv graph sum (point vertices, rank 1, eq:contbetagiv) versus the single Hodge vertex Z[g,0] = Fpot_g(q_0)
# versus beta_{2g-2}/(-e1 e2 E^3)^(g-1). Tests the reconciliation by PPZ Prop 2.8, 2.9.
import sys, math, re, itertools
from fractions import Fraction as Fr
from collections import Counter
from hodge_lib import *
from fc import theta_series, zmul
from gap_series import beta
from fc4_graphs import graphs4
def run(g,e1,e2,sgn=1):
    E=e1+e2; N=3*g+4
    h=zmul(zmul(theta_series(-1/E,N),theta_series(1/e1,N),N),theta_series(1/e2,N),N)
    hm=[h[k]*(-1)**k for k in range(N+1)] # h(-z)
    if sgn==-1: hm=h
    gs=[Fr(1)*0]*(N+1)
    geo=[Fr(1)/E**k for k in range(N+1)]  # 1/(1-z/E)
    prod=zmul(hm,geo,N)
    dil=lambda d: -prod[d-1]+(1 if d==1 else 0)
    def edge(a,b):
        # (1-h(-z)h(-w))/(z+w): 1-H(z,w) vanishes on z=-w as h(z)h(-z)=1. coefficient by telescoping
        # H = sum c_M c_b z^M w^b; (z^M-(-w)^M)/(z+w) = sum_{k<M} (-1)^k z^{M-1-k} w^k
        tot=Fr(0)
        for k in range(b+1):
            M=a+k+1; bb=b-k
            tot-=hm[M]*hm[bb]*(-1)**k
        return tot
    def vert(gv,hs):
        hs=tuple(sorted(hs)); dim=3*gv-3+len(hs)-sum(hs)
        if dim<0: return Fr(0)
        tot=Fr(0)
        def parts(rem,maxp):
            yield ()
            for p in range(min(rem,maxp),0,-1):
                for r in parts(rem-p,p): yield (p+1,)+r
        for ds in parts(dim,dim):
            val=wk(gv,tuple(sorted(hs+ds)))
            if val==0: continue
            m=1
            for c in Counter(ds).values(): m*=math.factorial(c)
            x=val/m
            for d in ds: x*=dil(d)
            tot+=x
        return tot
    glist=stable_graphs(g) if g<=3 else graphs4()
    total=Fr(0)
    for gen,edges,aut in glist:
        nv=len(gen); val=[0]*nv
        for u,v in edges:
            val[u]+=1; val[v]+=1
        bud=[3*gen[i]-3+val[i] for i in range(nv)]
        hs=[[] for _ in range(nv)]
        acc=[Fr(0)]
        def rec(k,pref):
            if k==len(edges):
                x=pref
                for i in range(nv):
                    x*=vert(gen[i],tuple(hs[i]))
                    if x==0: return
                acc[0]+=x; return
            u,v=edges[k]
            bu=bud[u]-sum(hs[u])
            for a in range(bu+1):
                bv=(bud[v]-sum(hs[v])) if u!=v else (bu-a)
                for b in range(bv+1):
                    hs[u].append(a); hs[v].append(b)
                    rec(k+1,pref*edge(a,b))
                    hs[u].pop(); hs[v].pop()
        rec(0,Fr(1))
        total+=acc[0]/aut
    # single Hodge vertex
    H=hodge_cache([E,-e1,-e2])
    Z=Fr(0); dim=3*g-3
    def parts(rem,maxp):
        yield ()
        for p in range(min(rem,maxp),0,-1):
            for r in parts(rem-p,p): yield (p+1,)+r
    for ds in parts(dim,dim):
        val=H(g,tuple(sorted(ds)))
        if val==0: continue
        m=1
        for c in Counter(ds).values(): m*=math.factorial(c)
        x=val/m
        for d in ds: x*=(-Fr(1)/E**(d-1))
        Z+=x
    b=beta(2*g-2,e1,e2)/(-e1*e2*E**3)**(g-1)
    return total,Z,b
for g in (2,3,4):
    for e1,e2 in ((Fr(3),Fr(7)),(Fr(1),Fr(1)),(Fr(5,2),Fr(-3,7))):
        T,Z,b=run(g,e1,e2)
        Tc,_,_=run(g,e1,e2,-1) if g<4 else (None,0,0)
        print(g,e1,e2,"graphsum-Z=",T-Z,"Z-beta/(..)=",Z-b,"Z=",Z,"control(h(+z)) differs:",(Tc!=Z) if g<4 else "-",flush=True)
