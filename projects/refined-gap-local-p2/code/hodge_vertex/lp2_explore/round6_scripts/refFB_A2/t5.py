# genus 4 via degeneration graph generator; tests z^7 constant (a_4), d=0,1 (mod nothing: exact)
from gsum import *
import gsum, engine, time
from itertools import permutations, combinations
from math import factorial

def canon(gv,E):
    V=len(gv)
    # E multiset of (u,v) u<=v
    from collections import Counter
    cnt=Counter((min(a,b),max(a,b)) for a,b in E)
    best=None; auts=0; key0=None
    # prune perms by genus/valence class
    val=[sum((a==v)+(b==v) for a,b in E) for v in range(V)]
    lab=[(gv[v],val[v],cnt.get((v,v),0)) for v in range(V)]
    order=sorted(range(V),key=lambda v:lab[v])
    sig=tuple(lab[v] for v in order)
    groups=[]
    i=0
    while i<V:
        j=i
        while j<V and lab[order[j]]==lab[order[i]]: j+=1
        groups.append(order[i:j]); i=j
    def perms(gs):
        if not gs:
            yield []; return
        for p in permutations(gs[0]):
            for q in perms(gs[1:]):
                yield list(p)+q
    pos0=None
    for p in perms(groups):
        idx={v:i for i,v in enumerate(p)}
        k=tuple(sorted(((min(idx[a],idx[b]),max(idx[a],idx[b])),m) for (a,b),m in cnt.items()))
        if best is None or k<best: best=k
    # aut count
    n=0
    for p in perms(groups):
        idx={v:i for i,v in enumerate(p)}
        k=tuple(sorted(((min(idx[a],idx[b]),max(idx[a],idx[b])),m) for (a,b),m in cnt.items()))
        if k==best: n+=1
    aut=n
    for (a,b),m in cnt.items():
        aut*=factorial(m)
        if a==b: aut*=2**m
    return (sig,best),aut

def gen(g):
    seen={}
    start=([g],[])
    k0,a0=canon([g],[])
    seen[k0]=dict(gv=[g],edges=[],aut=a0,legs=[0])
    frontier=[start]
    while frontier:
        new=[]
        for gv,E in frontier:
            V=len(gv)
            for v in range(V):
                # (i) loop
                if gv[v]>=1:
                    g2=list(gv); g2[v]-=1; E2=E+[(v,v)]
                    k,a=canon(g2,E2)
                    if k not in seen:
                        seen[k]=dict(gv=g2,edges=E2,aut=a,legs=[0]*V); new.append((g2,E2))
                # (ii) split: half-edges at v
                hes=[]
                for ei,(a_,b_) in enumerate(E):
                    if a_==v: hes.append((ei,0))
                    if b_==v: hes.append((ei,1))
                n=len(hes)
                for mask in range(2**n):
                    for g1 in range(gv[v]+1):
                        g2_=gv[v]-g1
                        n1=bin(mask).count('1'); n2=n-n1
                        if 2*g1-2+n1+1<=0 or 2*g2_-2+n2+1<=0: continue
                        gnew=list(gv); gnew[v]=g1; gnew.append(g2_)
                        Enew=[list(e) for e in E]
                        for t,(ei,side) in enumerate(hes):
                            if mask>>t&1:
                                pass
                            else:
                                Enew[ei][side]=V
                        Enew=[tuple(e) for e in Enew]+[(v,V)]
                        k,a=canon(gnew,Enew)
                        if k not in seen:
                            seen[k]=dict(gv=gnew,edges=Enew,aut=a,legs=[0]*(V+1)); new.append((gnew,Enew))
        frontier=new
    return list(seen.values())

if __name__=="__main__":
    for g in (2,3):
        G=gen(g); print(g,len(G),sum(Fr(1,x['aut']) for x in G), sum(Fr(1,x['aut']) for x in graphs(g,0)))
    t=time.time()
    G4=gen(4); print("g=4 graphs",len(G4),time.time()-t,flush=True)
    orig=gsum.graphs
    gsum.graphs=lambda g,n: G4
    al,ep=(Fr(1),Fr(2),Fr(5)),(Fr(3),Fr(7))
    c=setup(al,ep,1,10); gs=GS(c)
    F=gs.total(4); print('g4 graph sum',[str(x) for x in F],time.time()-t,flush=True)
    R=q_to_series(c,[loc(al,ep,4,d) for d in range(2)])
    print('loc        ',[str(x) for x in R],'MATCH',R==F,time.time()-t,flush=True)
