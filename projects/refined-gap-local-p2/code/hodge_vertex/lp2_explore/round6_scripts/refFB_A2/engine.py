# Independent engines: WK numbers, Mumford ch-integrals, stable graph generator.
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import permutations, product, combinations
from math import factorial, comb
from collections import Counter
import sys
sys.setrecursionlimit(100000)

def dfact(n):
    r=1
    while n>1:
        r*=n; n-=2
    return r

def stable(g,n): return g>=0 and 2*g-2+n>0

@lru_cache(maxsize=None)
def wk(g,a):
    a=tuple(sorted(a)); n=len(a)
    if not stable(g,n): return Fr(0)
    if sum(a)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and n==1: return Fr(1,24)
    k1=a[-1]; rest=a[:-1]; k=k1-1
    res=Fr(0)
    for j,dj in enumerate(rest):
        if dj+k<0: continue
        b=list(rest); b[j]=dj+k
        res+=Fr(dfact(2*k+2*dj+1),dfact(2*dj-1))*wk(g,tuple(sorted(b)))
    for r in range(0,k):
        s=k-1-r
        c=dfact(2*r+1)*dfact(2*s+1)
        if g>=1:
            res+=Fr(c,2)*wk(g-1,tuple(sorted(rest+(r,s))))
        idx=list(range(len(rest)))
        for m in range(len(rest)+1):
            for I in combinations(idx,m):
                J=[i for i in idx if i not in I]
                aI=tuple(rest[i] for i in I); aJ=tuple(rest[i] for i in J)
                for g1 in range(0,g+1):
                    g2=g-g1
                    if not stable(g1,len(aI)+1) or not stable(g2,len(aJ)+1): continue
                    if sum(aI)+r!=3*g1-3+len(aI)+1: continue
                    res+=Fr(c,2)*wk(g1,tuple(sorted(aI+(r,))))*wk(g2,tuple(sorted(aJ+(s,))))
    return res/dfact(2*k+3)

_bc={}
def bern(n):
    if n in _bc: return _bc[n]
    B=[Fr(0)]*(n+1); B[0]=Fr(1)
    for m in range(1,n+1):
        B[m]=-sum(comb(m+1,k)*B[k] for k in range(m))/(m+1)
    _bc[n]=B[n]
    return B[n]

@lru_cache(maxsize=None)
def chint(g,a,chs):
    """int_{Mbar_{g,n}} prod psi^a prod ch_k(E), chs tuple of odd k, via Mumford."""
    a=tuple(sorted(a)); chs=tuple(sorted(chs)); n=len(a)
    if not stable(g,n): return Fr(0)
    if sum(a)+sum(chs)!=3*g-3+n: return Fr(0)
    if not chs: return wk(g,a)
    if g==0: return Fr(0)
    m=chs[-1]; rest=chs[:-1]
    assert m%2==1
    K=m+1
    pref=bern(K)/factorial(K)
    res=chint(g,a+(K,),rest)
    for i in range(n):
        b=list(a); b[i]+=m
        res-=chint(g,tuple(b),rest)
    for r in range(0,K-1):
        s=K-2-r
        sg=(-1)**s
        res+=Fr(sg,2)*chint(g-1,a+(r,s),rest)
        idx=list(range(n)); cidx=list(range(len(rest)))
        for p in range(n+1):
            for I in combinations(idx,p):
                J=[i for i in idx if i not in I]
                aI=tuple(a[i] for i in I); aJ=tuple(a[i] for i in J)
                for q in range(len(rest)+1):
                    for CI in combinations(cidx,q):
                        CJ=[i for i in cidx if i not in CI]
                        cI=tuple(rest[i] for i in CI); cJ=tuple(rest[i] for i in CJ)
                        for g1 in range(0,g+1):
                            g2=g-g1
                            if not stable(g1,p+1) or not stable(g2,n-p+1): continue
                            if sum(aI)+r+sum(cI)!=3*g1-3+p+1: continue
                            res+=Fr(sg,2)*chint(g1,aI+(r,),cI)*chint(g2,aJ+(s,),cJ)
    return pref*res

def odd_partitions(D,maxp=None):
    if D==0:
        yield (); return
    if maxp is None: maxp=D if D%2==1 else D-1
    p=maxp
    while p>=1:
        if p<=D:
            for q in odd_partitions(D-p,p):
                yield (p,)+q
        p-=2

def hodge_norm(g,a,wts):
    """int prod_w [Lambda^vee(w)/w^g] psi^a = int exp(-sum (p-1)! p_p(1/w) ch_p) psi^a"""
    n=len(a); D=3*g-3+n-sum(a)
    if not stable(g,n) or D<0: return Fr(0)
    if g==0: return wk(0,tuple(sorted(a))) if D==0 else Fr(0)
    tot=Fr(0)
    for part in odd_partitions(D):
        c=Fr(1)
        cnt=Counter(part)
        for p,mu in cnt.items():
            P=sum(Fr(1)/w**p for w in wts)
            c*=(-factorial(p-1)*P)**mu/factorial(mu)
        tot+=c*chint(g,tuple(sorted(a)),tuple(sorted(part)))
    return tot

def lam_int_control():
    u,v,w=Fr(2),Fr(3),Fr(7)
    e3=1/(u*v*w)
    l111=Fr(1,2880); l12=Fr(1,5760)
    S=0
    for a in (u,v,w):
        for b in (u,v,w):
            if a!=b: S+=Fr(1)/(a*b*b)
    exact=-e3*l111-S*l12
    got=hodge_norm(2,(),(u,v,w))
    return exact,got

def graphs(g,n):
    """stable graphs genus g, n<=1 legs."""
    out={}
    for V in range(1,2*g-2+n+1):
        pairs=[(u,v) for u in range(V) for v in range(u+1,V)]
        for gv in product(range(g+1),repeat=V):
            sg=sum(gv)
            if sg>g: continue
            for legs in product(range(n+1),repeat=V):
                if sum(legs)!=n: continue
                need=g-sg+V-1
                slots=V+len(pairs)
                def comps(k,t):
                    if k==1:
                        yield (t,); return
                    for i in range(t+1):
                        for c in comps(k-1,t-i):
                            yield (i,)+c
                if need<0: continue
                for c in comps(slots,need):
                    loops=c[:V]; ms=c[V:]
                    adj={p:m for p,m in zip(pairs,ms)}
                    seen={0}; st=[0]
                    while st:
                        x=st.pop()
                        for (u,v),m in adj.items():
                            if m>0:
                                if u==x and v not in seen: seen.add(v); st.append(v)
                                if v==x and u not in seen: seen.add(u); st.append(u)
                    if len(seen)<V: continue
                    val=[2*loops[v]+legs[v]+sum(m for (a,b),m in adj.items() if v in (a,b)) for v in range(V)]
                    if any(2*gv[v]-2+val[v]<=0 for v in range(V)): continue
                    best=None; auts=0
                    key0=(tuple((gv[i],loops[i],legs[i]) for i in range(V)),tuple(adj[(i,j)] for i in range(V) for j in range(i+1,V)))
                    for perm in permutations(range(V)):
                        kv=tuple((gv[perm[i]],loops[perm[i]],legs[perm[i]]) for i in range(V))
                        ke=tuple(adj[(min(perm[i],perm[j]),max(perm[i],perm[j]))] for i in range(V) for j in range(i+1,V))
                        k=(kv,ke)
                        if k==key0: auts+=1
                        if best is None or k<best: best=k
                    if best in out: continue
                    aut=auts
                    for v in range(V): aut*=2**loops[v]*factorial(loops[v])
                    for m in ms: aut*=factorial(m)
                    kv,ke=best
                    edges=[]
                    it=iter(ke)
                    for i in range(V):
                        for j in range(i+1,V):
                            m=next(it)
                            edges+= [(i,j)]*m
                    for i in range(V):
                        edges+=[(i,i)]*kv[i][1]
                    out[best]=dict(gv=[x[0] for x in kv],legs=[x[2] for x in kv],edges=edges,aut=aut)
    return list(out.values())

if __name__=="__main__":
    print("wk checks",wk(2,(4,)),wk(3,(7,)),wk(2,(2,3)),wk(2,(2,2,2)),wk(1,(1,)),wk(4,(10,)))
    print("hodge control g=2",lam_int_control())
    print("lambda_1 on M1,1 via ch:",chint(1,(0,),(1,)))
    print("M2 ch1^3 -> lambda1^3=1/2880:",chint(2,(),(1,1,1)))
    print("M3 lambda1^6=1/90720:",chint(3,(),(1,1,1,1,1,1)))
    print("(2,2) lam1 psi1^2psi2^2 =5/576:",chint(2,(2,2),(1,)))
    print("#graphs (7,42,?,2,?)",len(graphs(2,0)),len(graphs(3,0)),len(graphs(2,1)),len(graphs(1,1)),len(graphs(3,1)))
    print("sum aut g=2",[g['aut'] for g in graphs(2,0)])
    # lambda_g lambda_{g-1} lambda_{g-2} (FP): g=3 1/1451520.
    # (2g-1)! ch_{2g-1}=(-1)^{g-1} lam_g lam_{g-1}; lam_{g-2}=lam_1 for g=3
    print("g=3: 5! ch5 ch1 * (+1) expected 1/1451520:",factorial(5)*chint(3,(),(1,5)))
