# R.omega graph sum (PPZ Defs 2.2, 2.13) with legs, in the Delta^(1/2)-normalised form P_a^{(i)}=Delta_i^(1/2) R_{ai}. Fresh code.
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from itertools import product,combinations_with_replacement,permutations
import math
from intersect import wk,ch1int
from gp import compositions

def ymul(a,b,D):
    r=[Fr(0)]*(D+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>D: break
            if y: r[i+j]+=x*y
    return r
def yadd(a,b): return [x+y for x,y in zip(a,b)]
def yscal(a,c): return [x*c for x in a]
def yinv(a,D):
    r=[Fr(0)]*(D+1);r[0]=1/a[0]
    for n in range(1,D+1):
        s=Fr(0)
        for k in range(1,n+1): s+=a[k]*r[n-k]
        r[n]=-s*r[0]
    return r
def ypow(a,n,D):
    if n<0: return ypow(yinv(a,D),-n,D)
    r=[Fr(1)]+[Fr(0)]*D
    for _ in range(n): r=ymul(r,a,D)
    return r

def stable_graphs(g,n):
    """list of (gs,edges,legs,1/|Aut|), edges=tuple of (u,v), u<=v. Aut in the half-edge sense."""
    res={}
    for V in range(1,max(2*g-2+n,1)+1):
        for b1 in range(0,g+1):
            E=V-1+b1
            pairs=[(u,v) for u in range(V) for v in range(u,V)]
            for ch in combinations_with_replacement(pairs,E):
                # connected
                seen={0};st=[0]
                while st:
                    x=st.pop()
                    for u,v in ch:
                        if u==x and v not in seen: seen.add(v);st.append(v)
                        if v==x and u not in seen: seen.add(u);st.append(u)
                if len(seen)!=V: continue
                for gs in compositions(g-b1,V):
                    for legs in product(range(V),repeat=n):
                        val=[0]*V
                        for u,v in ch: val[u]+=1;val[v]+=1
                        for l in legs: val[l]+=1
                        if any(2*gs[i]-2+val[i]<=0 for i in range(V)): continue
                        best=None;nauts=0
                        for p in permutations(range(V)):
                            gs2=[0]*V
                            for i in range(V): gs2[p[i]]=gs[i]
                            ee=tuple(sorted((min(p[u],p[v]),max(p[u],p[v])) for u,v in ch))
                            ll=tuple(p[x] for x in legs)
                            key=(tuple(gs2),ee,ll)
                            if best is None or key<best: best=key
                            if key==(tuple(gs),tuple(sorted(ch)),tuple(legs)): nauts+=1
                        # multiplicities
                        if best in res: continue
                        aut=nauts
                        cnt={}
                        for e in ch: cnt[e]=cnt.get(e,0)+1
                        for e,m in cnt.items():
                            aut*=math.factorial(m)
                            if e[0]==e[1]: aut*=2**m
                        res[best]=Fr(1,aut)
    return [(k[0],k[1],k[2],w) for k,w in res.items()]

class Rows:
    """ncol colours; P[a][j]: list over z-order of y-lists; Delta[j]: y-list; etai: matrix of Fractions."""
    def __init__(s,ncol,P,Delta,etai,D,K):
        s.ncol=ncol;s.P=P;s.Delta=Delta;s.etai=etai;s.D=D;s.K=K
        s.na=len(P)
        s.vc={};s.ec={}
    def z(s,a,j,n):
        L=s.P[a][j]
        return L[n] if n<len(L) else [Fr(0)]*(s.D+1)
    def dil(s,i,d):
        # Dil_i(d)=(-1)^d [z^{d-1}] P_0^{(i)}, d>=2
        v=s.z(0,i,d-1)
        return yscal(v,(-1)**d)
    def nMb(s,i,j,M,b):
        r=[Fr(0)]*(s.D+1)
        for c in range(s.na):
            for dd in range(s.na):
                e=s.etai[c][dd]
                if e==0: continue
                r=yadd(r,yscal(ymul(s.z(c,i,M),s.z(dd,j,b),s.D),e))
        return yscal(r,(-1)**(M+b))
    def E(s,i,j,a,b):
        key=(i,j,a,b)
        if key in s.ec: return s.ec[key]
        assert a+b+1<=s.K,'z-order too low for edge'
        r=[Fr(0)]*(s.D+1)
        for k in range(0,b+1):
            r=yadd(r,yscal(s.nMb(i,j,a+k+1,b-k),-(-1)**k))
        s.ec[key]=r
        return r
    def V(s,i,g,exps,rv=0):
        """Delta_i^{g-1} Sum_m (1/m!) Sum_{d_k>=2} <tau_exps tau_d>_g prod Dil_i(d_k)"""
        exps=tuple(sorted(exps))
        key=(i,g,exps,rv)
        if key in s.vc: return s.vc[key]
        n=len(exps)
        dim=3*g-3+n
        ex=dim-sum(exps)-rv
        r=[Fr(0)]*(s.D+1)
        if ex>=0:
            for part in partitions(ex):
                # parts p=d-1>=1
                ds=tuple(sorted(p+1 for p in part))
                w=ch1int(g,tuple(sorted(exps+ds)),rv)
                if w==0: continue
                mult={}
                for x in ds: mult[x]=mult.get(x,0)+1
                coef=Fr(1)
                for x,m in mult.items(): coef/=math.factorial(m)
                term=[Fr(1)]+[Fr(0)]*s.D
                for x in ds: term=ymul(term,s.dil(i,x),s.D)
                r=yadd(r,yscal(term,w*coef))
        r=ymul(r,ypow(s.Delta[i],g-1,s.D),s.D)
        s.vc[key]=r
        return r
    def integral(s,g,ins,legsign=-1,lam=0):
        """Int_{Mbar_{g,n}} prod psi_l^{k_l} (R.omega)_{g,n}(H^{a_1},...), ins=[(k,a)]. Leg series P_a(legsign psi).  y-list."""
        n=len(ins)
        tot=[Fr(0)]*(s.D+1)
        for gs,edges,legs,w in graphs(g,n):
            if lam==0:
                tot=yadd(tot,yscal(s.graph_val(gs,edges,legs,ins,legsign),w))
            else:
                pos=[v for v in range(len(gs)) if gs[v]>0]
                if not pos: continue
                for dist in compositions(lam,len(pos)):
                    rv=[0]*len(gs);mult=math.factorial(lam)
                    for v,r in zip(pos,dist):
                        rv[v]=r;mult//=math.factorial(r)
                    tot=yadd(tot,yscal(s.graph_val(gs,edges,legs,ins,legsign,rv),w*mult))
        return tot
    def graph_val(s,gs,edges,legs,ins,legsign,rv=None):
        V=len(gs);E=len(edges);n=len(legs)
        val=[0]*V
        for u,v in edges: val[u]+=1;val[v]+=1
        for l in legs: val[l]+=1
        if rv is None: rv=[0]*V
        dims=[3*gs[i]-3+val[i]-rv[i] for i in range(V)]
        if any(x<0 for x in dims): return [Fr(0)]*(s.D+1)
        tot=[Fr(0)]*(s.D+1)
        # half-edge variable list: edges (two) and legs (one)
        hes=[]
        for ei,(u,v) in enumerate(edges):
            hes.append(u);hes.append(v)
        for l,x in enumerate(legs): hes.append(x)
        base=[0]*V
        for l,x in enumerate(legs): base[x]+=ins[l][0]
        if any(base[i]>dims[i] for i in range(V)):
            return tot
        def rec(idx,used,exps):
            if idx==len(hes):
                yield tuple(exps);return
            x=hes[idx]
            for e in range(0,dims[x]-used[x]+1):
                used[x]+=e;exps.append(e)
                yield from rec(idx+1,used,exps)
                exps.pop();used[x]-=e
        for cols in product(range(s.ncol),repeat=V):
            for exps in rec(0,list(base),[]):
                term=[Fr(1)]+[Fr(0)]*s.D
                zero=False
                for ei,(u,v) in enumerate(edges):
                    f=s.E(cols[u],cols[v],exps[2*ei],exps[2*ei+1])
                    if not any(f): zero=True;break
                    term=ymul(term,f,s.D)
                if zero: continue
                for l,x in enumerate(legs):
                    m=exps[2*E+l]
                    f=s.z(ins[l][1],cols[x],m)
                    if legsign==-1 and m%2: f=yscal(f,-1)
                    if not any(f): zero=True;break
                    term=ymul(term,f,s.D)
                if zero: continue
                vex=[[] for _ in range(V)]
                for ei,(u,v) in enumerate(edges):
                    vex[u].append(exps[2*ei]);vex[v].append(exps[2*ei+1])
                for l,x in enumerate(legs):
                    vex[x].append(exps[2*E+l]+ins[l][0])
                for i in range(V):
                    f=s.V(cols[i],gs[i],tuple(vex[i]),rv[i])
                    if not any(f): zero=True;break
                    term=ymul(term,f,s.D)
                if zero: continue
                tot=yadd(tot,term)
        return tot

def partitions(n,m=None):
    if m is None: m=n
    if n==0:
        yield ();return
    for k in range(min(n,m),0,-1):
        for r in partitions(n-k,k):
            yield (k,)+r

_gcache={}
def graphs(g,n):
    if (g,n) not in _gcache: _gcache[(g,n)]=stable_graphs(g,n)
    return _gcache[(g,n)]
