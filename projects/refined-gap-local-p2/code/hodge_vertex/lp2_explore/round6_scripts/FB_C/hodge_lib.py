# Evidence only. Exact Fractions. Witten-Kontsevich numbers (DVV), Hodge integrals by Mumford/Faber recursion,
# stable graphs without legs, and Q(omega) arithmetic.
import math, itertools
from fractions import Fraction as Fr
from functools import lru_cache

def bern(nmax):
    B=[Fr(0)]*(nmax+1); B[0]=Fr(1)
    for m in range(1,nmax+1):
        B[m]=-sum(Fr(math.comb(m+1,k))*B[k] for k in range(m))/(m+1)
    return B
BN=bern(40)

def dfact(n):
    r=1
    while n>1:
        r*=n; n-=2
    return r

def subsets(t):
    n=len(t)
    for mask in range(1<<n):
        A=tuple(t[i] for i in range(n) if mask>>i&1)
        B=tuple(t[i] for i in range(n) if not mask>>i&1)
        yield A,B

@lru_cache(None)
def wk(g,ds):
    ds=tuple(sorted(ds)); n=len(ds)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(ds)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    if ds[0]==0:
        rest=ds[1:]; tot=Fr(0)
        for i,d in enumerate(rest):
            if d>0:
                new=list(rest); new[i]=d-1
                tot+=wk(g,tuple(sorted(new)))
        return tot
    # DVV with k+1 = largest
    k=ds[-1]-1; rest=ds[:-1]
    tot=Fr(0)
    for j,d in enumerate(rest):
        new=list(rest); new[j]=d+k
        tot+=Fr(dfact(2*k+2*d+1),dfact(2*d-1))*wk(g,tuple(sorted(new)))
    for r in range(0,k):
        s=k-1-r
        c=Fr(dfact(2*r+1)*dfact(2*s+1),2)
        tot+=c*wk(g-1,tuple(sorted(rest+(r,s))))
        for A,B in subsets(rest):
            for g1 in range(0,g+1):
                g2=g-g1
                v1=wk(g1,tuple(sorted(A+(r,))))
                if v1==0: continue
                tot+=c*v1*wk(g2,tuple(sorted(B+(s,))))
    return tot/dfact(2*k+3)

@lru_cache(None)
def ch_int(g,a,K):
    """Int_{Mbar_{g,n}} prod ch_{K_j}(E) prod psi^{a_i}, K_j odd."""
    a=tuple(sorted(a)); K=tuple(sorted(K)); n=len(a)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(a)+sum(K)!=3*g-3+n: return Fr(0)
    if not K: return wk(g,a)
    if g==0: return Fr(0)
    k=K[-1]; Kr=K[:-1]
    l=(k+1)//2
    pref=BN[2*l]/math.factorial(2*l)
    tot=ch_int(g,tuple(sorted(a+(2*l,))),Kr)
    for i in range(n):
        new=list(a); new[i]+=2*l-1
        tot-=ch_int(g,tuple(sorted(new)),Kr)
    for x in range(0,2*l-1):
        y=2*l-2-x
        sg=(-1)**x
        t=ch_int(g-1,tuple(sorted(a+(x,y))),Kr)
        for A,B in subsets(a):
            for K1,K2 in subsets(Kr):
                for g1 in range(0,g+1):
                    g2=g-g1
                    v1=ch_int(g1,tuple(sorted(A+(x,))),tuple(sorted(K1)))
                    if v1==0: continue
                    t+=v1*ch_int(g2,tuple(sorted(B+(y,))),tuple(sorted(K2)))
        tot+=Fr(sg,2)*t
    return pref*tot

def odd_partitions(D,maxp=None):
    if maxp is None: maxp=D
    if D==0:
        yield (); return
    for p in range(min(D,maxp),0,-1):
        if p%2==0: continue
        for r in odd_partitions(D-p,p):
            yield (p,)+r

def hodge_cache(wts):
    """returns fn(g,a) = Int Omega_w psi^a, Omega_w = prod Lambda^vee(w)/w^g = exp(-Sum (2k-2)! ch_{2k-1} p_{2k-1}(1/w))."""
    memo={}
    def p(k): return sum(Fr(1)/w**k for w in wts)
    def fn(g,a):
        a=tuple(sorted(a)); n=len(a)
        key=(g,a)
        if key in memo: return memo[key]
        D=3*g-3+n-sum(a)
        tot=Fr(0)
        if D>=0:
            for part in odd_partitions(D):
                c=Fr(1); mult={}
                for q in part:
                    mult[q]=mult.get(q,0)+1
                    c*=-math.factorial(q-1)*p(q)
                for m in mult.values(): c/=math.factorial(m)
                tot+=c*ch_int(g,a,part)
        memo[key]=tot
        return tot
    return fn

# ---- stable graphs without legs
def stable_graphs(g):
    """list of (genera tuple, edges list [(u,v)] with u<=v, |Aut|)."""
    out={}
    for nv in range(1,2*g-2+1):
        for gen in itertools.product(range(g+1),repeat=nv):
            if sum(gen)>g: continue
            h1=g-sum(gen)
            ne=h1+nv-1
            pairs=[(i,j) for i in range(nv) for j in range(i,nv)]
            for ed in itertools.combinations_with_replacement(pairs,ne):
                # connected
                par=list(range(nv))
                def f(x):
                    while par[x]!=x: x=par[x]
                    return x
                for u,v in ed: par[f(u)]=f(v)
                if len({f(i) for i in range(nv)})!=1: continue
                val=[0]*nv
                for u,v in ed:
                    val[u]+=1; val[v]+=1
                if any(2*gen[i]-2+val[i]<=0 for i in range(nv)): continue
                # canonical form
                best=None; naut=0
                forms={}
                for perm in itertools.permutations(range(nv)):
                    g2=[0]*nv
                    for i in range(nv): g2[perm[i]]=gen[i]
                    e2=sorted(tuple(sorted((perm[u],perm[v]))) for u,v in ed)
                    fm=(tuple(g2),tuple(e2))
                    forms[fm]=forms.get(fm,0)+1
                    if best is None or fm<best: best=fm
                if best in out: continue
                naut=forms[best]
                gg,ee=best
                aut=naut
                cnt={}
                for e in ee: cnt[e]=cnt.get(e,0)+1
                for e,m in cnt.items():
                    aut*=math.factorial(m)
                    if e[0]==e[1]: aut*=2**m
                out[best]=aut
    return [(k[0],list(k[1]),v) for k,v in sorted(out.items(),key=lambda kv:(len(kv[0][0]),kv[0]))]

# ---- Q(omega), omega^2 = -1-omega
class W:
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o):
        if not isinstance(o,W): return W(s.a+o,s.b)
        return W(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __neg__(s): return W(-s.a,-s.b)
    def __sub__(s,o):
        if not isinstance(o,W): return W(s.a-o,s.b)
        return W(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return W(o-s.a,-s.b)
    def __mul__(s,o):
        if not isinstance(o,W): return W(s.a*o,s.b*o)
        # (a+bw)(c+dw)=ac+(ad+bc)w+bd w^2 = ac-bd + (ad+bc-bd) w
        return W(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a-s.b*o.b)
    __rmul__=__mul__
    def inv(s):
        # norm = a^2-ab+b^2; inverse conj/norm, conj(a+bw)=a+b w^2 = a-b - b w
        n=s.a*s.a-s.a*s.b+s.b*s.b
        return W((s.a-s.b)/n,-s.b/n)
    def __truediv__(s,o):
        if not isinstance(o,W): return W(s.a/o,s.b/o)
        return s*o.inv()
    def __rtruediv__(s,o): return s.inv()*o
    def __pow__(s,n):
        if n<0: return s.inv()**(-n)
        r=W(1,0); x=s
        while n:
            if n&1: r=r*x
            x=x*x; n>>=1
        return r
    def iszero(s): return s.a==0 and s.b==0
    def __repr__(s): return f"({s.a}+{s.b}w)"
OM=W(0,1)

if __name__=="__main__":
    print("graph counts",[len(stable_graphs(g)) for g in (2,3)])
    print("sum 1/Aut (g=2)",sum(Fr(1,x[2]) for x in stable_graphs(2)))
    print("tau4",wk(2,(4,)),"tau7",wk(3,(7,)),"tau2^3",wk(2,(2,2,2)),"(expect 1/1152,1/82944,7/240)")
    # lambda checks: lam1=ch1, on M2: lam1^3=ch1^3 -> 1/2880; M3: lam1^6 -> 1/90720
    print("lam1^3 M2",ch_int(2,(),(1,1,1)),"expect 1/2880")
    print("lam1^6 M3",ch_int(3,(),(1,)*6),"expect 1/90720")
    # lam1 lam2 lam3 via Hodge: coefficient
    print("ch3 lam-g check M2: 3! ch3 = lam1 lam2 sign:", 6*ch_int(2,(),(3,)), "expect (-1)^(g-1) lam_2 lam_1 = -1/5760")
