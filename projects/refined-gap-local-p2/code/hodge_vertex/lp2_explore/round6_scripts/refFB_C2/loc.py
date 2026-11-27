# Direct virtual localisation for GW^T_{g,d} of Z = K_P2 x C^2 at generic alpha, eps (Graber-Pandharipande).
# Weights at P_i: {3a_i, a_j-a_i, a_k-a_i, -e1, -e2}. Evidence only. No mirror symmetry, no R-matrix.
import itertools, math, sys
from fractions import Fraction as Fr
from hod import Hodge, tau

class Loc:
    def __init__(s,alpha,e1,e2):
        s.a=[Fr(x) for x in alpha]; s.e1=Fr(e1); s.e2=Fr(e2)
        s.W=[]
        for i in range(3):
            j,k=[x for x in range(3) if x!=i]
            s.W.append([3*s.a[i],s.a[j]-s.a[i],s.a[k]-s.a[i],-s.e1,-s.e2])
        s.eT=[]
        for w in s.W:
            p=Fr(1)
            for x in w: p*=x
            s.eT.append(p)
        s.H=[Hodge(w) for w in s.W]
        s.vm={}
    def vertex(s,i,g,oms):
        """Int_{Mbar_{g,n}} prod Lambda^vee(w_a) prod 1/(om-psi) * eT^{n-1}, with unstable conventions"""
        n=len(oms)
        key=(i,g,tuple(sorted(oms)))
        if key in s.vm: return s.vm[key]
        if g==0 and n==1: r=oms[0]
        elif g==0 and n==2: r=s.eT[i]/(oms[0]+oms[1])
        else:
            dim=3*g-3+n
            tot=Fr(0)
            pre=s.eT[i]**g
            # 1/(om-psi)=Sum psi^k/om^{k+1}
            for ks in itertools.product(range(dim+1),repeat=n):
                if sum(ks)>dim: continue
                c=Fr(1)
                for o,k in zip(oms,ks): c/=o**(k+1)
                tot+=c*s.H[i](g,ks)
            r=tot*pre*s.eT[i]**(n-1)
        s.vm[key]=r
        return r
    def edge(s,i,j,d):
        """degree d cover of the line P_iP_j: e(H^1)/e(H^0) moving, times 1/d."""
        k=[x for x in range(3) if x not in (i,j)][0]
        a=s.a
        om=(a[j]-a[i])/d
        r=Fr(1,d)
        for m in range(1,d+1):
            r/=(m*om)*(-m*om)
        A=a[k]-a[i]; Bb=a[k]-a[j]
        for m in range(d+1):
            r/=(A-m*(A-Bb)/d)
        A3=3*a[i]; B3=3*a[j]
        for m in range(1,3*d):
            r*=A3+m*(B3-A3)/(3*d)
        r/=(s.e1*s.e2)
        return r,om
    def gw(s,g,d):
        if d==0:
            if g<2: raise ValueError
            return sum(s.eT[i]**(g-1)*s.H[i](g,()) for i in range(3))
        tot=Fr(0)
        # enumerate graphs: nv vertices with colours and genera, edges list (u,v,deg) u<v; brute-force with automorphism division
        for nv in range(2,d+2):
            for edges in edge_multisets(nv,d):
                h1=len(edges)-nv+1
                if h1<0 or h1>g: continue
                if not connected(nv,edges): continue
                gr=g-h1
                for gens in compositions(gr,nv):
                    for cols in itertools.product(range(3),repeat=nv):
                        if any(cols[u]==cols[v] for u,v,_ in edges): continue
                        # count this labelled configuration with weight 1/N, N = number of labelled configs in iso class * aut -> use nv!: orbit-stabiliser on labelled vertices
                        val=s.config(cols,gens,edges)
                        tot+=val*edge_mult_factor(edges)/math.factorial(nv)
        return tot
    def config(s,cols,gens,edges):
        nv=len(cols)
        oms=[[] for _ in range(nv)]
        val=Fr(1)
        for u,v,dg in edges:
            r,om=s.edge(cols[u],cols[v],dg)
            val*=r
            oms[u].append(om); oms[v].append(-om)
        for x in range(nv):
            val*=s.vertex(cols[x],gens[x],oms[x])
        return val

def edge_multisets(nv,d):
    """all multisets of edges (u<v,deg) with total degree d on labelled vertices"""
    pairs=[(u,v) for u in range(nv) for v in range(u+1,nv)]
    def rec(idx,rem):
        if rem==0:
            yield (); return
        if idx==len(pairs): return
        # choose multiset of degrees on pair idx
        for part in partitions_upto(rem):
            for r in rec(idx+1,rem-sum(part)):
                yield tuple((pairs[idx][0],pairs[idx][1],q) for q in part)+r
    return list(rec(0,d))
def partitions_upto(rem):
    # all partitions (non-increasing) of integers 0..rem
    def parts(n,mx):
        if n==0: yield (); return
        for p in range(min(n,mx),0,-1):
            for r in parts(n-p,p): yield (p,)+r
    for n in range(rem+1):
        yield from parts(n,n)
def connected(nv,edges):
    par=list(range(nv))
    def f(x):
        while par[x]!=x: x=par[x]
        return x
    for u,v,_ in edges: par[f(u)]=f(v)
    return len({f(i) for i in range(nv)})==1
def compositions(n,k):
    if k==1: yield (n,); return
    for i in range(n+1):
        for r in compositions(n-i,k-1): yield (i,)+r
def edge_mult_factor(edges):
    # identical parallel edges of equal degree: 1/m!
    c={}
    for e in edges: c[e]=c.get(e,0)+1
    r=Fr(1)
    for m in c.values(): r/=math.factorial(m)
    return r

def TF_limit(g,d,e1,e2,c=(1,2,-3),extra=2):
    """eps1 eps2 GW (or 3 e3 eps1 eps2 GW for d=0) is polynomial of degree <= 2g (2g+3) in alpha; interpolate along alpha=E/3+s c."""
    E=Fr(e1)+Fr(e2)
    def val(sv):
        a=[E/3+sv*x for x in c]
        L=Loc(a,e1,e2)
        v=L.gw(g,d)*Fr(e1)*Fr(e2)
        if d==0: v*=3*a[0]*a[1]*a[2]
        return v
    deg=2*g+(3 if d==0 else 0)
    pts=[Fr(k,7) for k in range(1,deg+2+extra)]
    vals=[val(p) for p in pts]
    # Newton divided differences; check top differences vanish
    n=len(pts); co=list(vals)
    for j in range(1,n):
        for i in range(n-1,j-1,-1):
            co[i]=(co[i]-co[i-1])/(pts[i]-pts[i-j])
    ok=all(x==0 for x in co[deg+1:])
    # evaluate at 0
    r=Fr(0)
    for i in range(deg,-1,-1):
        r=r*(0-pts[i])+co[i]
    r/=Fr(e1)*Fr(e2)
    if d==0: r/=3*(E/3)**3
    return r,ok

if __name__=="__main__":
    # genus-0 controls at generic alpha: GW_{0,d} eps1eps2 is constant: 3, -45/8, 244/9
    al=(Fr(2,3),Fr(5,7),Fr(-3,11))
    L=Loc(al,3,7)
    print("g0: eps1eps2 GW_{0,1},GW_{0,2},GW_{0,3}:",[L.gw(0,d)*21 for d in (1,2,3)],"expect 3,-45/8,244/9")
    # genus-1 d=1: expect refined: GW_{1,1}
    print("g1 d1 at generic alpha:",L.gw(1,1))
    # unrefined: e1=-e2 -> N_{g,1}: GV: 3*(2 sin(lam/2))^-2: N_{1,1}=1/4, N_{2,1}=1/80, N_{3,1}=1/2016
    for (g,d) in ((1,1),(2,1),(2,2),(3,1),(1,2)):
        r,ok=TF_limit(g,d,3,7)
        print("TF limit g,d",g,d,r,"poly-ok",ok)
