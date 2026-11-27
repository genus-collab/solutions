# Referee's own code (evidence only): T_F point-frame colour sum at genus g (2,3) with rows from the closed form,
# colours in Q(omega), stable graphs enumerated here. Compare with Hodge constant term + refined BPS at tiny y.
import sys,math,time,itertools,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from functools import lru_cache
from wick_ref import Theta,mul,S,rec
sys.setrecursionlimit(10000)
# ---- Q(omega)
class W:
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=Fr(a);s.b=Fr(b)
    def __add__(s,o):
        o=o if isinstance(o,W) else W(o); return W(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __neg__(s): return W(-s.a,-s.b)
    def __sub__(s,o):
        o=o if isinstance(o,W) else W(o); return W(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return W(o)-s
    def __mul__(s,o):
        if not isinstance(o,W): return W(s.a*o,s.b*o)
        # w^2=-1-w
        return W(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a-s.b*o.b)
    __rmul__=__mul__
    def inv(s):
        n=s.a*s.a-s.a*s.b+s.b*s.b
        return W((s.a-s.b)/n,-s.b/n)
    def __truediv__(s,o):
        o=o if isinstance(o,W) else W(o); return s*o.inv()
    def __rtruediv__(s,o): return W(o)*s.inv()
    def __pow__(s,n):
        if n<0: return s.inv()**(-n)
        r=W(1)
        for _ in range(n): r=r*s
        return r
    def iszero(s): return s.a==0 and s.b==0
    def __repr__(s): return "(%s+%sw)"%(s.a,s.b)
om=W(0,1)
# ---- WK via DVV
def df(n):
    r=1
    while n>1: r*=n;n-=2
    return r
@lru_cache(None)
def wk(g,ds):
    ds=tuple(sorted(ds)); n=len(ds)
    if g<0 or 2*g-2+n<=0: return Fr(0)
    if sum(ds)!=3*g-3+n: return Fr(0)
    if g==0 and n==3: return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    if ds[0]==0:
        rest=ds[1:];tot=Fr(0)
        for i,d in enumerate(rest):
            if d>0:
                nw=list(rest);nw[i]=d-1;tot+=wk(g,tuple(sorted(nw)))
        return tot
    # DVV on last index k+1=ds[-1]
    k=ds[-1]-1; rest=ds[:-1]
    tot=Fr(0)
    for j,d in enumerate(rest):
        nw=list(rest);nw[j]=d+k
        tot+=Fr(df(2*k+2*d+1),df(2*d-1))*wk(g,tuple(sorted(nw)))
    for r in range(0,k):
        s=k-1-r
        co=Fr(df(2*r+1)*df(2*s+1),2)
        tot+=co*wk(g-1,tuple(sorted((r,s)+rest)))
        m=len(rest)
        for mask in range(1<<m):
            I=tuple(rest[i] for i in range(m) if mask>>i&1)
            J=tuple(rest[i] for i in range(m) if not mask>>i&1)
            di=sum(I)+r; ni=len(I)+1
            if (di+3-ni)%3: continue
            g1=(di+3-ni)//3
            g2=g-g1
            if g1<0 or g2<0: continue
            tot+=co*wk(g1,tuple(sorted((r,)+I)))*wk(g2,tuple(sorted((s,)+J)))
    return tot/df(2*k+3)
assert wk(2,(4,))==Fr(1,1152) and wk(3,(7,))==Fr(1,82944) and wk(2,(2,3))==Fr(29,5760) and wk(2,(2,2,2))==Fr(7,240)
# ---- stable graphs
def graphs(g):
    res=[]
    for nv in range(1,2*g-1):
        for gs in itertools.product(range(g+1),repeat=nv):
            if list(gs)!=sorted(gs): continue
            sg=sum(gs)
            if sg>g: continue
            h1=g-sg
            ne=h1+nv-1
            pairs=[(u,v) for u in range(nv) for v in range(u,nv)]
            for edges in itertools.combinations_with_replacement(pairs,ne):
                par=list(range(nv))
                def fd(x):
                    while par[x]!=x: x=par[x]
                    return x
                for u,v in edges: par[fd(u)]=fd(v)
                if len(set(fd(i) for i in range(nv)))>1: continue
                val=[0]*nv
                for u,v in edges: val[u]+=1;val[v]+=1
                if any(2*gs[i]-2+val[i]<=0 for i in range(nv)): continue
                best=None
                auts=0
                for p in itertools.permutations(range(nv)):
                    gp=tuple(gs[p[i]] for i in range(nv))
                    inv=[0]*nv
                    for i in range(nv): inv[p[i]]=i
                    ep=tuple(sorted(tuple(sorted((inv[u],inv[v]))) for u,v in edges))
                    key=(gp,ep)
                    if best is None or key<best: best=key
                    if gp==tuple(gs) and ep==tuple(sorted(edges)): auts+=1
                res.append((best,auts,gs,edges))
    seen={}
    for best,auts,gs,edges in res:
        if best in seen: continue
        cnt={}
        for e in edges: cnt[e]=cnt.get(e,0)+1
        a=auts
        for (u,v),m in cnt.items():
            a*=math.factorial(m)
            if u==v: a*=2**m
        seen[best]=(gs,edges,a)
    return list(seen.values())
# ---- rows as Laurent polys in v: dict power->Fraction, for each z-order
class Rows:
    def __init__(s,e1,e2,Z):
        s.e1,s.e2,s.E,s.Z=e1,e2,e1+e2,Z
        E=s.E
        r=rec(Z+2)
        S0=[{} for _ in range(Z+1)]
        for l in range(Z+1):
            t=[Fr(0)]*(Z+1);t[0]=Fr(1)
            for j in range(1,l+2):
                c=Fr(j)/E
                u=[Fr(0)]*(Z+1);u[0]=t[0]
                for n in range(1,Z+1): u[n]=t[n]-c*u[n-1]
                t=u
            for n in range(Z+1-l):
                S0[n+l][l]=S0[n+l].get(l,0)+r[l]*t[n]
        h=mul(mul(Theta(-1/E,Z),Theta(1/e1,Z),Z),Theta(1/e2,Z),Z)
        P0=s.zmulc(h,S0)
        s.P0=P0
        s.Q=s.Top(P0)
        s.TQ=s.Top(s.Q)
    def zmulc(s,h,A):
        Z=s.Z;out=[{} for _ in range(Z+1)]
        for n in range(Z+1):
            for k in range(n+1):
                for p,c in A[n-k].items():
                    out[n][p]=out[n].get(p,0)+h[k]*c
        return out
    def Top(s,g):
        # Tbar g = (-E/3+1/v) g + z (1/E - v/3)(g' + g/v)
        E=s.E;Z=s.Z;out=[{} for _ in range(Z+1)]
        def add(n,p,c):
            if n<=Z: out[n][p]=out[n].get(p,0)+c
        for n in range(Z+1):
            for p,c in g[n].items():
                add(n,p,-E/3*c); add(n,p-1,c)
                co=(p+1)*c   # coefficient of v^{p-1} in g'+g/v
                add(n+1,p-1,co/E); add(n+1,p,-co/3)
        return out
def ev(d,vp):
    t=W(0)
    for p,c in d.items():
        if c!=0: t=t+vp(p)*c
    return t
class Total:
    def __init__(s,e1,e2,t,Fv,Z):
        s.R=Rows(e1,e2,Z); s.e1,s.e2,s.E=e1,e2,e1+e2; s.k=e1*e2; s.t=Fr(t); s.y=s.t**3; s.Dc=1+27*s.y; s.F=Fv;s.Z=Z
        E=s.E
        ti=[W(s.t),om*s.t,om*om*s.t]
        s.v=[(1+3*x)/(x*E) for x in ti]
        s.vinv=[x.inv() for x in s.v]
        s.pw=[{} for _ in range(3)]
        s.Del=[s.vp(i,-2)*(3*s.k*E) for i in range(3)]
        s.P0=[[ev(s.R.P0[n],lambda p:s.vp(i,p)) for n in range(Z+1)] for i in range(3)]
        s.Q=[[ev(s.R.Q[n],lambda p:s.vp(i,p)) for n in range(Z+1)] for i in range(3)]
        s.TQ=[[ev(s.R.TQ[n],lambda p:s.vp(i,p)) for n in range(Z+1)] for i in range(3)]
        s.eta00=s.k*E**3/9
        s._E={}
        s._V={}
    def vp(s,i,p):
        d=s.pw[i]
        if p in d: return d[p]
        r=(s.v[i]**p) if p>=0 else (s.vinv[i]**(-p))
        d[p]=r
        return r
    def pair(s,i,j,M,b):
        k=s.k;E=s.E;Dc=s.Dc;F=s.F
        T=s.P0[i][M]*s.P0[j][b]*s.eta00
        T=T-s.Q[i][M]*s.Q[j][b]*(3*E*k)
        T=T-(s.Q[i][M]*s.TQ[j][b]+s.TQ[i][M]*s.Q[j][b])*(3*k*Dc)
        if M>=1: T=T+s.Q[i][M-1]*s.Q[j][b]*(3*k*Dc*F)
        if b>=1: T=T+s.Q[i][M]*s.Q[j][b-1]*(3*k*Dc*F)
        return T
    def edge(s,i,j,a,b):
        key=(i,j,a,b)
        if key in s._E: return s._E[key]
        T=W(0)
        for kk in range(b+1):
            M=a+kk+1;bb=b-kk
            T=T-s.pair(i,j,M,bb)*Fr((-1)**(kk+M+bb))
        s._E[key]=T
        return T
    def dil(s,i,d): return s.P0[i][d-1]*Fr((-1)**d)
    def vert(s,i,g,hs):
        hs=tuple(sorted(hs))
        key=(i,g,hs)
        if key in s._V: return s._V[key]
        n=len(hs); dim=3*g-3+n-sum(hs)
        tot=W(0)
        if dim>=0:
            def parts(rem,mx):
                if rem==0:
                    yield ();return
                for p in range(min(rem,mx),0,-1):
                    for r in parts(rem-p,p): yield (p,)+r
            for pt in parts(dim,dim):
                ds=tuple(p+1 for p in pt)
                val=wk(g,tuple(sorted(hs+ds)))
                if val==0: continue
                mult=1
                for x in set(ds): mult*=math.factorial(ds.count(x))
                term=W(val/mult)
                for d in ds: term=term*s.dil(i,d)
                tot=tot+term
        if g==0: tot=tot*s.Del[i].inv()
        elif g>1: tot=tot*(s.Del[i]**(g-1))
        s._V[key]=tot
        return tot
    def graph(s,gs,edges):
        nv=len(gs)
        tot=W(0)
        he=[]
        for (u,v) in edges:
            he.append(u);he.append(v)
        nh=len(he)
        val=[he.count(i) for i in range(nv)]
        dims=[3*gs[i]-3+val[i] for i in range(nv)]
        for col in itertools.product(range(3),repeat=nv):
            def rec_(h,exps,used):
                if h==nh:
                    prod=W(1)
                    for e_idx in range(len(edges)):
                        u,v=edges[e_idx]
                        ee=s.edge(col[u],col[v],exps[2*e_idx],exps[2*e_idx+1])
                        prod=prod*ee
                        if prod.iszero(): return W(0)
                    for i in range(nv):
                        hs=[exps[q] for q in range(nh) if he[q]==i]
                        prod=prod*s.vert(col[i],gs[i],tuple(hs))
                        if prod.iszero(): return W(0)
                    return prod
                res=W(0)
                vx=he[h]
                for a in range(dims[vx]-used[vx]+1):
                    used[vx]+=a
                    res=res+rec_(h+1,exps+[a],used)
                    used[vx]-=a
                return res
            tot=tot+rec_(0,[],[0]*nv)
        return tot
    def total(s,g,gl):
        tot=W(0)
        for gs,edges,a in gl:
            tot=tot+s.graph(gs,edges)*Fr(1,a)
        return tot
# ---- comparison side
def hodge_const(g,e1,e2):
    # GW^{T_F}_{g,0}=Int_{Mbar_g}[9 l_g l_{g-2} P(0)-3 l_g l_{g-1} P'(0)], P(H)=Lv(E-3H)Lv(-e1)Lv(-e2)/((E-3H)e1 e2)
    E=e1+e2
    def Lv(u):
        d={}
        for i in range(g+1):
            e=[0]*g
            if i>0: e[i-1]=1
            d[tuple(e)]=Fr((-1)**i)*u**(g-i)
        return d
    def pm(a,b):
        r={}
        for m,c in a.items():
            for n,x in b.items():
                k=tuple(p+q for p,q in zip(m,n))
                r[k]=r.get(k,0)+c*x
        return r
    C=pm(Lv(-e1),Lv(-e2))
    P0=pm(C,dict((m,c/(E*e1*e2)) for m,c in Lv(E).items()))
    dL={}
    for i in range(g+1):
        e=[0]*g
        if i>0: e[i-1]=1
        p=g-i-1  # Lv(u)/u = Sum (-1)^i l_i u^p ; d/du: p u^(p-1)
        dL[tuple(e)]=Fr((-1)**i)*p*(E**(p-1) if p!=0 else 0)
    P1=pm(C,dict((m,-3*c/(e1*e2)) for m,c in dL.items()))
    def li(i):
        e=[0]*g
        if i>0: e[i-1]=1
        return {tuple(e):Fr(1)}
    A=pm(pm(li(g),li(g-2)),P0)
    Bb=pm(pm(li(g),li(g-1)),P1)
    def integ(d):
        tot=Fr(0)
        for m,c in d.items():
            deg=sum((i+1)*x for i,x in enumerate(m))
            if deg!=3*g-3 or c==0: continue
            tot+=c*mono(g,m)
        return tot
    return 9*integ(A)-3*integ(Bb)
def mono(g,m):
    if g==2:
        l1,l2=m
        if l2>=2: return Fr(0)
        if l2==1 and l1==1: return Fr(1,5760)
        if l2==0 and l1==3: return Fr(1,2880)
        return Fr(0)
    if g==3:
        l1,l2,l3=m
        val=Fr(1,1451520) # l1 l2 l3
        c=Fr(1)
        while True:
            if l3>=2: return Fr(0)
            if l1>=2:
                l1-=2;l2+=1;c*=2;continue
            if l2>=2:
                l2-=2;l1+=1;l3+=1;c*=2;continue
            break
        return c*val if (l1,l2,l3)==(1,1,1) else Fr(0)
def invsinh_series(a,b,N):
    def inv(c):
        ser=[Fr(0)]*(N+1)
        for n in range(0,N//2+1):
            ser[2*n]=(Fr(c)/2)**(2*n)/math.factorial(2*n+1)
        out=[Fr(0)]*(N+1);out[0]=Fr(1)
        for n in range(1,N+1):
            out[n]=-sum(ser[k]*out[n-k] for k in range(1,n+1))
        return out
    return mul(inv(a),inv(b),N)
def refined(g,e1,e2,dmax):
    # degree 2g-2 part of the d-th coefficient of the refined BPS sum
    E=e1+e2
    N=2*g
    BPS={1:[(0,1,1)],2:[(0,Fr(5,2),1)],3:[(0,3,1),(Fr(1,2),Fr(9,2),1)]}
    res={}
    for d in range(1,dmax+1):
        tot=Fr(0)
        for k in range(1,d+1):
            if d%k: continue
            for jl,jr,mult in BPS[d//k]:
                sign=(-1)**int(2*jl+2*jr)
                ex=[]
                for i in range(int(2*jl)+1):
                    for j in range(int(2*jr)+1):
                        ml=Fr(jl)-i; mr=Fr(jr)-j
                        ex.append(k*(ml*(e1-e2)+mr*E))
                ch=[sum(x**n for x in ex)/math.factorial(n) for n in range(N+1)]
                ins=invsinh_series(k*e1,k*e2,N)
                ser=mul(ch,ins,N)
                tot+=sign*mult*ser[2*g]/(k*k*k*e1*e2)
        res[d]=tot
    return res
def Aseries(N):
    return [Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d for d in range(N+1)]
def Fval(y,N=14):
    a=Aseries(N)
    A=sum(a[k]*y**k for k in range(N+1)); dA=sum(k*a[k]*y**k for k in range(N+1))
    return dA/A
def Qval(y,N=14):
    x=sum(3*Fr(math.factorial(3*d-1),math.factorial(d)**3)*(-y)**d for d in range(1,N+1))
    e=Fr(1);term=Fr(1)
    for n in range(1,N+1):
        term=term*x/n;e+=term
    return y*e
if __name__=="__main__":
    g=int(sys.argv[1])
    gl=graphs(g)
    print("genus",g,"stable graphs:",len(gl),flush=True)
    samples=[(Fr(3),Fr(7)),(Fr(1),Fr(1)),(Fr(5,2),Fr(-3,7))]
    Z=3*g-3+2
    for e1,e2 in samples:
        t=Fr(1,10**4)
        y=t**3
        t0=time.time()
        T=Total(e1,e2,t,Fval(y),Z)
        val=T.total(g,gl)
        c0=hodge_const(g,e1,e2)
        rf=refined(g,e1,e2,3)
        Q=Qval(y)
        print("eps",e1,e2,"imag part zero:",val.b==0,"time %.1f"%(time.time()-t0),flush=True)
        print("  total-GW0 =%.4e ; (total-GW0)/Q - GW1 = %.4e ; ((total-GW0)-Q GW1)/Q^2 - GW2 = %.4e ; ...-GW3 = %.4e"%(float(val.a-c0),float((val.a-c0)/Q-rf[1]),float(((val.a-c0)-Q*rf[1])/Q**2-rf[2]),float(((val.a-c0)-Q*rf[1]-Q*Q*rf[2])/Q**3-rf[3])),flush=True)
        print("  GW0",c0,"GW1",rf[1],"GW2",rf[2],flush=True)
