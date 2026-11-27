# Referee's own T_F graph sums. Rows from Def 1.3 with free c1 and free F (rows built, then pairing evaluated by eta^{ab}).
# Evidence only.
import itertools, math, sys, re
from fractions import Fraction as Fr
from collections import Counter
from hod import tau, Hodge, B

class Om:
    "a+b w, w^2=-1-w"
    __slots__=('a','b')
    def __init__(s,a=0,b=0):
        s.a=a if isinstance(a,Fr) else Fr(a); s.b=b if isinstance(b,Fr) else Fr(b)
    def _c(o): return o if isinstance(o,Om) else Om(o)
    def __add__(s,o):
        o=Om._c(o); return Om(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __sub__(s,o):
        o=Om._c(o); return Om(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return Om._c(o)-s
    def __neg__(s): return Om(-s.a,-s.b)
    def __mul__(s,o):
        o=Om._c(o)
        return Om(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a-s.b*o.b)
    __rmul__=__mul__
    def inv(s):
        n=s.a*s.a-s.a*s.b+s.b*s.b
        return Om((s.a-s.b)/n,-s.b/n)
    def __truediv__(s,o): return s*Om._c(o).inv()
    def __rtruediv__(s,o): return Om._c(o)*s.inv()
    def __pow__(s,n):
        if n<0: return s.inv()**(-n)
        r=Om(1)
        for _ in range(n): r=r*s
        return r
    def zero(s): return s.a==0 and s.b==0
    def __repr__(s): return f"({s.a}+{s.b}w)"
W=Om(0,1)

def rec_r(n):
    r=[Fr(1)]
    for m in range(n):
        rm1=r[m-1] if m>=1 else Fr(0)
        r.append((3*(3*m*m+3*m+1)*r[m]-m**3*rm1)/(27*(m+1)))
    return r

def theta_coeffs(c,N):
    """Theta(c z) = exp Sum_{m>=1} B_{2m} (c z)^{2m-1}/(2m(2m-1)) as list of z-coefficients"""
    lg=[Fr(0)]*(N+1)
    for m in range(1,N+2):
        if 2*m-1<=N: lg[2*m-1]=B[2*m]*c**(2*m-1)/(2*m*(2*m-1))
    out=[Fr(0)]*(N+1); out[0]=Fr(1)
    for n in range(1,N+1):
        out[n]=sum(k*lg[k]*out[n-k] for k in range(1,n+1))/n
    return out
def mul_ser(a,b,N):
    out=[Fr(0)]*(N+1)
    for i in range(N+1):
        if a[i]==0: continue
        for j in range(N+1-i):
            out[i+j]+=a[i]*b[j]
    return out

class Closed:
    """z-coefficients of S, Q=TS, TQ as Laurent polys in v, independent of colour."""
    def __init__(s,E,N):
        s.E=Fr(E); s.N=N; E=s.E
        r=rec_r(N+2)
        # S_k(v)=Sum_l r_l v^l [z^{k-l}] 1/prod_{j<=l+1}(1+jz/E)
        s.S=[{} for _ in range(N+1)]
        for l in range(N+1):
            ser=[Fr(0)]*(N+1); ser[0]=Fr(1)
            for j in range(1,l+2):
                for n in range(1,N+1): ser[n]-=Fr(j)/E*ser[n-1]
            for k in range(l,N+1):
                s.S[k][l]=s.S[k].get(l,Fr(0))+r[l]*ser[k-l]
        s.Q=s.T(s.S); s.TQ=s.T(s.Q)
    def T(s,f):
        E=s.E; N=s.N
        def theta(p):
            o={}
            for e,x in p.items():
                if e==0: continue
                o[e-1]=o.get(e-1,Fr(0))+x*e/E
                o[e]=o.get(e,Fr(0))-x*e/3
            return o
        def mv(p,k,c=Fr(1)):
            return {e+k:c*x for e,x in p.items()}
        def add(*ps):
            o={}
            for p in ps:
                for e,x in p.items(): o[e]=o.get(e,Fr(0))+x
            return o
        out=[]
        for n in range(N+1):
            g=add(mv(f[n],0,-E/3),mv(f[n],-1))   # Y=-E/3+1/v
            if n>=1:
                # z(theta + (theta v)/v) f_{n-1}; theta v/v = 1/(E v)-1/3
                g=add(g,theta(f[n-1]),mv(f[n-1],-1,1/E),mv(f[n-1],0,Fr(-1,3)))
            out.append(g)
        return out

def stable(g):
    """stable graphs: list (genera, edges, aut); own brute force"""
    res={}
    for nv in range(1,2*g-1):
        for gen in itertools.product(range(g+1),repeat=nv):
            if sum(gen)>g: continue
            ne=g-sum(gen)+nv-1
            prs=[(i,j) for i in range(nv) for j in range(i,nv)]
            for es in itertools.combinations_with_replacement(prs,ne):
                val=[0]*nv
                for u,v in es: val[u]+=1; val[v]+=1
                if any(2*gen[i]-2+val[i]<=0 for i in range(nv)): continue
                seen={0}
                ch=True
                while ch:
                    ch=False
                    for u,v in es:
                        if (u in seen)!=(v in seen):
                            seen|={u,v}; ch=True
                if len(seen)!=nv: continue
                orbit=Counter()
                for p in itertools.permutations(range(nv)):
                    gg=[0]*nv
                    for i in range(nv): gg[p[i]]=gen[i]
                    ee=tuple(sorted(tuple(sorted((p[u],p[v]))) for u,v in es))
                    orbit[(tuple(gg),ee)]+=1
                c=min(orbit)
                if c in res: continue
                aut=orbit[c]
                for e,m in Counter(c[1]).items():
                    aut*=math.factorial(m)*(2**m if e[0]==e[1] else 1)
                res[c]=aut
    return [(k[0],k[1],v) for k,v in res.items()]

class TF:
    """mode 'pt': point vertices, rows h^s Phat; mode 'hv': Omega_w vertices."""
    def __init__(s,cl,e1,e2,t,F,c1,mode='pt',hs=1,wts=None):
        s.e1=Fr(e1); s.e2=Fr(e2); s.E=s.e1+s.e2; E=s.E; s.kap=s.e1*s.e2
        s.t=Fr(t); s.F=Fr(F); s.c1=Fr(c1); s.mode=mode; N=cl.N; s.N=N
        y=s.t**3; s.Dc=1+27*y; Dc=s.Dc
        d2=E*(1/s.c1**2-1)/(2*Dc)
        if hs==0: h=[Fr(1)]+[Fr(0)]*N
        else:
            sg=hs
            h=mul_ser(mul_ser(theta_coeffs(-sg/E,N),theta_coeffs(sg/s.e1,N),N),theta_coeffs(sg/s.e2,N),N)
        s.h=h
        s.rows=[]; s.Del=[]
        for i in range(3):
            ti=(W**i)*s.t
            v=(1+3*ti)/(E*ti)
            def ev(lst):
                out=[]
                for p in lst:
                    x=Om(0)
                    for e,c in p.items(): x=x+(v**e)*c
                    out.append(x)
                return out
            S=ev(cl.S); Q=ev(cl.Q); TQ=ev(cl.TQ)
            P0=S
            P1=[x/s.c1 for x in Q]
            P2=[]
            for n in range(N+1):
                x=TQ[n]-Q[n]*d2
                if n>=1: x=x-Q[n-1]*s.F
                P2.append(x*(Dc*s.c1))
            def hm(a):
                out=[]
                for n in range(N+1):
                    x=Om(0)
                    for k in range(n+1):
                        if h[k]!=0: x=x+a[n-k]*h[k]
                    out.append(x)
                return out
            s.rows.append((hm(P0),hm(P1),hm(P2)))
            s.Del.append(3*s.kap*E/(v*v))
        k=s.kap
        s.eta={(0,0):k*E**3/9,(1,1):-3*k*E,(1,2):-3*k,(2,1):-3*k}
        s.hod=Hodge(wts) if mode=='hv' else None
        s.em={}; s.vm={}; s.pm={}
    def pair(s,i,j,M,b):
        key=(i,j,M,b)
        if key in s.pm: return s.pm[key]
        x=Om(0)
        for (a,c),e in s.eta.items():
            x=x+s.rows[i][a][M]*s.rows[j][c][b]*e
        s.pm[key]=x
        return x
    def edge(s,i,j,m,n):
        # E(psi,w)=[delta Del - <P(-psi),P(-w)>]/(psi+w), telescoping
        key=(i,j,m,n)
        if key in s.em: return s.em[key]
        x=Om(0)
        for k in range(n+1):
            M=m+k+1; b=n-k
            x=x+s.pair(i,j,M,b)*Fr(-1)*Fr((-1)**(M+b))*Fr((-1)**k)
        s.em[key]=x
        return x
    def dil(s,i,d):
        return s.rows[i][0][d-1]*Fr((-1)**d)
    def vert(s,i,g,hs):
        hs=tuple(sorted(hs)); key=(i,g,hs)
        if key in s.vm: return s.vm[key]
        n=len(hs); B_=3*g-3+n-sum(hs); tot=Om(0)
        if B_>=0:
            def multisets(rem,mx):
                yield ()
                for p in range(min(rem,mx),0,-1):
                    for r in multisets(rem-p,p): yield (p+1,)+r
            for ds in multisets(B_,B_):
                if s.mode=='pt':
                    if sum(d-1 for d in ds)!=B_: continue
                    val=tau(g,tuple(sorted(hs+ds)))
                else:
                    val=s.hod(g,tuple(sorted(hs+ds)))
                if val==0: continue
                sym=1
                for c in Counter(ds).values(): sym*=math.factorial(c)
                term=Om(val/sym)
                for d in ds: term=term*s.dil(i,d)
                tot=tot+term
        tot=tot*(s.Del[i]**(g-1))
        s.vm[key]=tot
        return tot
    def graph(s,gen,edges,colourings=None):
        nv=len(gen); val=[0]*nv
        for u,v in edges: val[u]+=1; val[v]+=1
        bud=[3*gen[i]-3+val[i] for i in range(nv)]
        total=Om(0)
        hs=[[] for _ in range(nv)]
        for col in (colourings or itertools.product(range(3),repeat=nv)):
            def go(k,pref):
                nonlocal total
                if k==len(edges):
                    x=pref
                    for i in range(nv):
                        x=x*s.vert(col[i],gen[i],tuple(hs[i]))
                        if x.zero(): return
                    total=total+x
                    return
                u,v=edges[k]
                ru=bud[u]-sum(hs[u])
                for a in range(ru+1):
                    rv=(bud[v]-sum(hs[v])) if u!=v else ru-a
                    for b in range(rv+1):
                        hs[u].append(a); hs[v].append(b)
                        go(k+1,pref*s.edge(col[u],col[v],a,b))
                        hs[u].pop(); hs[v].pop()
            go(0,Om(1))
        return total
    def total(s,g,per=False):
        t=Om(0); pl=[]
        for gen,ed,aut in stable(g):
            x=s.graph(gen,ed)*Fr(1,aut); pl.append(x); t=t+x
        return (t,pl) if per else t

def F2ref(e1,e2,y,F):
    X=1/(1+27*Fr(y)); F=Fr(F); E=Fr(e1)+Fr(e2); k=Fr(e1)*Fr(e2)
    f02=5*F**3/(8*X)+(5/(8*X)-Fr(1,2))*F**2+Fr(1,96)*(13*X+20/X-32)*F-(27*X**3-98*X**2+121*X-50)/(2160*X)
    f11=F**2/16+Fr(1,96)*(8-7*X)*F+(84*X**2-196*X+121)/5760
    f20=X/384*F-(21*X**2-29*X+8)/5760
    return E**4/k*f20+E**2*f11+k*f02

if __name__=="__main__":
    print("counts",[len(stable(g)) for g in (2,3)], "sum 1/aut g2",sum(Fr(1,a) for _,_,a in stable(2)))
    for (e1,e2,t,F,c1) in ((3,7,Fr(1,5),Fr(1,3),Fr(2,3)),(1,2,Fr(-2,7),Fr(5,4),Fr(-3)),(Fr(5,2),Fr(-3,7),Fr(3,11),Fr(-7,2),Fr(11,5))):
        E=Fr(e1)+Fr(e2)
        cl=Closed(E,10)
        wts=[E,-Fr(e1),-Fr(e2)]
        pt=TF(cl,e1,e2,t,F,c1,'pt',1); P=pt.total(2)
        hv=TF(cl,e1,e2,t,F,c1,'hv',0,wts); H=hv.total(2)
        hv2=TF(cl,e1,e2,t,F,Fr(c1)*7,'hv',0,wts); H2=hv2.total(2)
        print("g2",e1,e2,t,F,c1,"pt=",P,"hv-pt",H-P,"c1 indep",H2-H,"F2ref-pt",F2ref(e1,e2,Fr(t)**3,F)-P.a)
        # controls
        print("   ctrl hv(rows h)-pt nonzero:",not (TF(cl,e1,e2,t,F,c1,'hv',1,wts).total(2)-P).zero(), " pt(rows h^-1)-pt nonzero:",not (TF(cl,e1,e2,t,F,c1,'pt',-1).total(2)-P).zero())
