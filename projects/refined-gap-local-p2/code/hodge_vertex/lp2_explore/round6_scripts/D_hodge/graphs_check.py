# Direct evaluation of T(Gamma,c) of (0.1) from scratch. Evidence only.
import sys,time
from fractions import Fraction as Fr
from itertools import product
from math import factorial
from collections import Counter
from hodge import H, odd_multisets
from rows import *

class Vtx:
    def __init__(s,th):
        s.th=th
        w=[th.E,-th.e1,-th.e2]
        s.ck={k:-factorial(k-1)*sum(Fr(1)/x**k for x in w) for k in range(1,30,2)}
        s.cache={}
    def omega_int(s,g,exps):
        exps=tuple(sorted(exps)); key=(g,exps)
        if key in s.cache: return s.cache[key]
        n=len(exps); dim=3*g-3+n; rem=dim-sum(exps)
        tot=Fr(0)
        if rem>=0 and 2*g-2+n>0:
            for K in odd_multisets(rem,max(1,2*g-1)):
                co=Fr(1)
                for k,mu in Counter(K).items():
                    co*=s.ck[k]**mu/factorial(mu)
                tot+=co*H(g,exps,tuple(sorted(K)))
        s.cache[key]=tot
        return tot
    def decor(s,g,n_exps):
        n=len(n_exps); base=3*g-3+n-sum(n_exps)
        out=[]
        m=0
        while 2*m<=base+m:
            def rec(m_,bud):
                if m_==0:
                    yield ()
                    return
                for d in range(2,bud-2*(m_-1)+1):
                    for r in rec(m_-1,bud-d):
                        yield (d,)+r
            out.extend(rec(m,base+m))
            m+=1
        return out
    def vertex_nodelta(s,g,colour,exps):
        key=('Vn',g,colour,tuple(sorted(exps)))
        if key in s.cache: return s.cache[key]
        th=s.th
        acc=None
        for d in s.decor(g,exps):
            I=s.omega_int(g,tuple(exps)+d)
            if I==0: continue
            t=L.const(fp(I)*inv(factorial(len(d)))%P,th.prec)
            for x in d:
                t=t*th.dil[colour][x]
            acc=t if acc is None else acc+t
        if acc is None: acc=L.const(0,th.prec)
        s.cache[key]=acc
        return acc
    def vertex(s,g,colour,exps):
        key=('V',g,colour,tuple(sorted(exps)))
        if key in s.cache: return s.cache[key]
        th=s.th
        acc=s.vertex_nodelta(g,colour,exps)
        D=th.Delta[colour]
        if g>=1:
            pw=L.const(1,th.prec)
            for _ in range(g-1): pw=pw*D
        else:
            pw=D.inverse()
        r=acc*pw
        s.cache[key]=r
        return r

def graph_value(th,V,genera,edges,col):
    nh=[0]*len(genera)
    for (u,v) in edges:
        nh[u]+=1; nh[v]+=1
    dims=[3*genera[v]-3+nh[v] for v in range(len(genera))]
    he=[]
    for e,(u,v) in enumerate(edges):
        he.append(u); he.append(v)
    tot=None
    def rec(i,cur,used):
        nonlocal tot
        if i==len(he):
            term=None
            for v in range(len(genera)):
                ex=tuple(cur[k] for k in range(len(he)) if he[k]==v)
                f=V.vertex(genera[v],col[v],ex)
                term=f if term is None else term*f
            for e,(u,w) in enumerate(edges):
                term=term*th.edge(col[u],col[w],cur[2*e],cur[2*e+1])
            tot=term if tot is None else tot+term
            return
        v=he[i]
        for a in range(0,dims[v]-used[v]+1):
            used[v]+=a
            cur.append(a)
            rec(i+1,cur,used)
            cur.pop()
            used[v]-=a
    rec(0,[],[0]*len(genera))
    return tot

def recon(x):
    r0,r1=P,x%P; s0,s1=0,1
    while r1>(1<<30):
        q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    if s1==0: return None
    if s1<0: r1,s1=-r1,-s1
    return Fr(r1,s1)

GR3=[
 ("g2 loop",[2],[(0,0)]),
 ("g1 2 loops",[1],[(0,0),(0,0)]),
 ("g0 3 loops",[0],[(0,0),(0,0),(0,0)]),
 ("g2 - g1",[2,1],[(0,1)]),
 ("g2 - g0(loop)",[2,0],[(0,1),(1,1)]),
 ("g1 = g1",[1,1],[(0,1),(0,1)]),
 ("g1(loop) - g1",[1,1],[(0,1),(1,1)]),
 ("g0 =3= g1",[0,1],[(0,1),(0,1),(0,1)]),
 ("g1 - g1 - g1",[1,1,1],[(0,1),(1,2)]),
 ("g1 - g1 - g0(loop)",[1,1,0],[(0,1),(1,2),(2,2)]),
]
if __name__=="__main__":
    e1,e2=Fr(sys.argv[1]),Fr(sys.argv[2])
    mode=sys.argv[3] if len(sys.argv)>3 else "on"
    t0=time.time()
    kw={}
    if mode=="off": kw=dict(c1factor=[1,Fr(1,10)])
    if mode=="r2": kw=dict(rperturb={2:Fr(3,81)})
    th=Theory(e1,e2,N=8,K=9,**kw)
    V=Vtx(th)
    print("built",round(time.time()-t0,1),"mode",mode)
    print("Delta0 t^-2 coef:",recon(th.Delta[0].coef(-2)),"expected",Fr(3)*th.e1*th.e2*th.E**3)
    for i,j in [(0,0),(0,1),(1,2),(1,1)]:
        ok=True
        for n in range(0,8):
            acc=None
            for m in range(n+1):
                for (a,b),h in th.eta.items():
                    t=th.rows[i][a][m]*th.rows[j][b][n-m]*h*((-1)**(n-m))
                    acc=t if acc is None else acc+t
            exp=(th.Delta[i] if (i==j and n==0) else th.Delta[i]*0)
            diff=acc-exp
            for k in range(diff.v,3):
                if diff.coef(k)!=0: ok=False
        print("Prop 2.10 pair",i,j,ok)
    print("edge symmetry 0,1:",all((th.edge(0,1,m,n)-th.edge(1,0,n,m)).coef(k)==0 for m in range(3) for n in range(3) for k in range(-2,3)))
    if mode=="on":
        sv=V.vertex(3,0,())
        print("single vertex genus3 colour0 coefs t^-4..t^-1:",[recon(sv.coef(k)) for k in range(-4,0)])
        print("single vertex genus3 colour1 val",V.vertex(3,1,()).norm().v)
    for name,gen,ed in GR3:
        res={}
        t1=time.time()
        worst=10**9; wc=None
        for col in product(range(3),repeat=len(gen)):
            val=graph_value(th,V,gen,ed,col).norm()
            if val.v<worst: worst=val.v; wc=col
            res[col]=val
        print("%s: min valuation over %d colourings = %d at %s; time %.1fs"%(name,3**len(gen),worst,wc,time.time()-t1))
        sys.stdout.flush()
