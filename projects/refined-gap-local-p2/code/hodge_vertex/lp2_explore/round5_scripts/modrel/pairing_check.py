# Evidence only (exact). Checks Pi_ij = delta_ij for the closed form, as Laurent polynomials in w = 1/t over Q(omega),
# using ONLY the closed form (no BS). Frame B (no h_w). Includes negative controls.
from fractions import Fraction as Fr
def run(e1,e2,N,perturb=None,eta00factor=Fr(1)):
    E=e1+e2; eps=e1*e2
    r=[Fr(1)]
    for m in range(0,N+3):
        nxt=(3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>=1 else 0))/(27*(m+1))
        r.append(nxt)
    if perturb is not None:
        r[perturb]+=Fr(1,1000)
    ZERO=(Fr(0),Fr(0))
    def cm(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
    def ca(x,y): return (x[0]+y[0],x[1]+y[1])
    def cs(c,x): return (c*x[0],c*x[1])
    om=[(Fr(1),Fr(0)),(Fr(0),Fr(1)),(Fr(-1),Fr(-1))]
    def omp(k): return om[k%3]
    def Ladd(p,q):
        r_=dict(p)
        for k,v in q.items():
            r_[k]=ca(r_.get(k,ZERO),v)
        return {k:v for k,v in r_.items() if v!=ZERO}
    def Lmul(p,q):
        r_={}
        for k,v in p.items():
            for l,u in q.items():
                r_[k+l]=ca(r_.get(k+l,ZERO),cm(v,u))
        return {k:v for k,v in r_.items() if v!=ZERO}
    def Ls(c,p): return {k:cs(c,v) for k,v in p.items()} if c!=0 else {}
    def Lth(p): return {k:cs(Fr(-k,3),v) for k,v in p.items() if k!=0}
    one={0:(Fr(1),Fr(0))}
    def Zmul(a,b):
        r_=[{} for _ in range(N+1)]
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                if i+j<=N: r_[i+j]=Ladd(r_[i+j],Lmul(x,y))
        return r_
    def Zadd(a,b): return [Ladd(x,y) for x,y in zip(a,b)]
    def Zs(p,a): return [Lmul(p,x) for x in a]
    def Zth(a): return [{}]+[Lth(a[k]) for k in range(N)]
    def Zflip(a): return [Ls(Fr((-1)**k),x) for k,x in enumerate(a)]
    def v_of(i): return {1:cs(1/E,omp(-i)),0:cs(3/E,om[0])}
    def zfrac(l):
        s=[Fr(0)]*(N+1)
        if l<=N: s[l]=Fr(1)
        for j in range(1,l+2):
            for k in range(1,N+1):
                s[k]-=s[k-1]*Fr(j)/E
        return s
    def S_of(i):
        v=v_of(i); res=[{} for _ in range(N+1)]; vp=one
        for l in range(N+1):
            zf=zfrac(l)
            for k in range(N+1):
                if zf[k]!=0:
                    res[k]=Ladd(res[k],Ls(zf[k]*r[l],vp))
            vp=Lmul(vp,v)
        return res
    Delta_c={0:om[0],-3:(Fr(27),Fr(0))}
    def rows(i):
        v=v_of(i); S=S_of(i)
        P0=Zs(v,S)
        Q1=Zadd(Zadd([Ls(-E/3,x) for x in P0],S),Zth(P0))
        ti={-1:omp(i)}
        DY=Ladd(Ls(-E/3,Delta_c), Ls(E,Lmul(ti,Ladd(Ladd(one,Ls(Fr(-3),ti)),Ls(Fr(9),Lmul(ti,ti))))))
        T=Zadd(Zadd(Zs(DY,Q1),Zs(Delta_c,Zth(Q1))),[Ls(E/2,x) for x in Q1])
        return P0,Q1,T
    R={i:rows(i) for i in range(3)}
    bad=[]
    for i in range(3):
        for j in range(i,3):
            P,Q,T=R[i]; Pp,Qp,Tp=[Zflip(x) for x in R[j]]
            Pi=Zadd([Ls(eps*E**3/9*eta00factor,x) for x in Zmul(P,Pp)],
                    [Ls(-3*eps,x) for x in Zadd(Zmul(Q,Tp),Zmul(T,Qp))])
            for k in range(N+1):
                target={0:(3*eps*E,Fr(0))} if (i==j and k==0) else {}
                if Pi[k]!=target:
                    bad.append((i,j,k))
    return (len(bad)==0),bad[:4]
for (e1,e2,N) in [(Fr(1),Fr(1),6),(Fr(5),Fr(2),6),(Fr(7),Fr(-3),5),(Fr(1,2),Fr(3,5),5)]:
    print("eps",e1,e2,"N",N,run(e1,e2,N))
print("control: perturbed r_2 ->",run(Fr(5),Fr(2),5,perturb=2))
print("control: wrong eta^00 ->",run(Fr(5),Fr(2),4,eta00factor=Fr(2)))
