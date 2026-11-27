# Evidence only. 3d formal stationary phase of exp(W/z) dOmega at colour i, generic alpha, y-series through y^2.
# Test: omega(y,z)/z-free part = (rho/rho_0)(y,z) * prod_3 Theta(-z/a), through z^3; control with Theta(+z/a).
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2d1 import *
from functools import lru_cache

ZTOP=3
def mat_inv3(H):
    # H 3x3 of y-series; Gauss-Jordan over the series ring (pivots units at y=0 after pivoting by constant term)
    n=3
    A=[[list(H[i][j]) for j in range(n)]+[s0(1) if i==j else s0(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c][0]!=0)
        A[c],A[p]=A[p],A[c]
        iv=sinv(A[c][c])
        A[c]=[smul(iv,x) for x in A[c]]
        for r in range(n):
            if r!=c:
                f=A[r][c]
                A[r]=[ssub(x,smul(f,y)) for x,y in zip(A[r],A[c])]
    return [row[n:] for row in A]

def pmul(a,b,dmax=None):
    r={}
    for ma,ca in a.items():
        for mb,cb in b.items():
            m=(ma[0]+mb[0],ma[1]+mb[1],ma[2]+mb[2])
            if m in r: r[m]=sadd(r[m],smul(ca,cb))
            else: r[m]=smul(ca,cb)
    return r
def padd(a,b,c=Fr(1)):
    r=dict(a)
    for m,x in b.items():
        y=sscal(c,x)
        r[m]=sadd(r[m],y) if m in r else y
    return r

def run(al,e1,e2,i):
    T=Theory(al,e1,e2)
    rh,Y=T.rho(i)
    a=T.al
    q=[sadd(Y,s0(a[j])) for j in range(3)]
    P=sscal(Fr(-3),Y)
    # Hessian of W in log coordinates
    H=[[sscal(Fr(1,9),P) for _ in range(3)] for _ in range(3)]
    for j in range(3): H[j][j]=sadd(H[j][j],q[j])
    Hinv=mat_inv3(H)
    G=[[sscal(Fr(-1),Hinv[r][c]) for c in range(3)] for r in range(3)]  # covariance for phi=-W
    def T_m(m):
        out={}
        for j in range(3):
            e=[0,0,0]; e[j]=m
            out[tuple(e)]=sscal(Fr(1,math.factorial(m)),q[j])
        # P*(sum xi/3)^m/m!
        for e0 in range(m+1):
            for e1_ in range(m+1-e0):
                e2_=m-e0-e1_
                mult=Fr(math.factorial(m),math.factorial(e0)*math.factorial(e1_)*math.factorial(e2_))
                c=sscal(mult*Fr(1,3**m)/math.factorial(m),P)
                k=(e0,e1_,e2_)
                out[k]=sadd(out[k],c) if k in out else c
        return out
    Ts={m:T_m(m) for m in range(3,2*ZTOP+3)}
    Es={0:{(0,0,0):s0(1)}}
    for n in range(1,2*ZTOP+1):
        cur={}
        for k in range(1,n+1):
            cur=padd(cur,pmul(Ts[k+2],Es[n-k]),Fr(k,n))
        Es[n]=cur
    @lru_cache(None)
    def mom(b):
        if sum(b)==0: return tuple(s0(1))
        if sum(b)%2: return tuple(s0(0))
        c=next(j for j in range(3) if b[j]>0)
        tot=s0(0)
        for d in range(3):
            nb=list(b); nb[c]-=1
            if nb[d]==0: continue
            mult=nb[d]; nb[d]-=1
            tot=sadd(tot,sscal(Fr(mult),smul(G[c][d],list(mom(tuple(nb))))))
        return tuple(tot)
    om={}
    for n in range(0,ZTOP+1):
        tot=s0(0)
        for m,c in Es[2*n].items():
            tot=sadd(tot,smul(c,list(mom(m))))
        om[n]=tot
    # expected: (rho/rho0)*prod Theta(-z/a)
    w3=[3*a[i]]+[a[j]-a[i] for j in range(3) if j!=i]
    def kappa(sign):
        h={0:Fr(1)}
        for w in w3:
            h=lmul(h,theta_exp(sign,w,ZTOP),ZTOP)
        return h
    rho0=[rh[d].get(0,Fr(0)) for d in range(NY+1)]
    irho0=sinv(rho0)
    rr={p:smul(irho0,[rh[d].get(p,Fr(0)) for d in range(NY+1)]) for p in range(ZTOP+1)}
    def expected(sign):
        k=kappa(sign); out={}
        for p in range(ZTOP+1):
            t=s0(0)
            for q_ in range(p+1):
                t=sadd(t,sscal(k.get(q_,Fr(0)),rr[p-q_]))
            out[p]=t
        return out
    ex=expected(-1); exc=expected(+1)
    ok=all(om[p]==ex[p] for p in range(ZTOP+1))
    okc=all(om[p]==exc[p] for p in range(ZTOP+1))
    # z^0 prefactor check: rho_0^2 * W(y)/W(0)=1
    W=T.Wfn(Y)
    chk=smul(smul(rho0,rho0),sscal(1/W[0],W))==s0(1)
    return ok,okc,chk,om[1][:2]

if __name__=="__main__":
    for al,e1,e2 in [((1,2,5),3,7),((2,-3,7),1,4)]:
        for i in range(3):
            print(al,"colour",i,"(match,control Theta(+z/a),rho0^2 W/W0=1,om_1):",run(al,e1,e2,i))
