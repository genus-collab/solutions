# (1) Feynman expansion of the oscillating integral vs rows derived from I (A3.5, Stage B Phi)
# (2) Hodge-vertex form of the graph sum vs point form
from gsum import *
from engine import hodge_norm
from functools import lru_cache
import sys,time

class YS:
    __slots__=['c']
    D=2
    def __init__(s,c): s.c=c
    def __add__(s,o): return YS([x+w for x,w in zip(s.c,o.c)])
    def __mul__(s,o): return YS(ymul(s.c,o.c,YS.D))
    def scal(s,k): return YS([k*x for x in s.c])
    def iszero(s): return all(x==0 for x in s.c)

def const(k): return YS([Fr(k)]+[Fr(0)]*YS.D)

def phi_test(al,ep,D,N):
    YS.D=D
    c=setup(al,ep,D,N+2)
    res=[]
    for i in range(3):
        Y=YS(c.Y[i])
        d=[Y+const(a) for a in al]
        # H=diag(d)-(Y/3)J ; inverse via adjugate
        Y3=Y.scal(Fr(1,3))
        H=[[ (d[a] if a==b else const(0)) + Y3.scal(Fr(-1)) for b in range(3)] for a in range(3)]
        def det2(a,b,c_,d_): return a*d_+(b*c_).scal(Fr(-1))
        cof=[[None]*3 for _ in range(3)]
        for a in range(3):
            for b in range(3):
                r=[x for x in range(3) if x!=a]; s=[x for x in range(3) if x!=b]
                m=det2(H[r[0]][s[0]],H[r[0]][s[1]],H[r[1]][s[0]],H[r[1]][s[1]])
                cof[a][b]=m.scal(Fr((-1)**(a+b)))
        det=H[0][0]*cof[0][0]+H[0][1]*cof[0][1]+H[0][2]*cof[0][2]
        # det H=W/3 check
        assert [3*x for x in det.c]==c.Wval[i], "detH!=W/3"
        dinv=YS(yinv(det.c,D))
        G=[[ (cof[b][a]*dinv).scal(Fr(-1)) for b in range(3)] for a in range(3)]
        # polys in xi: dict (b0,b1,b2)->YS ; T_m
        from math import comb
        def Tm(m):
            P={}
            for j in range(3):
                e=[0,0,0]; e[j]=m
                P[tuple(e)]=d[j].scal(Fr(1,factorial(m)))
            # -3Y (sum xi/3)^m/m!
            for b0 in range(m+1):
                for b1 in range(m+1-b0):
                    b2=m-b0-b1
                    co=Fr(factorial(m),factorial(b0)*factorial(b1)*factorial(b2))/3**m/factorial(m)*(-3)
                    k=(b0,b1,b2)
                    t=Y.scal(co)
                    P[k]=P[k]+t if k in P else t
            return P
        S=2*N  # powers of s=sqrt z
        Ts={m:Tm(m) for m in range(3,S+3)}
        # E: list by s-power of poly
        def pmul(A,B):
            R={}
            for k1,v1 in A.items():
                for k2,v2 in B.items():
                    k=(k1[0]+k2[0],k1[1]+k2[1],k1[2]+k2[2])
                    t=v1*v2
                    R[k]=R[k]+t if k in R else t
            return R
        # exp via recursion in s: n E_n=sum_{m} (m-2) T_m E_{n-(m-2)}
        E=[{} for _ in range(S+1)]; E[0]={(0,0,0):const(1)}
        for n in range(1,S+1):
            acc={}
            for mm in range(1,n+1):
                m=mm+2
                t=pmul(Ts[m],E[n-mm])
                for k,v in t.items():
                    v=v.scal(Fr(mm))
                    acc[k]=acc[k]+v if k in acc else v
            E[n]={k:v.scal(Fr(1,n)) for k,v in acc.items()}
        @lru_cache(maxsize=None)
        def mom(b):
            if sum(b)%2: return const(0)
            if sum(b)==0: return const(1)
            a=next(x for x in range(3) if b[x]>0)
            bb=list(b); bb[a]-=1
            tot=const(0)
            for j in range(3):
                if bb[j]>0:
                    b2=list(bb); b2[j]-=1
                    tot=tot+(G[a][j]*mom(tuple(b2))).scal(Fr(bb[j]))
            return tot
        out=[]
        for n in range(0,N+1):
            tot=const(0)
            for k,v in E[2*n].items():
                # xi has covariance z G -> s^{|b|}; E uses T_m s^{m-2}*xi'^m with xi=s xi'
                tot=tot+v*mom(k)
            out.append(tot)
        # compare with Pt_0 h_{3a,tau,tau'}
        hh=h_series(c.wts[i][:3],N+1)
        hz=[{k:x for k,x in enumerate(hh) if x!=0}]+[{} for _ in range(D)]
        Pp=yzmul(c.Pt[i][0],hz,D,N+1)
        ok=True
        for n in range(N+1):
            tgt=[Pp[y].get(n,Fr(0)) for y in range(D+1)]
            if tgt!=out[n].c: ok=False
        res.append(ok)
        # control: Theta with other sign
        hh2=h_series(c.wts[i][:3],N+1,sign=-1)
        hz2=[{k:x for k,x in enumerate(hh2) if x!=0}]+[{} for _ in range(D)]
        Pp2=yzmul(c.Pt[i][0],hz2,D,N+1)
        res.append(('ctrl',all([Pp2[y].get(n,Fr(0)) for y in range(D+1)]==out[n].c for n in range(N+1))))
    return res

def hodge_form(al,ep,g,D):
    """graph sum with rows Pt (no h) and vertex Delta^{g_v-1} Hodge_wt"""
    c=setup(al,ep,D,3*g)
    # replace P by Pt
    c.P=c.Pt
    gs=GS(c)
    def L(i,gv,h):
        key=(i,gv,h)
        if key in gs.Lcache: return gs.Lcache[key]
        n=len(h); res=[Fr(0)]*(D+1)
        if stable(gv,n):
            rem=3*gv-3+n-sum(h)
            if rem>=0:
                from gsum import partitions_min2
                # leaves: total psi+Hodge degree free: leaves d_k>=2 with sum(d_k-1)=r<=rem
                for r in range(0,rem+1):
                    for part in partitions_min2(r):
                        ds=tuple(p+1 for p in part)
                        v=hodge_norm(gv,tuple(sorted(h+ds)),tuple(c.wts[i]))
                        if v==0: continue
                        t=[v]+[Fr(0)]*D
                        for dd,mu in Counter(ds).items():
                            for _ in range(mu): t=ymul(t,gs.Dil[i][dd],D)
                            t=[x/factorial(mu) for x in t]
                        res=yadd(res,t)
        res=ymul(res,gs.Dpow[(i,gv-1)],D)
        gs.Lcache[key]=res
        return res
    gs.L=L
    return gs.total(g)

if __name__=="__main__":
    al,ep=(Fr(1),Fr(2),Fr(5)),(Fr(3),Fr(7))
    t=time.time()
    print("Phi vs Pt0*h_{3a,tau,tau'}, z^<=3, y^<=2 :",phi_test(al,ep,2,3),time.time()-t,flush=True)
    al2,ep2=(Fr(2),Fr(-3),Fr(8)),(Fr(1),Fr(4))
    print(phi_test(al2,ep2,2,3),flush=True)
    for g in (2,3):
        F=hodge_form(al,ep,g,2)
        c=setup(al,ep,2,3*g); G=GS(c).total(g)
        print("Hodge-vertex form == point form, g=",g,F==G,flush=True)
