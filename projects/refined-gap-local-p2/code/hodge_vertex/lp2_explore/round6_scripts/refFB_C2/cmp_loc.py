# Direct test: T_F graph sum (no stored data, F = theta log A from the hypergeometric series) against
# virtual localisation at generic alpha, interpolated to alpha = E/3. Evidence only.
import sys, time, math
from fractions import Fraction as Fr
from loc import Loc, TF_limit
from tf import Closed, TF, F2ref

def series(N):
    """A,F,log Q - log y as y-coefficients up to y^N"""
    a=[Fr(1)]
    for d in range(1,N+1):
        a.append(Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d)
    th=[d*a[d] for d in range(N+1)]
    # F=thA/A
    F=[Fr(0)]*(N+1)
    for n in range(N+1):
        F[n]=th[n]-sum(F[k]*a[n-k] for k in range(n))
    g=[Fr(0)]+[3*Fr(math.factorial(3*d-1),math.factorial(d)**3)*(-1)**d for d in range(1,N+1)]
    return a,F,g
def evalser(c,y):
    r=Fr(0)
    for x in reversed(c): r=r*y+x
    return r
def Qofy(y,g,ORD=12):
    # exp of g(y) truncated (exact rational truncation; y small)
    x=evalser(g,y); s=Fr(1); term=Fr(1)
    for k in range(1,ORD+1):
        term=term*x/k; s+=term
    return y*s

if __name__=="__main__":
    g=int(sys.argv[1]); Dmax=int(sys.argv[2])
    eps=[(3,7),(1,2),(Fr(5,2),Fr(-3,7))]
    if len(sys.argv)>3: eps=eps[:int(sys.argv[3])]
    a,F,gs=series(30)
    for e1,e2 in eps:
        t0=time.time()
        vals=[];oks=[]
        for d in range(0,Dmax+1):
            r,ok=TF_limit(g,d,e1,e2)
            vals.append(r); oks.append(ok)
        print(f"g={g} eps=({e1},{e2}) GW_(g,d), d=0..{Dmax}:",[str(v) for v in vals],"polynomial in alpha:",oks,f"[{time.time()-t0:.0f}s]",flush=True)
        E=Fr(e1)+Fr(e2)
        cl=Closed(E,3*g+1)
        wts=[E,-Fr(e1),-Fr(e2)]
        for t in (Fr(1,50),Fr(-1,50),Fr(1,500),Fr(-1,500)):
            y=t**3
            Fv=evalser(F,y); Av=evalser(a,y)
            pt=TF(cl,e1,e2,t,Fv,Av,'pt',1).total(g)
            hv=TF(cl,e1,e2,t,Fv,Fr(1),'hv',0,wts).total(g)
            Q=Qofy(y,gs)
            S=sum(vals[d]*Q**d for d in range(Dmax))   # up to d=Dmax-1
            rat=(pt.a-S)/Q**Dmax
            print(f"   t={t}: im={pt.b}, pt-hv={(pt-hv).a}, (Tbar-Sum_(d<{Dmax}))/Q^{Dmax} ={float(rat):.10f}; GW_(g,{Dmax})={float(vals[Dmax]):.10f}",flush=True)
        # controls: wrong h exponent on rows
        t=Fr(1,500); y=t**3; Fv=evalser(F,y); Av=evalser(a,y); Q=Qofy(y,gs)
        ctrl=TF(cl,e1,e2,t,Fv,Av,'pt',-1).total(g)
        ctrl0=TF(cl,e1,e2,t,Fv,Av,'pt',0).total(g)
        ctrlF=TF(cl,e1,e2,t,Fv+Fr(1,9),Av,'pt',1).total(g)
        print("   controls at t=1/500: rows h^-1, GW0 mismatch:",float(ctrl.a-vals[0]),"; rows h^0:",float(ctrl0.a-vals[0]),"; F+1/9:",float(ctrlF.a-vals[0]),flush=True)
