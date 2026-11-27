import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cylib import *
import cylib
# ---- Hodge integrals: known values
print("Hodge checks")
checks=[("l1 M11",int_lam(1,[1],(0,)),Fr(1,24)),("l1 psi^3 M21",int_lam(2,[1],(3,)),Fr(1,480)),("l2 psi^2 M21",int_lam(2,[2],(2,)),Fr(7,5760)),
("l1l2 M2",int_lam(2,[1,2],()),Fr(1,5760)),("ch3 M2",hi(2,(),(3,)),Fr(-1,34560)),("l1^3 M2",int_lam(2,[1,1,1],()),Fr(1,2880)),
("l3 psi^4 M31",int_lam(3,[3],(4,)),Fr(31,967680)),("l3l2 psi^2 M31",int_lam(3,[3,2],(2,)),Fr(1,120960)),("l1l2l3 M3",int_lam(3,[1,2,3],()),Fr(1,1451520)),
("l2^3 M3",int_lam(3,[2,2,2],()),Fr(1,725760))]
for n,a,b in checks: print(" ",n,a,b,"OK" if a==b else "FAIL")
cylib.USE_CUTOFF[0]=False
cylib._hi.cache_clear()
print(" raw (no cutoff) l1^6 M3 =",int_lam(3,[1]*6,()),"expect 1/90720")
print(" raw ch5 M3 =",hi(3,(),(5,1)),"; ch_7 on M_{3,1}, vanish check:", hi(3,(0,),(7,)), hi(2,(0,),(4-1+0,)) )
print(" raw lambda_1^2 on M_{1,2}:",hi(1,(0,0),(1,1)))
print(" raw vanish l1^4 M_{2,1}:",hi(2,(0,),(1,1,1,1)))
cylib.USE_CUTOFF[0]=True
cylib._hi.cache_clear()
# ---- Decorated graphs counts
for g in (2,3):
    t0=time.time()
    for hodge in (False,True):
        d=Decorated(g,hodge); print("g",g,"hodge",hodge,d.counts(),round(time.time()-t0,1),"s")
# ---- QRR constants at Q=0
print("Q->0 normalisation")
def wick_from_W(Wk,zmax):
    out=[1]
    W2i=inv(Wk[2])
    for m in range(1,zmax+1):
        tot=0
        for cf,n,pr in wick_terms(m):
            t=fr2p(cf)*pow(W2i,n,P)%P
            for k,mu in pr: t=t*pow(Wk[k],mu,P)%P
            tot+=t
        out.append(tot%P)
    return out
def peval(pc,x):
    r=0
    for c in reversed(pc): r=(r*x+c)%P
    return r
zmax=9
for (a1,a2,e1,e2) in [(3,4,5,11),(7,13,2,29),(17,-5,19,-8)]:
    W=cyw(a1,a2,e1,e2)
    a0,a1_,a2_,e1_,e2_=W['a0'],W['a1'],W['a2'],W['e1'],W['e2']
    Ps=thetapolys(2*zmax+3)
    ffc=ffc_closed(W,zmax)
    # root 1: U0=-a2/a0
    U0=(-a2_*inv(a0))%P; X2=U0*inv(1-U0)%P
    Wk={k:(-(a0+a2_)*peval(Ps[k],X2))%P for k in range(2,2*zmax+3)}
    w1=wick_from_W(Wk,zmax)
    # root 2: X1=a0/a1, X2=-1
    X1=a0*inv(a1_)%P
    Wk={k:(-a1_*peval(Ps[k],X1)-(a0+a2_)*peval(Ps[k],P-1))%P for k in range(2,2*zmax+3)}
    w2=wick_from_W(Wk,zmax)
    def conv(a,b): return [sum(a[i]*b[m-i] for i in range(m+1))%P for m in range(zmax+1)]
    t1=qrr_target([a0,a1_,a2_,-e1_%P,-e2_%P],zmax)
    t2=qrr_target([-a0%P,(a0+a1_)%P,(a0+a2_)%P,-e1_%P,-e2_%P],zmax)
    print(" ",(a1,a2,e1,e2),"root1 ok:",conv(w1,ffc)==t1,"root2 ok:",conv(w2,ffc)==t2)
    # wrong sign of constant fails
    ffc2=ffc[:]; 
    for i in range(1,zmax+1,2): pass
    print("   wrong constant (inverse) gives equality:",conv(w1,series_exp([0]+[(-x)%P for x in [0]*zmax],zmax))==t1)
