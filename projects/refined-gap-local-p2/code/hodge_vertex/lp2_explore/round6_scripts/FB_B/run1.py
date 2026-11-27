# Evidence only. E1: Phi(alpha*,Y;z) = Theta(-z/E) S(v,z). E2: T_F total vs BS eq:F2ref (independent of BS: BS's formula is only a comparison).
# E2b: pairing with explicit rows (A arbitrary) equals the A-free formula. Symplectic defect at z^1..3.
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from trace_total import *

def Sz_in_C(R,E,ZT):
    v=R.inv(R.add(R.scal(E/3,R.one),R.Y))  # v=3/(E+3Y)=1/(E/3+Y)
    r=[Fr(0),Fr(1)]
    for m in range(0,ZT+2): r.append((3*(3*m*m+3*m+1)*r[m+1]-m**3*r[m])/(27*(m+1)))
    Sz=[R.zero for _ in range(ZT+1)]
    vp=R.one
    for l in range(ZT+1):
        den=[Fr(1)]+[Fr(0)]*ZT
        for j in range(1,l+2):
            den=zmul(den,[Fr(1),Fr(j)/E]+[Fr(0)]*(ZT-1),ZT)
        inv=[Fr(0)]*(ZT+1); inv[0]=Fr(1)
        for n in range(1,ZT+1):
            inv[n]=-sum(den[k]*inv[n-k] for k in range(1,n+1))
        for n in range(ZT+1-l):
            Sz[n+l]=R.add(Sz[n+l],R.scal(r[l+1]*inv[n],vp))
        vp=R.mul(vp,v)
    return Sz

print("E1: Phi at T_F vs Theta(-z/E) S(v,z)",flush=True)
for (e1,e2,y) in [(Fr(3),Fr(7),Fr(1,5)),(Fr(1),Fr(1),Fr(-2,7)),(Fr(5,2),Fr(-3,7),Fr(3,11))]:
    E=e1+e2; R=Ring([E/3]*3,y)
    ZT=3
    J,Phi,W0,iW0=wick_rows(R,ZT)
    Sz=Sz_in_C(R,E,ZT)
    th=theta_series(-1/E,ZT); thc=theta_series(1/E,ZT)
    RHS=[R.zero for _ in range(ZT+1)]; RHSc=[R.zero for _ in range(ZT+1)]
    for n in range(ZT+1):
        for k in range(n+1):
            RHS[n]=R.add(RHS[n],R.scal(th[k],Sz[n-k]))
            RHSc[n]=R.add(RHSc[n],R.scal(thc[k],Sz[n-k]))
    print(" eps",e1,e2,"y",y,"match z^0..3:",[Phi[n][0]==RHS[n] for n in range(ZT+1)],"control Theta(+z/E):",[Phi[n][0]==RHSc[n] for n in range(ZT+1)],flush=True)

print("E2: T_F total vs BS eq:F2ref (F indeterminate, values chosen freely)",flush=True)
for (e1,e2,y) in [(Fr(3),Fr(7),Fr(1,5)),(Fr(1),Fr(1),Fr(-2,7)),(Fr(5,2),Fr(-3,7),Fr(3,11))]:
    E=e1+e2
    t0=time.time()
    T=Total([E/3]*3,e1,e2,y)
    print(" setup",round(time.time()-t0,1),"s; symplectic defect z^1,2,3:",[T.sympl_defect(n,Fr(5,3)) for n in (1,2,3)],flush=True)
    A0=Fr(7,4); F0=Fr(5,3)
    same=all(T.pairing(M,b,F0)==T.pairing_explicit(M,b,F0,A0) for M in range(4) for b in range(4))
    print(" A-free pairing = explicit rows at A=7/4:",same,flush=True)
    for F in [Fr(0),Fr(1),Fr(-2,3),Fr(5,3),Fr(11,7)]:
        val=T.F2(F); bs=bs_F2ref(e1,e2,y,F)
        print("  F",F,"total == BS eq:F2ref:",val==bs, "" if val==bs else (float(val),float(bs)),flush=True)
    # control: wrong sign in Theta(z/eps)
    Tc=Total([E/3]*3,e1,e2,y,hsign=-1)
    print(" control hsign=-1, F=1:",Tc.F2(Fr(1))==bs_F2ref(e1,e2,y,Fr(1)),flush=True)
