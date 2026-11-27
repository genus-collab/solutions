# Prop 2.1: Schur side vs Hodge engine (MV-free) + vanishing of the other powers of x; and small exact s.8.2 tests.
from schur_side import *
from hodge import *
import sys,time
K=Fld(None)
taus=[Fr(2,7),Fr(-5,3),Fr(-1,2),Fr(-1,3),Fr(13,9)]
tot=0;bad=0;bad2=0;ctot=0
for tau in taus:
    sig=tau*(1+tau); z=z_of_tau(tau)
    for R,amax in [((),5),((2,),4),((3,),3),((2,2),3),((2,3),2),((4,),2),((2,2,2),1)]:
        P=sum(R)+amax+ 2*3+ len(R)+amax+2  # up to g=3
        F,items=compute_F(R,amax,P,tau,K)
        for ev,a in F:
            sub=()
            for (v,_),e in zip(items,ev): sub+=(v,)*e
            l=len(sub)+a
            W=sum(sub)+a
            mu=sub+(1,)*a
            for i,val in enumerate(F[(ev,a)]):
                pw=-W+i
                gg=pw+2-l
                if gg%2 or gg<0:
                    ctot+=1
                    if val!=0:
                        bad2+=1; print("nonzero bad power",tau,mu,pw,val)
                    continue
                g=gg//2
                if g>3: continue
                if g==1 and l<1: continue
                if g>=4: continue
                if 2*g-2+l<=0 and g!=0: continue
                # Hodge
                Mv=peval(M_poly(g,mu),z) if not (g==0) else Fr(sum(mu))**(l-3)
                PiM=Fr(1)
                for j in mu: PiM*=Pi(j,tau)
                exp=sig**(g+l-1)*PiM*Mv/aut(mu)
                tot+=1
                if exp!=val:
                    bad+=1
                    print("MISMATCH",tau,mu,g,val,exp)
    print("tau",tau,"so far: compared",tot,"bad",bad,"other powers checked",ctot,"bad",bad2,flush=True)
