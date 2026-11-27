from mk import *
import sys,time
Nmax=int(sys.argv[1]); kmax=int(sys.argv[2])
taus=[Fr(2,7),Fr(-5,3),Fr(13,9),Fr(-1,2),Fr(-1,3),Fr(-2,3),Fr(-1,4),Fr(-1,5),Fr(1,1)]
Rs=[(),(2,),(3,),(4,),(5,),(2,2),(2,3),(3,3),(2,4),(2,2,2),(2,2,3)]
def perturbX1(g,U,x):
    return x+Fr(1,1000) if g==1 else x
def perturbX0(g,U,x):
    return x*(1+Fr(1,1000)) if (g==0 and len(U)==1) else x
def perturbX2(g,U,x):
    return x+Fr(1,1000) if (g==2) else x
controls={
 'c_odd=+1':Params(c_odd=Fr(1)),
 'c_even0=1':Params(c_even0=Fr(1)),
 'c_even=2':Params(c_even=Fr(2)),
 'beta->-beta':Params(beta_sign=-1),
 'rising':Params(rising=True),
 'sigma shift':Params(sig_shift=1),
 'X1 perturbed':Params(X_perturb=perturbX1),
 'X0(r=1) perturbed':Params(X_perturb=perturbX0),
 'X2 perturbed':Params(X_perturb=perturbX2),
 'c_(k) perturbed':Params(c_top=Fr(1,7)),
}
tot=0;fail=0
ctlfail=Counter();ctltot=Counter()
for tau in taus:
    t0=time.time()
    for k0 in range(2,kmax+1):
        for R in Rs:
            # budget: keep runtime modest for large profiles
            N=Nmax
            if len(R)>=2 and k0>=4: N=min(Nmax,6)
            if len(R)>=3 and k0>=3: N=min(Nmax,6)
            v=relation(k0,R,tau,N)
            for n,x in enumerate(v):
                tot+=1
                if x!=0:
                    fail+=1
                    print("FAIL (5.1)",tau,k0,R,n,x,flush=True)
            if tau in (Fr(2,7),Fr(-1,2)):
                for name,P in controls.items():
                    w=relation(k0,R,tau,N,P)
                    for n,x in enumerate(w):
                        ctltot[name]+=1
                        if x!=0: ctlfail[name]+=1
    print(f"tau={tau}: cumulative relations {tot}, fails {fail}, t={time.time()-t0:.1f}",flush=True)
print("controls: relations broken / total (2 tau values, includes relations with empty content)")
for k in controls: print(" ",k,ctlfail[k],"/",ctltot[k])
print("TOTAL",tot,"fails",fail)
