# s.6: term-by-term mechanical run of the induction with REAL values, plus target response.
# For T in range at genus G, each part k0, R=T-k0, n=2G: enumerate every term of (5.1); check
#  (a) Count identity Sum_B(delta_B-k_B)=D+1-2G,
#  (b) each term's product has a factor classified per s.6 and vanishes, except target (which vanishes by Theorem L),
#  (c) response: replacing X_G(T) by X_G(T)+1 changes the LHS by exactly -(k0-1) Pi_{k0} sigma^G tau delta_G(T).
from mk import *
from t4_lists import multisets
import sys,time
Gmax=int(sys.argv[1]); maxr=int(sys.argv[2])
taus=[Fr(2,7),Fr(-5,3)]
def enum_budget(cfg,n):
    """yield tuple of (g,k) per block with Sum n_B=n."""
    blocks=list(cfg)
    def rec(i,left):
        if i==len(blocks):
            if left==0: yield ()
            return
        B,RB=blocks[i]
        b=len(B)
        g=0
        while 2*g-2+b<=left:
            kmax=left-(2*g-2+b)
            for k in range(0,kmax+1):
                if b%2==1 and k==0: continue
                nb=2*g-2+b+k
                for rest in rec(i+1,left-nb):
                    yield ((g,k),)+rest
            g+=1
    return rec(0,n)

stats=Counter()
for tau in taus:
    z=z_of_tau(tau); sig=tau*(1+tau)
    for G in range(2,Gmax+1):
        for D in range(1,2*G-2):
            for T in multisets(D,maxr=maxr):
                for k0 in sorted(set(T)):
                    if k0>7: stats['skipped k0>7']+=1; continue
                    R=list(T); R.remove(k0); R=tuple(sorted(R))
                    n=2*G
                    D_R=sum(R)-len(R)
                    LHS=Fr(0); nz_terms=0; targets=0
                    for alpha in partitions(k0):
                        ca=c_alpha(alpha,k0)
                        for cfg,mult in block_decomp(alpha,R).items():
                            for gk in enum_budget(cfg,n):
                                stats['terms']+=1
                                su=0; deficient=[]; pr=ca*mult;
                                for (B,RB),(g,k) in zip(cfg,gk):
                                    U=tuple(sorted(RB+B)); DU=sum(U)-len(U); delta=DU-2*g+2
                                    su+=delta-k
                                    if delta<k: deficient.append((B,RB,g,k,U,delta))
                                    cc=Fr(-1) if len(B)%2 else (Fr(2) if k==0 else Fr(1))
                                    pref=Fr(1)
                                    for i in B: pref*=i*Pi(i,tau)
                                    pr*=cc*pref*sig**(len(B)+g-1)*tau**k/factorial(k)*falling(delta,k)
                                    pr*= Xval(g,U,z) if pr else 1
                                assert su==D+1-2*G,(su,D,G)
                                assert deficient
                                # classification
                                cls=None
                                for (B,RB,g,k,U,delta) in deficient:
                                    if delta>=0: cls='falling'; break
                                    assert g>=2
                                    if g<G: cls='H1'; assert not Xi_poly(g,U); break
                                    # g==G
                                    assert len(B)+k==2
                                    if len(alpha)>=2:
                                        assert sum(U)-len(U)<D
                                        assert not Xi_poly(g,U)
                                        cls='H2'; break
                                    assert U==tuple(sorted(T)) and alpha==(k0,) and k==1
                                    cls='target'; targets+=1
                                    break
                                stats[cls]+=1
                                if pr: nz_terms+=1
                                LHS+=pr
                    stats['cases']+=1
                    assert targets==1,(T,k0,targets)
                    assert nz_terms==0,(T,k0,nz_terms)
                    assert LHS==0
                    # response: analytic coefficient of X_G(T): count terms containing X_G(T) as a factor with its coefficient
                    coef=Fr(0)
                    for alpha in partitions(k0):
                        ca=c_alpha(alpha,k0)
                        for cfg,mult in block_decomp(alpha,R).items():
                            for gk in enum_budget(cfg,n):
                                pos=[i for i,((B,RB),(g,k)) in enumerate(zip(cfg,gk)) if g==G and tuple(sorted(RB+B))==tuple(sorted(T))]
                                if not pos: continue
                                pr=ca*mult
                                for i,((B,RB),(g,k)) in enumerate(zip(cfg,gk)):
                                    U=tuple(sorted(RB+B)); DU=sum(U)-len(U); delta=DU-2*g+2
                                    cc=Fr(-1) if len(B)%2 else (Fr(2) if k==0 else Fr(1))
                                    pref=Fr(1)
                                    for j in B: pref*=j*Pi(j,tau)
                                    pr*=cc*pref*sig**(len(B)+g-1)*tau**k/factorial(k)*falling(delta,k)
                                    if i not in pos[:1]: pr*=Xval(g,U,z)
                                coef+=pr
                    exp=-(k0-1)*Pi(k0,tau)*sig**G*tau*(D+2-2*G)
                    assert coef==exp,(T,k0,coef,exp)
                    stats['target coefficient ok']+=1
        print("tau",tau,"G",G,dict(stats),flush=True)
print("DONE",dict(stats))
