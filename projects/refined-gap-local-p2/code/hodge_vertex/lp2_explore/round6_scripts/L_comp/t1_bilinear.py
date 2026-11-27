# Test (i): (4.2) on Z_0, Z_1 as polynomial identities in A and p'.
# Z_n = Sum_nu A^{f2(nu)+n|nu|} s_nu(p') s_nu(p); t_k = p_k/k; H_k Hirota in t_k.
import sys, itertools, time, random
from common import *
from collections import defaultdict
PR=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]
def code(lam):
    c=1
    for x in lam: c*=PR[x-1]
    return c

def build(Wmax, nshift=(0,1)):
    """Z[n][N][mu] = dict (Aexp, code lam) -> Fr"""
    Z={n:{} for n in nshift}
    for n in nshift:
        for N in range(Wmax+1):
            Z[n][N]={}
            pars=list(partitions(N))
            for mu in pars:
                d=defaultdict(Fr)
                for lam in pars:
                    cc=code(lam); zz=zee(mu)*zee(lam)
                    for nu in pars:
                        cf=chi(nu,mu)*chi(nu,lam)
                        if cf:
                            d[(contents_sum(nu)+n*N,cc)]+=Fr(cf,zz)
                Z[n][N][mu]={k:v for k,v in d.items() if v}
    return Z

def deriv(mu, part):
    """d/dt_part p_mu = part * d/dp_part p_mu -> (coef, mu')"""
    m=mu.count(part)
    if m==0: return None
    l=list(mu); l.remove(part)
    return (part*m, tuple(l))

def derivs(mu, parts):
    c=1
    for p in parts:
        r=deriv(mu,p)
        if r is None: return None
        c*=r[0]; mu=r[1]
    return c,mu

def merge(mu1,mu2):
    return tuple(sorted(mu1+mu2,reverse=True))

def c_alpha(alpha,k,perturb=None):
    c=Fr(1)
    for i,m in Counter(alpha).items():
        c*=Fr(-1,i)**m/factorial(m)
    if alpha==(k,):
        c+=1
    if perturb and alpha==(k,):
        c+=perturb
    return c

def hirota_poly(Zf,Zg,alpha,d):
    """H^alpha Zf.Zg, coefficient of p-monomials of weight d. returns dict mono-> dict (Ae,code)->Fr"""
    k=sum(alpha); W=d+k
    out=defaultdict(lambda: defaultdict(Fr))
    pos=list(range(len(alpha)))
    # group subsets by (S parts) to save work
    splits=Counter()
    for mask in range(1<<len(alpha)):
        S=tuple(sorted(alpha[i] for i in pos if mask>>i&1))
        T=tuple(sorted(alpha[i] for i in pos if not mask>>i&1))
        splits[(S,T)]+=(-1)**len(T)
    for (S,T),sg in splits.items():
        if sg==0: continue
        for a in range(sum(S),W-sum(T)+1):
            b=W-a
            for mu1,P1 in Zf[a].items():
                r1=derivs(mu1,S)
                if r1 is None: continue
                for mu2,P2 in Zg[b].items():
                    r2=derivs(mu2,T)
                    if r2 is None: continue
                    mono=merge(r1[1],r2[1])
                    cf=sg*r1[0]*r2[0]
                    tgt=out[mono]
                    for (e1,c1),v1 in P1.items():
                        for (e2,c2),v2 in P2.items():
                            tgt[(e1+e2,c1*c2)]+=cf*v1*v2
    return out

def Pk_value(Zf,Zg,k,d,perturb=None):
    tot=defaultdict(lambda: defaultdict(Fr))
    for alpha in partitions(k):
        c=c_alpha(alpha,k,perturb)
        if c==0: continue
        res=hirota_poly(Zf,Zg,alpha,d)
        for m,P in res.items():
            t=tot[m]
            for key,v in P.items():
                t[key]+=c*v
    return tot

def nz(tot):
    n=0;terms=0
    for m,P in tot.items():
        nzt=sum(1 for v in P.values() if v)
        terms+=nzt
        if nzt: n+=1
    return n,terms

if __name__=="__main__":
    Wmax=int(sys.argv[1]) if len(sys.argv)>1 else 9
    t0=time.time()
    Z=build(Wmax)
    print("built Z up to weight",Wmax,time.time()-t0,flush=True)
    # Check examples of P_k vs (4.3)
    print("c_alpha 2P2:",{a:2*c_alpha(a,2) for a in partitions(2)})
    print("c_alpha 6P3:",{a:6*c_alpha(a,3) for a in partitions(3)})
    print("c_alpha 24P4:",{a:24*c_alpha(a,4) for a in partitions(4)})
    print("c_alpha P1:",{a:c_alpha(a,1) for a in partitions(1)})
    summary=[]
    for k in range(2,7):
        for d in range(0,Wmax-k+1):
            t1=time.time()
            tot=Pk_value(Z[0],Z[1],k,d)
            nm,nt=nz(tot)
            # number of p-monomials examined
            print(f"k={k} d={d}: monomials touched={len(tot)} nonzero monomials={nm} nonzero terms={nt} time={time.time()-t1:.1f}",flush=True)
            summary.append((k,d,nm))
            if d in (0,1,2):
                # controls: swapped order, perturbed c_(k), wrong charge
                sw=Pk_value(Z[1],Z[0],k,d)
                pt=Pk_value(Z[0],Z[1],k,d,perturb=Fr(1,7))
                print(f"   controls: swap Z1.Z0 nonzero={nz(sw)}  perturbed c_(k) nonzero={nz(pt)}",flush=True)
    print("TOTAL nonzero identities:",sum(s[2] for s in summary),"of cases",len(summary))
