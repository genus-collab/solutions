# all members of (4.1): coefficient of eps^beta in Sum_j Sch_j(2eps)Sch_j(-H~)exp(eps.H) - exp(-eps.H) annihilates Z0.Z1
from t1_bilinear import *
from itertools import product
W=9
Z=build(W)
def sch_poly(j):
    # Sch_j(-H~) as dict alpha->coef
    return {al:(lambda c:c)(c_alpha(al,-1)) for al in partitions(j)}
def betas(maxw):
    for n in range(1,maxw+1):
        for b in partitions(n):
            yield b
def member(beta):
    """dict alpha->coef"""
    out=defaultdict(Fr)
    cb=Counter(beta); vals=sorted(cb)
    fb=Fr(1)
    for v in vals: fb/=factorial(cb[v])
    out[tuple(sorted(beta,reverse=True))]-=(-1)**len(beta)*fb
    for ev in product(*[range(cb[v]+1) for v in vals]):
        gam=[];rest=[]
        cg=Fr(1)
        for v,e in zip(vals,ev):
            gam+=[v]*e; rest+=[v]*(cb[v]-e)
            cg*=Fr(2)**e/factorial(e)/factorial(cb[v]-e)
        j=sum(gam)
        for al,c in sch_poly(j).items():
            key=tuple(sorted(list(al)+rest,reverse=True))
            out[key]+=cg*c
    return {k:v for k,v in out.items() if v}
tot=0;bad=0
for beta in betas(5):
    m=member(beta)
    if not m: continue
    wts={sum(k) for k in m}
    # inhomogeneous: only homogeneous components by weight separately? operator total weight varies; each weight part acts separately because of p' grading
    for w in sorted(wts):
        mw={k:v for k,v in m.items() if sum(k)==w}
        for d in range(0,W-w+1):
            T=defaultdict(lambda: defaultdict(Fr))
            for al,c in mw.items():
                res=hirota_poly(Z[0],Z[1],al,d)
                for mono,P in res.items():
                    t=T[mono]
                    for key,v in P.items(): t[key]+=c*v
            n=nz(T); tot+=1
            if n[0]: bad+=1; print("NONZERO",beta,w,d,n,flush=True)
print("eps members tested",tot,"bad",bad)
print("example member beta=(1,1):",member((1,1)))
print("example member beta=(2,):",member((2,)))
