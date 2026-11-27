# Theorem M_k (5.1): own implementation from the document's formula, with Hodge-side X values (MV-free).
from hodge import *
from itertools import product
from collections import Counter
from functools import lru_cache

def set_partitions(items):
    items=list(items)
    if not items:
        yield []
        return
    first=items[0]
    for sp in set_partitions(items[1:]):
        yield [[first]]+sp
        for i in range(len(sp)):
            yield sp[:i]+[[first]+sp[i]]+sp[i+1:]

def c_alpha(alpha,k0):
    c=Fr(1)
    for i,m in Counter(alpha).items():
        c*=Fr(-1,i)**m/factorial(m)
    if alpha==(k0,): c+=1
    return c

def X0(U):
    J=sum(U); r=len(U)
    if r==1: return Fr(1,J*(J+1))
    if r==2: return Fr(1,J)
    x=Fr(1)
    for t in range(1,r-2):
        x*=(J-t)
    return x

_Xc={}
def Xval(g,U,z):
    """X_g(U) at z; U sorted tuple"""
    if g==0: return X0(U)
    key=(g,U)
    if key not in _Xc:
        _Xc[key]=Xi_poly(g,U)
    return peval(_Xc[key],z)

def block_decomp(alpha,R):
    """Counter of configurations: tuple of (B sorted, R_B sorted) sorted -> number of (set partition, assignment) giving it"""
    out=Counter()
    pos=list(range(len(alpha)))
    r=len(R)
    for sp in set_partitions(pos):
        blocks=[tuple(sorted(alpha[i] for i in b)) for b in sp]
        for assign in product(range(len(blocks)),repeat=r):
            RB=[[] for _ in blocks]
            for part,b in zip(R,assign):
                RB[b].append(part)
            cfg=tuple(sorted((blocks[i],tuple(sorted(RB[i]))) for i in range(len(blocks))))
            out[cfg]+=1
    return out

class Params:
    def __init__(self,**kw):
        self.c_odd=Fr(-1); self.c_even0=Fr(2); self.c_even=Fr(1)
        self.beta_sign=1; self.rising=False; self.sig_shift=0; self.X_perturb=None
        self.c_top=Fr(0)
        self.__dict__.update(kw)

def cfun(b,k,P):
    if b%2==1:
        return Fr(0) if k==0 else P.c_odd
    return P.c_even0 if k==0 else P.c_even

def phi_table(B,RB,tau,Nmax,P,z,cache):
    """list: nb -> Sum_{g,k: 2g-2+|B|+k=nb} phi_B"""
    key=(B,RB)
    if key in cache: return cache[key]
    sig=tau*(1+tau)
    U=tuple(sorted(RB+B))
    b=len(B)
    D=sum(U)-len(U)
    pref=Fr(1)
    for i in B: pref*=i*Pi(i,tau)
    tab=[Fr(0)]*(Nmax+1)
    for nb in range(Nmax+1):
        tot=Fr(0)
        g=0
        while 2*g-2+b<=nb:
            k=nb-(2*g-2+b)
            cc=cfun(b,k,P)
            if cc:
                delta=D-2*g+2
                fal=rising(delta,k) if P.rising else falling(delta,k)
                if fal:
                    x=Xval(g,U,z)
                    if P.X_perturb is not None:
                        x=P.X_perturb(g,U,x)
                    if x:
                        tot+=cc*pref*sig**(b+g-1+P.sig_shift)*(P.beta_sign*tau)**k/factorial(k)*fal*x
            g+=1
        tab[nb]=tot
    cache[key]=tab
    return tab

def relation(k0,R,tau,Nmax,P=None):
    """returns list n-> LHS of (5.1), n<=Nmax"""
    if P is None: P=Params()
    z=z_of_tau(tau)
    R=tuple(sorted(R))
    cache={}
    tot=[Fr(0)]*(Nmax+1)
    for alpha in partitions(k0):
        ca=c_alpha(alpha,k0)
        if alpha==(k0,): ca+=P.c_top
        if not ca: continue
        for cfg,mult in block_decomp(alpha,R).items():
            pr=[Fr(0)]*(Nmax+1); pr[0]=Fr(1)
            for (B,RB) in cfg:
                t=phi_table(B,RB,tau,Nmax,P,z,cache)
                nw=[Fr(0)]*(Nmax+1)
                for i,x in enumerate(pr):
                    if x:
                        for j in range(Nmax+1-i):
                            if t[j]:
                                nw[i+j]+=x*t[j]
                pr=nw
            for n in range(Nmax+1):
                tot[n]+=ca*mult*pr[n]
    return tot
