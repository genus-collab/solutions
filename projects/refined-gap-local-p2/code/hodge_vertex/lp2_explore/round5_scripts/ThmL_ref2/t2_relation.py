# Relation (4.1) of the report, evaluated with X_g(U) computed on the HODGE side (Faber engine),
# with no MV input and no partition-size limit. So this tests MV + P_k + dilaton calculus together,
# at |mu| far beyond the author's |mu| <= 20.
import sys,time
from hodge import *
from math import factorial

def Pi(k,tau):
    r=Fr(1)
    for a in range(1,k): r*=(k*tau+a)
    return r/factorial(k-1)

def X0(U):
    J=sum(U); r=len(U)
    # Gamma(J)/Gamma(J-r+3)
    return Fr(factorial(J-1),factorial(J-r+2))

def Xg(g,U,tau):
    if g==0: return X0(U)
    return X(g,U,tau)

def falling(a,k):
    r=Fr(1)
    for i in range(k): r*=(a-i)
    return r

def c_blk(sizeB,k):
    if sizeB%2==1:
        return 0 if k==0 else -1
    return 2 if k==0 else 1

def c_alpha(alpha,k0):
    # alpha tuple partition of k0
    r=Fr(1)
    for i in set(alpha):
        m=alpha.count(i)
        r*=Fr(-1,i)**m/factorial(m)
    if alpha==(k0,): r+=1
    return r

def set_partitions(lst):
    if not lst:
        yield []
        return
    first,rest=lst[0],lst[1:]
    for sp in set_partitions(rest):
        yield [[first]]+sp
        for i in range(len(sp)):
            yield sp[:i]+[[first]+sp[i]]+sp[i+1:]

def assignments(R,nb):
    # labelled assignments of positions in R to nb blocks
    if not R:
        yield [[] for _ in range(nb)]
        return
    for a in assignments(R[1:],nb):
        for i in range(nb):
            b=[list(x) for x in a]
            b[i].append(R[0])
            yield b

class Opts:
    fall=True; beta_sign=1; cwrong=False; evenwrong=False; pert=None

def phi(sB,T,g,k,tau,sigma,opts):
    U=tuple(sorted(T+sB))
    D=sum(U)-len(U)
    delta=D-2*g+2
    c=c_blk(len(sB),k)
    if opts.evenwrong and len(sB)%2==0 and k==0: c=1
    if c==0: return Fr(0)
    if opts.fall: f=falling(delta,k)
    else:
        f=Fr(1)
        for i in range(k): f*=(delta+i)
    if f==0: return Fr(0)
    w=Fr(c)
    for i in sB: w*=i*Pi(i,tau)
    w*=sigma**(len(sB)+g-1)*(opts.beta_sign*tau)**k/factorial(k)
    x=Xg(g,U,tau)
    if opts.pert and (g,U)==opts.pert: x+=1
    return w*f*x

def relation(k0,n,R,tau,opts=Opts(),gmaxavail=99):
    sigma=tau*(1+tau)
    total=Fr(0)
    for alpha in partitions(k0):
        ca=c_alpha(tuple(alpha),k0)
        if opts.cwrong and tuple(alpha)==(k0,): ca=Fr(1,k0)
        for sp in set_partitions(list(alpha)):
            nb=len(sp)
            for asg in assignments(list(R),nb):
                # per block enumerate (g,k,e)
                opts_blocks=[]
                dead=False
                for B,T in zip(sp,asg):
                    lst=[]
                    for e in range(0,n+1):
                        for g in range(0,e//2+1):
                            k=e-(2*g-2+len(B))
                            if k<0: continue
                            if g>gmaxavail: continue
                            v=phi(tuple(B),tuple(T),g,k,tau,sigma,opts)
                            if v!=0: lst.append((e,v))
                    opts_blocks.append(lst)
                # convolve on e
                poly={0:Fr(1)}
                for lst in opts_blocks:
                    Q={}
                    for e0,v0 in poly.items():
                        for e,v in lst:
                            if e0+e<=n:
                                Q[e0+e]=Q.get(e0+e,0)+v0*v
                    poly=Q
                total+=ca*poly.get(n,0)
    return total

if __name__=="__main__":
    pass
