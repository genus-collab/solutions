from fractions import Fraction as Fr
# --- Part 1: PF4 in u-form acting on Laurent polys in u with log: represent f = sum_n (a_n + b_n*ell) u^n
# D = -(1-u) d/du ; D ell = 1 - 1/u
from collections import defaultdict
def D(f):
    # f: dict (n,k)-> coeff, meaning u^n ell^k, k in {0,1}
    g = defaultdict(Fr)
    for (n,k),c in f.items():
        # d/du (u^n ell^k) = n u^(n-1) ell^k + k u^n ell^(k-1) * (1/u) [since theta ell =1 => d ell/du = 1/u]
        # D = -(1-u) d/du
        if n!=0:
            g[(n-1,k)] += -c*n
            g[(n,k)]   += c*n
        if k==1:
            g[(n-1,0)] += -c
            g[(n,0)]   += c
    return {key:v for key,v in g.items() if v!=0}
def add(f,g,cf=1,cg=1):
    h=defaultdict(Fr)
    for k,v in f.items(): h[k]+=cf*v
    for k,v in g.items(): h[k]+=cg*v
    return {k:v for k,v in h.items() if v!=0}
def mulu(f,m=1):
    return {(n+m,k):c for (n,k),c in f.items()}
def PF4(f):
    d1=D(f); d2=D(d1); d3=D(d2); d4=D(d3)
    inner = add(add(add({k:2*v for k,v in d3.items()}, {k:Fr(7,5)*v for k,v in d2.items()}), {k:Fr(2,5)*v for k,v in d1.items()}), {k:Fr(24,625)*v for k,v in f.items()})
    # u D^4 - (1-u) inner
    res = add(mulu(d4), inner, 1, -1)
    res = add(res, mulu(inner))
    return res
# P_j(n): PF4 u^n = sum_j P_j(n) u^{n+j}
def P(j,n):
    r=PF4({(n,0):Fr(1)})
    return r.get((n+j,0),Fr(0))
for n in range(-3,8):
    x=Fr(n)
    assert P(-3,n)==x*(x-1)**2*(x-2), (n,P(-3,n))
    assert P(-2,n)==-x*(x-1)*(4*x*x-6*x+Fr(17,5)), (n,P(-2,n))
print("P_-3,P_-2 ok; P_-1(0)=",P(-1,0))
r=PF4({(0,1):Fr(1)})
print("PF4(ell) =",sorted(r.items()))
# solutions f_k = u^k + O(u^3), k=0,1,2 up to order N
N=12
def solve(b0,b1,b2):
    b=[Fr(b0),Fr(b1),Fr(b2)]
    # coefficient of u^M in PF4 g: sum_j P_j(M-j) b_{M-j}; solve b_{M+3}
    for M in range(0,N-3):
        s=Fr(0)
        for j in (-2,-1,0,1):
            n=M-j
            if 0<=n<len(b): s+=P(j,n)*b[n]
        b.append(-s/P(-3,M+3))
    return b
f0=solve(1,0,0); f1=solve(0,1,0); f2=solve(0,0,1)
print("f0",f0[:6]); print("f1",f1[:6]); print("f2",f2[:6])
w1=[f1[i]+Fr(7,10)*f2[i] for i in range(N)]
print("varpi1",w1[:6])
# tc = varpi1/f0 in frame (0,0)
def inv(a):
    b=[1/a[0]]
    for n in range(1,len(a)):
        b.append(-sum(a[k]*b[n-k] for k in range(1,n+1))/a[0])
    return b
def mul(a,b):
    n=min(len(a),len(b)); return [sum(a[k]*b[i-k] for k in range(i+1)) for i in range(n)]
tc=mul(w1,inv(f0))
print("tc",tc[:6])
# obstruction r(g): coefficients u^-3,u^-2,u^-1 of PF4(g ell) for g=f0,f1,f2
for name,g in (("f0",f0),("f1",f1),("f2",f2)):
    r=PF4({(n,1):g[n] for n in range(6) if g[n]!=0})
    print(name,[r.get((m,0),Fr(0)) for m in (-3,-2,-1)], "log part low:", [r.get((m,1),Fr(0)) for m in range(-3,2)])
