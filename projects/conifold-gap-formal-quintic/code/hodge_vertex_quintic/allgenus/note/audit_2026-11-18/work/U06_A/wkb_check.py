# Check (lem:gau-wkb (i), thm:gau-R0 (iii)): in Q[[q]] with L=(1-5^5 q)^(-1/5),
#  (a) coefficient of z^1 of PF5^{(0)} applied to L*g is 5*D g;
#  (b) A = sum_n L^{4n+1} ph_n(L^{-5}) z^n (ph_n from the HPF recursion) satisfies PF5^{(0)} A = 0.
from fractions import Fraction as Fr
import importlib.util, sys, io, contextlib
spec=importlib.util.spec_from_file_location("hpf","hpf_check.py")
m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(m)
K=9   # q-adic precision
def smul(a,b):
    r=[Fr(0)]*K
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>=K: break
            r[i+j]+=x*y
    return r
def sadd(a,b,c=Fr(1)): return [x+c*y for x,y in zip(a,b)]
def D(a): return [i*x for i,x in enumerate(a)]
def binom_series(alpha,c):
    # (1+c q)^alpha
    r=[Fr(1)]
    for n in range(1,K):
        r.append(r[-1]*(alpha-(n-1))/n*c)
    return r
L=binom_series(Fr(-1,5),Fr(-5**5))
Linv=binom_series(Fr(1,5),Fr(-5**5))
assert smul(L,Linv)==[Fr(1)]+[Fr(0)]*(K-1)
def spow(a,n):
    r=[Fr(1)]+[Fr(0)]*(K-1)
    base=a if n>=0 else None
    if n<0:
        base=Linv if a is L else None
        n=-n
    for _ in range(n): r=smul(r,base)
    return r
qser=[Fr(0),Fr(1)]+[Fr(0)]*(K-2)
# series in z: list of q-series
def Nop(A):   # (zD+L)A
    r=[[Fr(0)]*K for _ in range(len(A)+1)]
    for n,a in enumerate(A):
        r[n]=sadd(r[n],smul(L,a))
        r[n+1]=sadd(r[n+1],D(a))
    return r
def zmul(A,c):  # c*z*A
    return [[Fr(0)]*K]+[[c*x for x in a] for a in A]
def zadd(A,B,c=Fr(1)):
    n=max(len(A),len(B)); r=[]
    for i in range(n):
        a=A[i] if i<len(A) else [Fr(0)]*K
        b=B[i] if i<len(B) else [Fr(0)]*K
        r.append(sadd(a,b,c))
    return r
def PF5_0(A):
    X=A
    for _ in range(5): X=Nop(X)
    Y=A
    for k in range(1,6):
        Y=zadd([[5*x for x in a] for a in Nop(Y)], zmul(Y,Fr(k)))
    Y=[smul(qser,a) for a in Y]
    return zadd(zadd(X,A,Fr(-1)),Y,Fr(-1))
# (a)
g=[Fr(3),Fr(-2),Fr(7,5),Fr(1,3),Fr(-4),Fr(2),Fr(1),Fr(5),Fr(-1)][:K]
out=PF5_0([smul(L,g)])
print("z^0 coefficient zero:", all(x==0 for x in out[0]))
print("z^1 coefficient == 5 D g:", out[1]==[5*x for x in D(g)])
# (b)
N=8
for free in ({}, {1:Fr(7,3)}):
    ph,_=m.hpf(N,free)
    A=[]
    for n in range(N+1):
        a=[Fr(0)]*K
        for k,c in ph[n].items():
            a=sadd(a,spow(L,4*n+1-5*k),c)
        A.append(a)
    out=PF5_0(A)
    ok=all(all(x==0 for x in out[mm]) for mm in range(N+2))
    print("free=",dict(free)," PF5^(0)A = 0 through z^%d mod q^%d:"%(N+1,K), ok)
