# Small exact checks of displayed identities in a3_uconv.tex (e3), (f1), (f2).
from fractions import Fraction as Fr
from math import comb, factorial
N=12  # truncation order in u
def mul(a,b):
    c=[Fr(0)]*N
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j<N: c[i+j]+=x*y
    return c
def inv(a):
    # a[0]!=0
    b=[Fr(0)]*N; b[0]=1/a[0]
    for n in range(1,N):
        s=sum(a[k]*b[n-k] for k in range(1,n+1))
        b[n]=-s/a[0]
    return b
def pw(a,k):
    r=[Fr(1)]+[Fr(0)]*(N-1)
    for _ in range(k): r=mul(r,a)
    return r
# (f1): kappa(0)=1, c=1-1/kappa ; sum_{m1} binom(Nn+m1-1,m1) c^m1 == kappa^Nn
kappa=[Fr(1),Fr(3,7),Fr(-2,5),Fr(1,3)]+[Fr(0)]*(N-4)
c=[ -x for x in inv(kappa)]; c[0]+=1
for Nn in [1,2,5]:
    s=[Fr(0)]*N
    for m1 in range(N+1):
        t=pw(c,m1); b=comb(Nn+m1-1,m1)
        s=[s[i]+b*t[i] for i in range(N)]
    assert s==pw(kappa,Nn), Nn
print("f1 binomial ok")
# (e3): sum_j delta^j/j! * d^j/dv^j (v^{2-2g}) == (v+delta)^{2-2g}, with v=u*unit, delta in u^2
# work with w = v/u (unit) and x=delta/v in u Q[[u]]: check (1+x)^{-k} = sum_j binom(-k,j) x^j
def binom_gen(a,j):
    r=Fr(1)
    for i in range(j): r*= Fr(a-i, 1)
    return r/factorial(j)
x=[Fr(0),Fr(2,3),Fr(-1,4),Fr(5,9)]+[Fr(0)]*(N-4)
for g in [2,3,4]:
    k=2*g-2
    s=[Fr(0)]*N
    for j in range(N+1):
        t=pw(x,j); b=binom_gen(2-2*g,j)
        s=[s[i]+b*t[i] for i in range(N)]
    onepx=[Fr(1)]+x[1:]
    assert s==inv(pw(onepx,k)), g
print("e3 binomial ok")
# (f2): for c0!=0 the u^o coefficient sum_{m2 in S} c0^{m-m2} a_{m2} binom(N_{m2}+m-m2-1, N_{m2}-1) is nonzero for large m
c0=Fr(-1,2); S={0:Fr(3),1:Fr(-6),2:Fr(1,2)}; base=1  # N_{m2}=base+m2
vals=[sum(c0**(m-m2)*a*comb(base+m2+m-m2-1,base+m2-1) for m2,a in S.items()) for m in range(2,40)]
print("f2 first values nonzero:", all(v!=0 for v in vals[5:]))
