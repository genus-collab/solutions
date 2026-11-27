# (f1)/(f2)/(e3) formal checks with truncated power series over Q (exact).
from fractions import Fraction as F
from math import comb, factorial
import random
ORD = 14
def mul(a,b):
    r=[F(0)]*ORD
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>=ORD: break
            r[i+j]+=x*y
    return r
def add(a,b): return [x+y for x,y in zip(a,b)]
def scal(c,a): return [c*x for x in a]
def inv(a):
    r=[F(0)]*ORD; r[0]=1/a[0]
    for n in range(1,ORD):
        r[n] = -sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
    return r
def powi(a,n):
    r=[F(1)]+[F(0)]*(ORD-1)
    if n<0: a=inv(a); n=-n
    for _ in range(n): r=mul(r,a)
    return r
random.seed(1)
def rnd(c0=None):
    a=[F(random.randint(-3,3),random.randint(1,3)) for _ in range(ORD)]
    if c0 is not None: a[0]=F(c0)
    return a
one=[F(1)]+[F(0)]*(ORD-1)

# (f1): kappa(0)=1, c=1-1/kappa; A_{m2} arbitrary series (stand for <X y^m2>_h / m2!), m2<=Dmax, N_{m2}=2h-2+n+m2
for (h,n) in [(0,3),(1,1),(2,0),(2,2),(3,1)]:
    Dmax = 3*h-3+n
    kappa = rnd(1)
    c = add(one, scal(-1,inv(kappa)))
    assert c[0]==0
    A = [rnd() for _ in range(Dmax+1)]
    # left side: sum_m T_m, T_m = sum_{m1+m2=m} c^{m1} binom(N_{m2}+m1-1,m1) A_{m2}; c^{m1}=O(u^{m1}) so m1<ORD suffices
    left=[F(0)]*ORD
    cp=[one]
    for m1 in range(1,ORD): cp.append(mul(cp[-1],c))
    for m2 in range(Dmax+1):
        N = 2*h-2+n+m2
        for m1 in range(ORD):
            left = add(left, scal(comb(N+m1-1,m1), mul(cp[m1],A[m2])))
    right=[F(0)]*ORD
    for m2 in range(Dmax+1):
        N = 2*h-2+n+m2
        right = add(right, mul(powi(kappa,N),A[m2]))
    assert left==right,(h,n)
print("(f1) resummation identity: OK")

# (f2): kappa(0)=2 -> c0=1/2 ; coefficient of u^o of T_m is a nonzero polynomial in m times c0^m
h,n=2,1; Dmax=3*h-3+n
kappa=rnd(2); c=add(one,scal(-1,inv(kappa))); c0=c[0]; assert c0!=0
A=[rnd(1) for _ in range(Dmax+1)]   # all of order 0 -> o=0, S = all m2
vals=[]
for m in range(Dmax, Dmax+12):
    t = sum(c0**(m-m2)*comb(2*h-2+n+m2+m-m2-1, m-m2)*A[m2][0] for m2 in range(Dmax+1))
    # direct: [u^0] T_m
    vals.append(t)
assert all(v!=0 for v in vals[3:]), vals
print("(f2) [u^o]T_m nonzero for large m (sample): OK", [float(v) for v in vals[:4]])

# (e3): sum_j binom(2-2g,j) (delta/v)^j v^{2-2g} = (v+delta)^{2-2g}, v=u*unit, delta in u^2 Q[[u]]
for g in (2,3,4):
    vunit=rnd(1); dl=[F(0),F(0)]+rnd()[:ORD-2]
    x = mul([F(0)]+dl[2:]+[F(0)], inv(vunit))[:ORD]      # delta/v = u * (delta/u^2)/unit
    # delta/v: delta = u^2 * D, v = u*U  -> delta/v = u*D/U
    D = dl[2:]+[F(0),F(0)]
    x = [F(0)]+mul(D,inv(vunit))[:ORD-1]
    lhs=[F(0)]*ORD; xp=one
    def gbin(a,j):
        r=F(1)
        for i in range(j): r*=F(a-i,i+1)
        return r
    for j in range(ORD):
        lhs=add(lhs,scal(gbin(2-2*g,j),xp)); xp=mul(xp,x)
    rhs = powi(add(one,x),2-2*g)
    assert lhs==rhs
print("(e3) binomial resummation: OK")

# redundancy of delta in u^2: v-exponent of Pi_0(z^m H^b) is m - floor(m/5) >= 2 for m>=2
assert all(m-m//5>=2 for m in range(2,200))
print("Pi_0(z^2 H_+) in v^2 Q[v]: OK")
