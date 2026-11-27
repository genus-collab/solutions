# Independent sanity check of prop:set-cyclicLR (i): kappa^LR = (1, L^5/(C0^2 C1^3), 1, C0/C1, C0/C1)
# exact arithmetic, truncated power series in q, polynomials in the formal log l.
from fractions import Fraction as Fr
from math import factorial
N = 7   # work mod q^N
W = 5   # need w^0..w^4

def smul(a,b):
    c=[Fr(0)]*N
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>=N: break
            c[i+j]+=x*y
    return c
def sadd(a,b): return [x+y for x,y in zip(a,b)]
def ssub(a,b): return [x-y for x,y in zip(a,b)]
def sinv(a):
    b=[Fr(0)]*N; b[0]=1/a[0]
    for n in range(1,N):
        s=Fr(0)
        for k in range(1,n+1): s+=a[k]*b[n-k]
        b[n]=-s/a[0]
    return b
def sD(a): return [n*a[n] for n in range(N)]
one=[Fr(1)]+[Fr(0)]*(N-1)

# polynomials in w truncated at degree W (list of Fractions)
def pmul(a,b):
    c=[Fr(0)]*W
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<W: c[i+j]+=x*y
    return c
def pinv(a):
    b=[Fr(0)]*W; b[0]=1/a[0]
    for n in range(1,W):
        s=Fr(0)
        for k in range(1,n+1): s+=a[k]*b[n-k]
        b[n]=-s/a[0]
    return b
# P(w,q) = sum_d q^d prod_{k=1}^{5d}(5w+k)/prod_{k=1}^d((w+k)^5-w^5)
ptil=[[Fr(0)]*N for _ in range(W)]   # ptil[m][d]
for d in range(N):
    num=[Fr(1)]+[Fr(0)]*(W-1)
    for k in range(1,5*d+1):
        num=pmul(num,[Fr(k),Fr(5)]+[Fr(0)]*(W-2))
    den=[Fr(1)]+[Fr(0)]*(W-1)
    for k in range(1,d+1):
        # (w+k)^5 - w^5 = 5k w^4 + 10k^2 w^3 + 10 k^3 w^2 + 5 k^4 w + k^5
        den=pmul(den,[Fr(k**5),Fr(5*k**4),Fr(10*k**3),Fr(10*k**2),Fr(5*k)])
    t=pmul(num,pinv(den))
    for m in range(W): ptil[m][d]=t[m]
C0=ptil[0]
assert C0[:4]==[1,120,113400,168168000]
# I1tilde check
I1t=[Fr(0)]*N
for d in range(1,N):
    I1t[d]=5*Fr(factorial(5*d),factorial(d)**5)*sum(Fr(1,m) for m in range(d+1,5*d+1))
assert I1t==ptil[1]
# elements of Q[[q]][l]: dict j -> series
def lD(f):
    g={}
    for j,a in f.items():
        g[j]=sadd(g.get(j,[Fr(0)]*N),sD(a))
        if j>=1: g[j-1]=sadd(g.get(j-1,[Fr(0)]*N),[j*x for x in a])
    return g
def lscal(s,f): return {j:smul(s,a) for j,a in f.items()}
def frak_p(k): return {j:[x/factorial(j) for x in ptil[k-j]] for j in range(k+1)}
C0inv=sinv(C0)
T=lscal(C0inv,frak_p(1))
DT=lD(T)
assert all(x==0 for x in DT.get(1,[0]))
C1=DT[0]
C1inv=sinv(C1)
def dT(f): return lscal(C1inv,lD(f))
def iszero(a): return all(x==0 for x in a)
def clean(f): return {j:a for j,a in f.items() if not iszero(a)}
# L^5 = 1/(1-5^5 q)
L5=[Fr(3125)**n for n in range(N)]
p2=lscal(C0inv,frak_p(2)); p3=lscal(C0inv,frak_p(3)); p4=lscal(C0inv,frak_p(4))
k1=clean(dT(dT(p2)))
assert list(k1.keys())==[0]
k1=k1[0]
target=smul(L5,sinv(smul(smul(C0,C0),smul(C1,smul(C1,C1)))))
print("kappa1 == L^5/(C0^2 C1^3) mod q^%d:"%N, k1==target)
k1inv=sinv(k1)
k2=clean(dT(lscal(k1inv,dT(dT(p3)))))
print("kappa2 == 1:", list(k2.keys())==[0] and k2[0]==one)
k3=clean(dT(dT(lscal(k1inv,dT(dT(p4))))))
print("kappa3 == C0/C1:", list(k3.keys())==[0] and k3[0]==smul(C0,C1inv))
# note: each D loses no q-order (q d/dq), so truncation mod q^N is exact.
# PF4 p_k = 0 (k<=3), PF4 p_4 = 1
def PF4(f):
    d1=lD(f); d2=lD(d1); d3=lD(d2); d4=lD(d3)
    out={}
    def acc(c,g,shift):
        for j,a in g.items():
            b=[Fr(0)]*N
            for n in range(N-shift): b[n+shift]=c*a[n]
            out[j]=sadd(out.get(j,[Fr(0)]*N),b)
    acc(1,d4,0); acc(-3125,d4,1); acc(-6250,d3,1); acc(-4375,d2,1); acc(-1250,d1,1); acc(-120,f,1)
    return clean(out)
for k in range(5):
    print("PF4 p_%d ="%k, {j:a[:3] for j,a in PF4(frak_p(k)).items()})
