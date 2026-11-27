# Independent sanity check of prop:set-cyclicLR (i), cor:set-cyclicPF5 (a), def:set-genus0
# Exact arithmetic with fractions; series in q truncated at order N, polynomials in l (formal log q).
from fractions import Fraction as Fr
from math import factorial
N = 6

def s_zero(): return [Fr(0)]*(N+1)
def s_one():
    a = s_zero(); a[0] = Fr(1); return a
def s_add(a,b): return [x+y for x,y in zip(a,b)]
def s_sub(a,b): return [x-y for x,y in zip(a,b)]
def s_scal(c,a): return [c*x for x in a]
def s_mul(a,b):
    c = s_zero()
    for i in range(N+1):
        if a[i]==0: continue
        for j in range(N+1-i):
            c[i+j] += a[i]*b[j]
    return c
def s_inv(a):
    assert a[0]!=0
    b = s_zero(); b[0] = 1/a[0]
    for n in range(1,N+1):
        s = sum(a[k]*b[n-k] for k in range(1,n+1))
        b[n] = -s/a[0]
    return b
def s_D(a): return [n*a[n] for n in range(N+1)]
def s_pow(a,k):
    r = s_one()
    for _ in range(k): r = s_mul(r,a)
    return r

# l-polynomials: list of series, index = power of l
def p_trim(p):
    while len(p)>1 and all(x==0 for x in p[-1]): p = p[:-1]
    return p
def p_add(p,q):
    n = max(len(p),len(q)); r=[]
    for i in range(n):
        a = p[i] if i<len(p) else s_zero(); b = q[i] if i<len(q) else s_zero()
        r.append(s_add(a,b))
    return p_trim(r)
def p_scal_series(s,p): return p_trim([s_mul(s,x) for x in p])
def p_D(p):
    # D(l^j f) = j l^{j-1} f + l^j D f
    r = [s_zero() for _ in p]
    for j,f in enumerate(p):
        r[j] = s_add(r[j], s_D(f))
        if j>0: r[j-1] = s_add(r[j-1], s_scal(Fr(j),f))
    return p_trim(r)
def p_mul(p,q):
    r = [s_zero() for _ in range(len(p)+len(q)-1)]
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            r[i+j] = s_add(r[i+j], s_mul(a,b))
    return p_trim(r)

# P(w,q) coefficients of w^m, m=0..4 (truncated in w at degree 4)
M = 5
def w_mul(a,b):
    c=[Fr(0)]*M
    for i in range(M):
        for j in range(M-i): c[i+j]+=a[i]*b[j]
    return c
def w_inv(a):
    b=[Fr(0)]*M; b[0]=1/a[0]
    for n in range(1,M):
        b[n] = -sum(a[k]*b[n-k] for k in range(1,n+1))/a[0]
    return b
pt = [s_zero() for _ in range(M)]
for d in range(N+1):
    num=[Fr(1)]+[Fr(0)]*(M-1)
    for k in range(1,5*d+1): num = w_mul(num,[Fr(k),Fr(5)]+[Fr(0)]*(M-2))
    den=[Fr(1)]+[Fr(0)]*(M-1)
    for k in range(1,d+1):
        # (w+k)^5 - w^5 = 5k w^4 + 10k^2 w^3 + 10 k^3 w^2 + 5 k^4 w + k^5
        den = w_mul(den,[Fr(k**5),Fr(5*k**4),Fr(10*k**3),Fr(10*k**2),Fr(5*k)])
    t = w_mul(num,w_inv(den))
    for m in range(M): pt[m][d] += t[m]
C0 = pt[0]
# check C0 and I1tilde against def:set-genus0
C0chk = [Fr(factorial(5*d), factorial(d)**5) for d in range(N+1)]
assert C0 == C0chk
I1t = [Fr(0)]+[5*Fr(factorial(5*d),factorial(d)**5)*sum(Fr(1,m) for m in range(d+1,5*d+1)) for d in range(1,N+1)]
assert pt[1] == I1t, (pt[1], I1t)
print("C0, I1tilde agree with p~_0, p~_1")

# p_k = sum_j l^j/j! p~_{k-j}
pk = []
for k in range(5):
    pk.append(p_trim([s_scal(Fr(1,factorial(j)), pt[k-j]) for j in range(k+1)]))
# PF4 check
q = s_zero(); q[1]=Fr(1)
def PF4(p):
    D1=p_D(p); D2=p_D(D1); D3=p_D(D2); D4=p_D(D3)
    u = s_sub(s_one(), s_scal(Fr(5**5),q))
    r = p_scal_series(u,D4)
    r = p_add(r, p_scal_series(s_scal(Fr(-6250),q),D3))
    r = p_add(r, p_scal_series(s_scal(Fr(-4375),q),D2))
    r = p_add(r, p_scal_series(s_scal(Fr(-1250),q),D1))
    r = p_add(r, p_scal_series(s_scal(Fr(-120),q),p))
    return r
for k in range(5):
    r = PF4(pk[k])
    print("PF4 p_%d ="%k, [ [str(x) for x in ser] for ser in r])

# Birkhoff
C0inv = s_inv(C0)
T = p_trim([s_mul(I1t,C0inv), s_one()])  # T = l + I1t/C0
C1p = p_D(T); assert len(C1p)==1; C1 = C1p[0]
C1inv = s_inv(C1)
def dT(p): return p_scal_series(C1inv, p_D(p))
P = [ [p_scal_series(C0inv, pk[k]) for k in range(5)] ]
kap = []
for a in range(4):
    row = P[a]
    kp = dT(row[a+1]); assert len(kp)==1, "kappa has l-dependence"
    kap.append(kp[0])
    kinv = s_inv(kp[0])
    new = [None]*5
    for k in range(a+1,5): new[k] = p_scal_series(kinv, dT(row[k]))
    P.append(new)
# L = (1-5^5 q)^{-1/5}
L = s_zero()
for n in range(N+1):
    c = Fr(1)
    for j in range(n): c *= (Fr(-1,5)-j)
    c /= factorial(n)
    L[n] = c*(-5**5)**n
L5 = s_pow(L,5)
k1 = s_mul(L5, s_inv(s_mul(s_pow(C0,2), s_pow(C1,3))))
k3 = s_mul(C0, C1inv)
print("kappa0 == 1:", kap[0]==s_one())
print("kappa1 == L^5/(C0^2 C1^3):", kap[1]==k1)
print("kappa2 == 1:", kap[2]==s_one())
print("kappa3 == C0/C1:", kap[3]==k3)
# D L = (L^6 - L)/5
print("DL == (L^6-L)/5:", s_D(L)==s_scal(Fr(1,5), s_sub(s_pow(L,6),L)))
print("kappa1 series:", [str(x) for x in kap[1]])
print("kappa3 series:", [str(x) for x in kap[3]])
