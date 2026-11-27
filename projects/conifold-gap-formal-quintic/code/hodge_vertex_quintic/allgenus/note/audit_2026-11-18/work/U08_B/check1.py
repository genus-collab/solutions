# Independent sanity checks for unit U08 (lem:tra-conv, Theorem R), exact arithmetic.
from fractions import Fraction as F
from math import factorial, comb

nu  = [0,-1,3,2,1]
eps = [1,-1,-1,1,1]

# --- vanishing solution: coefficient of z^{-n} in YC_b is  c(n,b) * v^{n - floor(n/5)}
def yc(n,b):
    """returns (coefficient, v-exponent) of [z^-n] YC_b  (closed form of lem:gau-vanshape)"""
    if n % 5 != (nu[b]+1) % 5: return (F(0),0)
    m = n//5
    return (F(eps[b]*(-1)**n, 4**m*factorial(m)*factorial(n-m)), n-m)

# independent derivation from the residue definition: IC_m = sum_i (-1/(4z))^i/i! (-v/z)^{m+4i+1}/(m+4i+1)!
def yc_res(n,b):
    mm = nu[b]
    tot = F(0); ex=None
    for i in range(0, n+2):
        j = mm+4*i+1
        if j < 0: continue
        if i + j == n:
            tot += F((-1)**i,4**i*factorial(i)) * F((-1)**j, factorial(j)) * eps[b]
            ex = j
    return (tot, ex if ex is not None else 0)

for n in range(0,40):
    for b in range(5):
        a = yc(n,b); c = yc_res(n,b)
        assert a[0]==c[0] and (a[0]==0 or a[1]==c[1]), (n,b,a,c)
        if a[0]!=0:
            assert a[1] == n - n//5 and 5*a[1] >= 4*n
print("closed form of [z^-n]YC_b = residue definition, exponent n-floor(n/5) >= 4n/5 : OK (n<40)")

# model system z d/dv YC_b = bm_b YC_{b+1},  bm = (1,-1/v,1,-1,-1)
for n in range(0,40):
    for b in range(5):
        c,e = yc(n,b)            # term c v^e z^-n ; z d/dv -> c e v^{e-1} z^{-(n-1)}
        lhs = (c*e, e-1) if c!=0 and e>0 else (F(0),0)
        if n>=1:
            c2,e2 = yc(n-1,(b+1)%5)
            if b==1: rhs = (-c2, e2-1)
            elif b in (0,2): rhs=(c2,e2)
            else: rhs=(-c2,e2)
            if c2==0: rhs=(F(0),0)
            assert lhs[0]==rhs[0] and (lhs[0]==0 or lhs[1]==rhs[1]), (n,b,lhs,rhs)
print("YC solves the model system: OK")

# --- Pi_0(z^m H^b) = (-1)^{m+1} [z^-m] YC_b   -> polynomial c * v^e
def Pi0(m,b):
    c,e = yc(m,b)
    return ((-1)**(m+1)*c, e)
assert Pi0(0,1)==(F(1),0)
assert Pi0(1,0)==(F(-1),1)
assert Pi0(6,0)==(F(-1,480),5)
# value (i): Pi_0(z^{5r} H) = pi_r v^{4r}
for r in range(0,6):
    assert Pi0(5*r,1)==(F(1,4**r*factorial(r)*factorial(4*r)),4*r)

# --- frames by the recursion f_i = -Pi_0(z^{k-i}H^a) - sum_{r>=1} pi_r v^{4r} f_{i+5r}
def pi_r(r): return F(1,4**r*factorial(r)*factorial(4*r))
def frame(k,a):
    f = {}   # j -> (coef, v-exponent)
    for i in range(k-1,-1,-1):
        c,e = Pi0(k-i,a)
        tot = {}
        if c!=0: tot[e] = tot.get(e,0) - c
        r=1
        while i+5*r < k:
            cc,ee = f[i+5*r]
            if cc!=0: tot[ee+4*r] = tot.get(ee+4*r,0) - pi_r(r)*cc
            r+=1
        tot = {e_:c_ for e_,c_ in tot.items() if c_!=0}
        assert len(tot)<=1
        f[i] = (list(tot.values())[0], list(tot.keys())[0]) if tot else (F(0),0)
    return f

def Pi_i_of_frame(i,k,a,f):
    # Pi_i(z^k H^a + sum_j f_j z^j H) = Pi_0(z^{k-i}H^a) + sum_{j>=i} f_j Pi_0(z^{j-i} H)
    tot={}
    if k>=i:
        c,e = Pi0(k-i,a)
        if c!=0: tot[e]=tot.get(e,0)+c
    for j,(cc,ee) in f.items():
        if j>=i and cc!=0:
            c,e = Pi0(j-i,1)
            if c!=0: tot[e+ee]=tot.get(e+ee,0)+c*cc
    return {e_:c_ for e_,c_ in tot.items() if c_!=0}

trivial=[]
for k in range(0,26):
    for a in (0,2,3,4):
        f = frame(k,a)
        for i in range(0,k+2):
            assert Pi_i_of_frame(i,k,a,f)=={}, (k,a,i)
        nz = {j:ce for j,ce in f.items() if ce[0]!=0}
        for j,(c,e) in nz.items():
            assert 5*e == 4*(k-j)+nu[a]+1 and e>=1, (k,a,j,e)
        if not nz: trivial.append((k,a))
print("frames lie in K_0 (all Pi_i vanish), exponents (4(k-j)+nu_a+1)/5 >= 1 : OK (k<26)")
print("trivial frames:", sorted(trivial))
print("frame(1,0):",frame(1,0)); print("frame(2,4):",frame(2,4)); print("frame(3,3):",frame(3,3)); print("frame(4,2):",frame(4,2)); print("frame(6,0):",frame(6,0))
# weights: wt(k,a) for a!=1,k>=2 >= 4/5 ; wt(k,b) k>=2 >= 3/5
assert min(F(4*k+nu[a]-4,5) for k in range(2,10) for a in (0,2,3,4))==F(4,5)
assert min(F(4*k+nu[a]-4,5) for k in range(2,10) for a in range(5))==F(3,5)
print("weights OK")
