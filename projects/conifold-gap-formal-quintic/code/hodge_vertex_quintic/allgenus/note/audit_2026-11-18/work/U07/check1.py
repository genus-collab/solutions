# Independent sanity checks for unit U07 (transport section). Exact arithmetic.
from fractions import Fraction as Fr
from math import factorial
import random
random.seed(7)

nu  = [0,-1,3,2,1]
eps = [1,-1,-1,1,1]
NV = 9   # work modulo v^NV

# ---------- polynomials in v mod v^NV : lists of Fractions
def pz(): return [Fr(0)]*NV
def padd(a,b): return [x+y for x,y in zip(a,b)]
def pscal(c,a): return [c*x for x in a]
def pmul(a,b):
    r=[Fr(0)]*NV
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j<NV and y!=0: r[i+j]+=x*y
    return r
def pmono(c,e):
    r=pz()
    if e<NV: r=list(r); r[e]=Fr(c)
    return r
def piszero(a): return all(x==0 for x in a)

# ---------- vanishing solution from the residue definition (def:gau-van)
# I_m = Res_{p=0} p^m exp(-p^4/(4z)) exp(-v/(pz)) dp : sum_{i,j} with m+4i-j=-1
def IC(m):
    # returns dict n -> poly in v  (coefficient of z^{-n})
    out={}
    for i in range(0,40):
        j=m+4*i+1
        if j<0: continue
        if j>=NV: continue
        c=Fr((-1)**i,4**i*factorial(i))*Fr((-1)**j,factorial(j))
        n=i+j
        out[n]=padd(out.get(n,pz()),pmono(c,j))
    return out
YC=[ {n:pscal(eps[b],c) for n,c in IC(nu[b]).items()} for b in range(5)]

# check closed form of lem:gau-vanshape
for b in range(5):
    for n in range(0,12):
        r=n//5
        if n%5==(nu[b]+1)%5 and n-r<NV:
            exp=pmono(Fr(eps[b]*(-1)**n,4**r*factorial(r)*factorial(n-r)),n-r)
        else: exp=pz()
        got=YC[b].get(n,pz())
        if n-r<NV or n%5!=(nu[b]+1)%5:
            assert got==exp,(b,n)
print("closed form of YC ok")

# ---------- Pi_i from the definition: -[z^i] <YC(-z), q(z)>, q = dict (k,b)->poly in v
def Pi(i,q):
    tot=pz()
    for (k,b),c in q.items():
        n=k-i   # need coefficient of z^{-(k-i)} of YC_b(-z)
        if n<0: continue
        y=YC[b].get(n)
        if y is None: continue
        tot=padd(tot,pmul(c,pscal(-(-1)**n,y)))   # YC_b(-z): z^{-n} -> (-1)^n
    return tot

def Pi0mono(m,b): return Pi(0,{(m,b):pmono(1,0)})
# eq:tra-values and examples
assert Pi0mono(0,1)==pmono(1,0)
assert Pi0mono(1,0)==pmono(-1,1)
assert Pi0mono(6,0)==pmono(Fr(-1,480),5)
for m in range(0,10):
    for b in range(5):
        r=m//5
        exp = pmono(Fr(-eps[b],4**r*factorial(r)*factorial(m-r)),m-r) if (b+m)%5==1 else pz()
        assert Pi0mono(m,b)==exp,(m,b)
print("eq:tra-values ok")

# ---------- frames via Theta^{-1}: solve Pi_i(z^k H^a + sum f_j z^j H)=0 by downward recursion on i (independent: uses Pi only)
def frame(k,a):
    q={(k,a):pmono(1,0)}
    for i in range(k-1,-1,-1):
        # Pi_i(q + f_i z^i H) = Pi_i(q) + f_i  (Pi_0(H)=1), higher f already in q
        val=Pi(i,q)
        q[(i,1)]=pscal(-1,val)
    return q
def show(q):
    s=[]
    for (k,b),c in sorted(q.items()):
        for e,x in enumerate(c):
            if x!=0: s.append("%s v^%d z^%d H^%d"%(x,e,k,b))
    return " + ".join(s)
for (k,a) in [(1,0),(2,4),(3,3),(4,2),(6,0)]:
    print("frame",(k,a),":",show(frame(k,a)))
triv=[]
for a in [0,2,3,4]:
    for k in range(0,8):
        f=frame(k,a)
        # check all Pi_i vanish (i up to k+1) modulo v^NV
        for i in range(0,k+2): assert piszero(Pi(i,f)),(k,a,i)
        if all(piszero(c) for (kk,bb),c in f.items() if (kk,bb)!=(k,a)): triv.append((k,a))
        # exponent check
        for (j,bb),c in f.items():
            if (j,bb)==(k,a): continue
            for e,x in enumerate(c):
                if x!=0: assert 5*e==4*(k-j)+nu[a]+1 and e>=1,(k,a,j,e)
print("trivial frames:",sorted(triv))
assert sorted(triv)==sorted([(0,0),(0,2),(1,2),(2,2),(3,2),(0,3),(1,3),(2,3),(0,4),(1,4)])
print("frames ok")
