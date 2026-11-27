# Independent sanity check of the frame formulas used by lem:tra-conv (a),(b2),(e2).
# Polynomials in v with rational coefficients are dicts {exponent: Fraction}.
from fractions import Fraction as F
from math import factorial
nu  = [0,-1,3,2,1]
eps = [1,-1,-1,1,1]

def YC_coef(n,b):
    """[z^{-n}] Y^C_b as (coef, v-exponent) from the residue definition (def:gau-van), not the closed form."""
    # I^C_m = sum_i (-1/(4z))^i/i! * (-v/z)^(m+4i+1)/(m+4i+1)!, m = nu_b ; power z^{-(i + m+4i+1)}
    m = nu[b]
    tot = {}
    for i in range(0, n+1):
        j = m+4*i+1
        if j < 0: continue
        if i + j == n:
            c = F((-1)**i, 4**i*factorial(i)) * F((-1)**j, factorial(j))
            tot[j] = tot.get(j,0) + eps[b]*c
    return tot

def Pi0(m,b):
    """Pi_0(z^m H^b) = (-1)^{m+1} [z^{-m}] Y^C_b   (eq:tra-shift)."""
    return {e: (-1)**(m+1)*c for e,c in YC_coef(m,b).items()}

def padd(p,q,s=1):
    r = dict(p)
    for e,c in q.items():
        r[e] = r.get(e,0)+s*c
        if r[e]==0: del r[e]
    return r
def pmul(p,q):
    r={}
    for e1,c1 in p.items():
        for e2,c2 in q.items():
            r[e1+e2]=r.get(e1+e2,0)+c1*c2
    return {e:c for e,c in r.items() if c!=0}

def Pi(i, vec):
    """vec: dict {(k,b): poly}; Pi_i(vec) = sum_{k>=i} vec[(k,b)] * Pi_0(z^{k-i} H^b)."""
    tot={}
    for (k,b),p in vec.items():
        if k>=i:
            tot = padd(tot, pmul(p, Pi0(k-i,b)))
    return tot

def frame(k,a):
    """solve Pi_i(z^k H^a + sum_{j<k} f_j z^j H)=0, i=k-1..0, top-down."""
    vec = {(k,a): {0:F(1)}}
    for i in range(k-1,-1,-1):
        r = Pi(i, vec)          # f_i enters with Pi_0(H)=1
        if r:
            vec[(i,1)] = {e:-c for e,c in r.items()}
    return vec

# values quoted in the text
assert Pi0(0,1)=={0:F(1)}
assert Pi0(1,0)=={1:F(-1)}
assert Pi0(6,0)=={5:F(-1,480)}
assert Pi0(5,1)=={4:F(1,96)}
# closed form eq:tra-values
for b in range(5):
    for m in range(0,40):
        p = Pi0(m,b)
        if (b+m)%5==1:
            r=m//5
            assert p=={m-r: F(-eps[b], 4**r*factorial(r)*factorial(m-r))}, (m,b,p)
        else:
            assert p=={}, (m,b,p)
# examples of lem:tra-frames
assert frame(1,0)=={(1,0):{0:F(1)}, (0,1):{1:F(1)}}
assert frame(2,4)=={(2,4):{0:F(1)}, (0,1):{2:F(1,2)}}
assert frame(3,3)=={(3,3):{0:F(1)}, (0,1):{3:F(1,6)}}
assert frame(4,2)=={(4,2):{0:F(1)}, (0,1):{4:F(-1,24)}}
assert frame(6,0)=={(6,0):{0:F(1)}, (5,1):{1:F(1)}, (0,1):{5:F(-1,120)}}
triv=[]
for a in [0,2,3,4]:
    for k in range(0,26):
        fr = frame(k,a)
        # all Pi_i vanish
        for i in range(0,k+3):
            assert Pi(i,fr)=={}, (k,a,i)
        # exponent of f_j is (4(k-j)+nu_a+1)/5 >= 1, single monomial
        for (j,b),p in fr.items():
            if (j,b)==(k,a): continue
            assert b==1 and j<k and len(p)==1
            e=list(p)[0]
            assert 5*e==4*(k-j)+nu[a]+1 and e>=1
        if len(fr)==1: triv.append((k,a))
print("trivial frames:",sorted(triv, key=lambda t:(t[1],t[0])))
# weight bookkeeping: wt(k,a)=(4k+nu_a-4)/5 ; frame component f_j z^j H has v-exp + wt(j,1) = wt(k,a)
wt=lambda k,b: F(4*k+nu[b]-4,5)
for a in [0,2,3,4]:
    for k in range(0,26):
        for (j,b),p in frame(k,a).items():
            e=list(p)[0]
            assert e+wt(j,b)==wt(k,a)
assert min(wt(k,b) for k in range(2,10) for b in range(5))==F(3,5)
assert min(wt(k,b) for k in range(2,10) for b in [0,2,3,4])==F(4,5)
assert -1 not in H if False else True
# -z1 and H in K ; Pi_0(-z1)=v
m1={(1,0):{0:F(-1)}}
assert Pi(0,m1)=={1:F(1)} and all(Pi(i,m1)=={} for i in range(1,6))
Hv={(0,1):{0:F(1)}}
assert Pi(0,Hv)=={0:F(1)} and all(Pi(i,Hv)=={} for i in range(1,6))
print("all frame checks passed")
