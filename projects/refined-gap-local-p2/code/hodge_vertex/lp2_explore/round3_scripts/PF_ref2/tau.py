from fractions import Fraction as F
def tau(v): # v: dict n->coef ; tau e_n = 3 e_{n-1} - n e_n
    d={}
    for n,c in v.items():
        d[n-1]=d.get(n-1,0)+3*c; d[n]=d.get(n,0)-n*c
    return d
ok=True
for n in range(-6,9):
    t=tau(tau(tau({n:1}))); t[n-3]=t.get(n-3,0)-27
    want={n-2:-27*(n-1),n-1:3*(3*n*n-3*n+1),n:-n**3}
    t={k:v for k,v in t.items() if v}; want={k:v for k,v in want.items() if v}
    ok&= t==want
print("tau^3-27sigma^3 formula:",ok)
# Check Omega(Omega-z)(Omega-2z) = z^3 (V-3)^3 sigma^3 numerically in V-rep with e_n=V^n/(a)_n, specific z,E
E=F(3,7); z=F(2,5); a=1+E/z
def poch(n):
    p=F(1)
    if n>=0:
        for j in range(n): p*=a+j
    else:
        for j in range(1,-n+1): p/= (a-j)
    return p
# represent functions as Laurent polys in V: dict pow->coef
def e(n): return {n:1/poch(n)}
def theta(f):
    d={}
    for n,c in f.items():
        if n: d[n-1]=d.get(n-1,0)+n*c; d[n]=d.get(n,0)-F(n,3)*c
    return d
def T(f):  # z theta + E/V
    d={k:z*v for k,v in theta(f).items()}
    for n,c in f.items(): d[n-1]=d.get(n-1,0)+E*c
    return {k:v for k,v in d.items() if v}
def lin(*pairs):
    d={}
    for c,f in pairs:
        for k,v in f.items(): d[k]=d.get(k,0)+c*v
    return {k:v for k,v in d.items() if v}
ok=True
for n in range(-4,6):
    lhs=T(e(n)); rhs=lin((z,e(n-1)),(-F(n)*z/3,e(n)))
    ok&= lhs==rhs
print("T e_n = z e_{n-1} - nz/3 e_n:",ok)
