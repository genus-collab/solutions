from fractions import Fraction as F
import sys
K=14
def recnums(n,pert=None):
    r=[F(1)]
    for m in range(n):
        r.append((3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m else 0))/(27*(m+1)))
    if pert: r[pert[0]]+=pert[1]
    return r
def add(a,b,c=1):
    d=dict(a)
    for k,v in b.items(): d[k]=d.get(k,0)+c*v
    return {k:v for k,v in d.items() if v!=0 and k[0]<=K}
def theta(f): # theta V^n = n V^{n-1} - n/3 V^n
    d={}
    for (zp,n),x in f.items():
        if n: 
            d[(zp,n-1)]=d.get((zp,n-1),0)+n*x
            d[(zp,n)]=d.get((zp,n),0)-F(n,3)*x
    return d
def zop(f,Ecoef,cst,zc): # (z theta + Ecoef/V + cst + zc*z) f
    d={}
    for (zp,n),x in theta(f).items(): d[(zp+1,n)]=d.get((zp+1,n),0)+x
    for (zp,n),x in f.items():
        d[(zp,n-1)]=d.get((zp,n-1),0)+Ecoef*x
        d[(zp,n)]=d.get((zp,n),0)+cst*x
        d[(zp+1,n)]=d.get((zp+1,n),0)+zc*x
    return {k:v for k,v in d.items() if v!=0 and k[0]<=K}
def series_inv_prod(l,E): # 1/prod_{j=1}^{l+1}(1+jz/E) as z-series up to K
    s=[F(0)]*(K+1); s[0]=F(1)
    for j in range(1,l+2):
        c=F(j)/E
        for k in range(1,K+1): s[k]-=c*s[k-1]
    return s
def build(E,r):
    f={}
    for l in range(K+1):
        ip=series_inv_prod(l,E)
        for k in range(K+1-l):
            if ip[k]: f[(l+k,l+1)]=f.get((l+k,l+1),0)+r[l]*ip[k]/E**(l+1)
    return f
def mulpoly(f,poly):
    d={}
    for (zp,n),x in f.items():
        for e,c in poly.items(): d[(zp,n+e)]=d.get((zp,n+e),0)+c*x
    return {k:v for k,v in d.items() if v!=0}
def check(E,r,alpha=None,ysign=1):
    if alpha is None: alpha=E/3
    f=build(E,r)
    # T = z theta + Y + alpha, Y = E/V - E/3
    T=lambda g: zop(g,E,alpha-E/3,0)
    A=T(T(T(f)))
    B=f
    for j in range(3): B=zop(B,E,-E/3,F(j,3))
    cube={3:F(1),2:F(-9),1:F(27),0:F(-27)}  # (V-3)^3
    R=add(mulpoly(A,cube),B,27*ysign)
    bad=sorted(set(k[0] for k in R if k[0]<=K-3))
    return bad
r=recnums(K+2)
for E in [F(1),F(5,2),F(-7,3)]:
    print("E",E,"nonzero z-orders:",check(E,r))
print("control pert r5:",check(F(1),recnums(K+2,(5,F(1,1000)))))
print("control alpha=E/3+1/5:",check(F(1),r,alpha=F(1,3)+F(1,5))[:5])
print("control ysign=-1:",check(F(1),r,ysign=-1)[:5])
