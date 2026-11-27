# Graded setting: frames from the definition pr_{K0}(z^k H^a), checked against Pi_i computed
# directly from the residue series I^C_m of def:gau-van (independent of eq:tra-values).
from fractions import Fraction as F
from math import factorial
nu=[0,-1,3,2,1]; eps=[1,-1,-1,1,1]
NMAX=40
def IC(m,n):  # [z^{-n}] I^C_m as (coef, vpower); from the explicit double sum of def:gau-van
    out={}
    for i in range(0,n+1):
        j=m+4*i+1
        if j<0: continue
        if i+j!=n: continue
        c=F((-1)**i,4**i*factorial(i))*F((-1)**j,factorial(j))
        out[j]=out.get(j,0)+c
    return out
def YC(b,n): return {e:eps[b]*c for e,c in IC(nu[b],n).items()}
def Pi0_mono(m,b):  # Pi_0(z^m H^b) = -[z^0] z^m YC_b(-z) = (-1)^{m+1}[z^{-m}]YC_b
    return {e:(-1)**(m+1)*c for e,c in YC(b,m).items()}
def padd(a,b,s=1):
    r=dict(a)
    for k,c in b.items():
        r[k]=r.get(k,0)+s*c
        if r[k]==0: del r[k]
    return r
def pmulmono(a,e0): return {e+e0:c for e,c in a.items()}
def Pi(i,q):  # q: dict (k,b)->poly in v
    s={}
    for (k,b),p in q.items():
        if k<i: continue
        for e,c in Pi0_mono(k-i,b).items():
            for e2,c2 in p.items(): s=padd(s,{e+e2:c*c2})
    return s
# closed form (eq:tra-values)
ok=True
for m in range(0,NMAX):
    for b in range(5):
        val=Pi0_mono(m,b)
        if (b+m)%5==1:
            r=m//5
            exp={m-r:F(-eps[b],4**r*factorial(r)*factorial(m-r))}
        else: exp={}
        if val!=exp: ok=False; print("values mismatch",m,b,val,exp)
print("eq:tra-values ok:",ok)
# frames by solving Pi_i(z^k H^a + sum_{j<k} f_j z^j H)=0, i in [0,k) via back-substitution (top-down)
def frame(k,a):
    f={}
    for i in range(k-1,-1,-1):
        q={(k,a):{0:F(1)}}
        for j,p in f.items(): q[(j,1)]=p
        r=Pi(i,q)          # current value with f_i unknown (f_i enters as + f_i * Pi_0(H) = f_i)
        f[i]={e:-c for e,c in r.items()}
    return {j:p for j,p in f.items() if p}
allok=True; trivial=[]
for a in [0,2,3,4]:
    for k in range(0,16):
        fr=frame(k,a)
        q={(k,a):{0:F(1)}}
        for j,p in fr.items(): q[(j,1)]=p
        for i in range(0,k+6):
            if Pi(i,q): allok=False; print("not in K0",k,a,i)
        for j,p in fr.items():
            for e in p:
                if 5*e!=4*(k-j)+nu[a]+1 or e<1: allok=False; print("exponent",k,a,j,e)
        if not fr: trivial.append((k,a))
print("frames in K0 and exponents ok:",allok)
print("trivial frames:",sorted(trivial))
print("fr(6,0):",frame(6,0)," fr(4,2):",frame(4,2)," fr(2,4):",frame(2,4)," fr(3,3):",frame(3,3)," fr(1,0):",frame(1,0))
