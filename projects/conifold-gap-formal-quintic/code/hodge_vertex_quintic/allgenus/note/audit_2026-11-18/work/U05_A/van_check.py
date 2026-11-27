# Sanity checks: (EP) for v = u tau^(-5/4); Y^C_{tilde A} (iterated integrals) == Y^C(v(u),z) (closed formula);
# leading-coefficient recursion of Lemma gau-vanq; (Mod1),(Mod2) for the residues I^C_m.
from fractions import Fraction as Fr
from math import factorial
N=14  # truncation order in u (keep coefficients u^0..u^N)
def mul(a,b):
    r=[Fr(0)]*(N+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>N: break
            r[i+j]+=x*y
    return r
def add(a,b,c=1): return [x+c*y for x,y in zip(a,b)]
def inv(a):
    r=[Fr(0)]*(N+1); r[0]=1/a[0]
    for n in range(1,N+1):
        r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
    return r
def powfrac(a,e):  # a in 1+u Q[[u]], exponent Fraction e ; via (a^e)' a = e a' a^e
    assert a[0]==1
    r=[Fr(0)]*(N+1); r[0]=Fr(1)
    for n in range(1,N+1):
        r[n]=sum((e*k-(n-k))*a[k]*r[n-k] for k in range(1,n+1))/n
    return r
def der(a): return [a[k+1]*(k+1) for k in range(N)]+[Fr(0)]
def integ(a): return [Fr(0)]+[a[k]/(k+1) for k in range(N)]
def shift_down(a):  # a/u, requires a[0]==0
    assert a[0]==0
    return a[1:]+[Fr(0)]
def shift_up(a): return [Fr(0)]+a[:N]
one=[Fr(1)]+[Fr(0)]*N
u=[Fr(0),Fr(1)]+[Fr(0)]*(N-1)
tau=inv([Fr(4,5*n+4) for n in range(N+1)])
w=powfrac(tau,Fr(-5,4))          # v = u*w
v=shift_up(w)
vp=der(v)
# (EP): v'^5 * u * (1-u)^5 == v   (valid through u^(N-1) because der loses one order)
lhs=one
for _ in range(5): lhs=mul(lhs,vp)
omu=add(one,u,-1)
for _ in range(5): lhs=mul(lhs,omu)
lhs=shift_up(lhs)
print("(EP) v'^5 u (1-u)^5 == v through u^%d:"%(N-1), lhs[:N]==v[:N], "  tau =",[str(x) for x in tau[:3]])
# tilde b: (v', -v'/v, v', -v', -v');  u*tilde b_1 = -v'/w
ub1=[-x for x in mul(vp,inv(w))]
tb={0:vp,2:vp,3:[-x for x in vp],4:[-x for x in vp]}
# iterated integrals y^(n)
nu=[0,-1,3,2,1]; eps=[1,-1,-1,1,1]
y=[[Fr(0)]*(N+1) for _ in range(5)]; y[1]=[Fr(-1)]+[Fr(0)]*N
ys=[y]
NZ=12
for n in range(NZ):
    new=[]
    for a in range(5):
        src=ys[-1][(a+1)%5]
        if a==1:
            integrand=mul(ub1,shift_down(src))   # asserts no pole/log
        else:
            integrand=mul(tb[a],src)
        new.append(integ(integrand))
    ys.append(new)
# closed formula
def vpow(k):
    r=one
    for _ in range(k): r=mul(r,v)
    return r
ok=True; okord=True
M=N-2  # safe comparison range (derivative/shift lose orders)
for n in range(NZ+1):
    for a in range(5):
        if (n-(nu[a]+1))%5==0:
            m=n//5
            c=Fr(eps[a]*(-1)**n, 4**m*factorial(m)*factorial(n-m))
            closed=[c*x for x in vpow(n-m)]
            e=(4*n+nu[a]+1)//5
            assert (4*n+nu[a]+1)%5==0 and e==n-m
            # exact order
            lead=[k for k,x in enumerate(ys[n][a]) if x!=0]
            if lead and lead[0]!=e: okord=False
        else:
            closed=[Fr(0)]*(N+1)
        if ys[n][a][:M]!=closed[:M]: ok=False
print("Y^C_{tilde A} == Y^C(v(u),z) for z^-n, n<=%d, through u^%d:"%(NZ,M-1), ok, "; exact orders (4n+nu_a+1)/5:", okord)
# leading-coefficient recursion with b^o=(1,-1,1,-1,-1)
bo=[1,-1,1,-1,-1]
c={(0,1):Fr(-1)}
okc=True
for n in range(0,25):
    for a in range(5):
        if (n+1-(1-a))%5==0:
            e=Fr(4*(n+1)+nu[a]+1,5)
            c[(n+1,a)]=bo[a]*c[(n,(a+1)%5)]/e
            m=(n+1)//5
            if c[(n+1,a)]!=Fr(eps[a]*(-1)**(n+1),4**m*factorial(m)*factorial(n+1-m)): okc=False
print("leading coefficients c^(n)_a == eps_a(-1)^n/(4^m m!(n-m)!) for n<=25:", okc)
# (Mod1),(Mod2) on the explicit residues, as polynomials in v, zi=1/z up to total v-degree VD
VD=30
def IC(m):
    r={}
    i=0
    while True:
        j=m+4*i+1
        if j>VD: break
        if j>=0:
            r[(j,i+j)]=Fr((-1)**i*(-1)**j, 4**i*factorial(i)*factorial(j))  # v^j zi^(i+j)
        i+=1
    return r
def dv_times_z(p): return {(j-1,k-1):x*j for (j,k),x in p.items() if j>0}
def padd(a,b,c=1):
    r=dict(a)
    for kk,x in b.items():
        r[kk]=r.get(kk,0)+c*x
        if r[kk]==0: del r[kk]
    return r
def trunc(p,d): return {kk:x for kk,x in p.items() if kk[0]<=d}
ok1=all(trunc(dv_times_z(IC(m)),VD-6)==trunc({kk:-x for kk,x in IC(m-1).items()},VD-6) for m in range(-6,8))
def zmul(p): return {(j,k-1):x for (j,k),x in p.items()}
def vmul(p): return {(j+1,k):x for (j,k),x in p.items()}
ok2=True
for m in range(-6,8):
    t=padd(padd({kk:m*x for kk,x in zmul(IC(m-1)).items()}, IC(m+3), -1), vmul(IC(m-2)))
    if trunc(t,VD-6)!={}: ok2=False
print("(Mod1) for m=-6..7:", ok1, "  (Mod2) for m=-6..7:", ok2)
for m in range(-6,8):
    t=padd(padd({kk:m*x for kk,x in zmul(IC(m-1)).items()}, IC(m+3), -1), vmul(IC(m-2)))
    t=trunc(t,VD-6)
    if t!={}: print("m=",m, sorted(t.items())[:4])
ok2=True
for m in range(-9,12):
    t=padd(padd({kk:m*x for kk,x in zmul(IC(m-1)).items() if m*x!=0}, IC(m+3), -1), vmul(IC(m-2)))
    if trunc(t,VD-6)!={}: ok2=False
print("(Mod2) recheck (zero entries filtered) for m=-9..11:", ok2)
