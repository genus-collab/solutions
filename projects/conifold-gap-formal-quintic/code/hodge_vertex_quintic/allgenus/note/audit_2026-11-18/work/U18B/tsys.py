from fractions import Fraction as F
N=6
def T(k,f):  # f list of coeffs in t
    out=[F(0)]*(N+1)
    for n,c in enumerate(f):
        out[n]-=c
        if n+1<=N: out[n+1]+=F(k-4*n,5)*c
    return out
def neg(f): return [-c for c in f]
# unknown F0 = 1 + sum a_n t^n ; solve order by order via linearity
F0=[F(0)]*(N+1); F0[0]=F(1)
def closing(F0):
    F1=T(-1,F0); F2=neg(T(-2,F1)); F3=T(2,F2); F4=neg(T(1,F3)); G=neg(T(0,F4))
    return [g-f for g,f in zip(G,F0)],(F1,F2,F3,F4)
for n in range(1,N+1):
    # residual at t^{n+1} depends on a_n linearly with pivot 4n (zeta^0)
    F0[n]=F(0); r0,_=closing(F0)
    F0[n]=F(1); r1,_=closing(F0)
    piv=r1[n+1-1+1] - r0[n+1-1+1] if n+1<=N else None
    if n+1>N: break
    a=-r0[n+1]/(r1[n+1]-r0[n+1]); F0[n]=a
res,(F1,F2,F3,F4)=closing(F0)
print("F0:",F0[:N])
print("residual:",res[:N])
print("F_b(0):",F0[0],F1[0],F2[0],F3[0],F4[0])
print(float(F0[3]), float(F0[4]))
