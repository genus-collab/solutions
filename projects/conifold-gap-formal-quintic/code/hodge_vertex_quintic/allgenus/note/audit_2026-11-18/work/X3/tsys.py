from fractions import Fraction as Fr
N=8
# colour-0 T-system: T_k = -1 + (t/5)(k - 4 t d/dt); series as lists
def T(k,F):
    G=[Fr(0)]*(N+1)
    for n,c in enumerate(F):
        if c==0: continue
        G[n]+= -c
        if n+1<=N: G[n+1]+= Fr(k-4*n,5)*c
    return G
def neg(F): return [-c for c in F]
# unknown F0 coefficients solved order by order
a=[Fr(1)]+[Fr(0)]*N
for n in range(1,N+1):
    # closing at order t^{n}: pivot 4(n-1) times a_{n-1}... determine a_{n-1}? do linear solve on a_n by trying
    pass
# simpler: solve via unknown symbolic linear: closing relation coefficient of t^{m+1} = 4m a_m + known
def closing(a):
    F0=a
    F1=T(-1,F0); F2=neg(T(-2,F1)); F3=T(2,F2); F4=neg(T(1,F3)); G=neg(T(0,F4))
    return [G[i]-F0[i] for i in range(N+1)], (F0,F1,F2,F3,F4)
for m in range(1,N):
    a[m]=Fr(0)
    r0,_=closing(a)
    a[m]=Fr(1)
    r1,_=closing(a)
    # residual at t^{m+1} is linear in a_m
    c0=r0[m+1]; c1=r1[m+1]-r0[m+1]
    a[m]=-c0/c1
r,F=closing(a)
print("residual",r[:N])
for b in range(5): print("Fs_%d"%b,[str(x) for x in F[b][:5]])
