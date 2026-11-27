from fractions import Fraction as Fr
N=6
def T(k,F):  # T_k F = -F + t/5 (k F - 4 t F')
    out=[Fr(0)]*(N+1)
    for n,c in enumerate(F):
        out[n]-=c
        if n+1<=N: out[n+1]+=Fr(k-4*n,5)*c
    return out
def neg(F): return [-c for c in F]
# unknown F0 coefficients solved order by order via closing relation
F0=[Fr(1)]+[Fr(0)]*N
for n in range(1,N):
    # choose F0[n] so that closing relation holds at order n+1
    def resid(F0):
        F1=T(-1,F0); F2=neg(T(-2,F1)); F3=T(2,F2); F4=neg(T(1,F3)); G=neg(T(0,F4))
        return [G[i]-F0[i] for i in range(N+1)]
    r0=resid(F0)[n+1]
    F0[n]+=1; r1=resid(F0)[n+1]; F0[n]-=1
    F0[n]-= r0/(r1-r0)
F1=T(-1,F0); F2=neg(T(-2,F1)); F3=T(2,F2); F4=neg(T(1,F3))
print("F0",F0[:4]); print("F1",F1[:3]);print("F2",F2[:3]);print("F3",F3[:3]);print("F4",F4[:3])
