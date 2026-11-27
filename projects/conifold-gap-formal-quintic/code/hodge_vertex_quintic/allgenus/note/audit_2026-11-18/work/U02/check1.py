from fractions import Fraction as Fr
import cmath
N=12
# T-system colour 0: T_k F = -F + (t/5)(k F - 4 t F')
def T(k,F):
    out=[Fr(0)]*(N+1)
    for n,c in enumerate(F):
        if c==0: continue
        out[n]-=c
        if n+1<=N: out[n+1]+=Fr(k-4*n,5)*c
    return out
def neg(F): return [-c for c in F]
# solve closing relation F0 = -T0 T1 T2 T-2 T-1 F0 by recursion on coefficients
F0=[Fr(0)]*(N+1); F0[0]=Fr(1)
def chain(F0):
    F1=T(-1,F0); F2=neg(T(-2,F1)); F3=T(2,F2); F4=neg(T(1,F3)); G=neg(T(0,F4))
    return F1,F2,F3,F4,G
for n in range(1,N):
    # coefficient of t^{n+1} of (O-1)F0 is 4n*c_n + lower; set c_n so that it vanishes
    F0[n]=Fr(0)
    *_,G=chain(F0)
    r=G[n+1]-F0[n+1]
    F0[n]=-r/(4*n)
F1,F2,F3,F4,G=chain(F0)
print("Fs0:",F0[:5])
print("closing residual (should be 0 up to t^N-1):",[G[k]-F0[k] for k in range(N)])
Fs=[F0,F1,F2,F3,F4]
print("Fs_b(0):",[F[0] for F in Fs])
# symplecticity N_ij(t) numerically
z=cmath.exp(2j*cmath.pi/5)
def Rhat(i,b,coefsign=1):
    # coefficients of t^k of zeta^{i(b+1)} Fs_b(zeta^{-i} t)
    return [z**(i*(b+1))*complex(Fs[b][k])*z**(-i*k) for k in range(N)]
maxerr=0
for i in range(5):
  for j in range(5):
    for n in range(N-1):
      s=0
      for a in range(5):
        A=Rhat(i,a); B=Rhat(j,(3-a)%5)
        for k in range(n+1):
          s+=A[k]*(-1)**k*B[n-k]
      s/=5
      target=1 if (i==j and n==0) else 0
      maxerr=max(maxerr,abs(s-target))
print("max |N_ij - delta| over orders <",N-1,":",maxerr)
