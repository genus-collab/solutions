# Independent check: L_pf conjugated, applied to v*S, as Laurent series in V (z-degree <= K, V-degree <= NV)
from fractions import Fraction as F
import sys
K=8; NV=14
def rec(L):
    r=[F(1)]; rm1=F(0)
    for m in range(L):
        prev = r[m-1] if m>=1 else F(0)
        r.append((3*(3*m*m+3*m+1)*r[m]-m**3*prev)/(27*(m+1)) + (F(1,1000) if m==3 else 0))
    return r
def run(E):
    r=rec(NV+2)
    # 1/(a)_n = (z/E)^n / prod_{j=1}^n (1+j z/E), series in z
    f={}
    for n in range(1,NV+1):
        ser=[F(0)]*(K+2); ser[0]=F(1)
        for j in range(1,n+1):
            c=F(-j)/E  # multiply by 1/(1+jz/E) = sum (-j z/E)^k
            new=[F(0)]*(K+2)
            for a in range(K+2):
                if ser[a]==0: continue
                p=F(1)
                for b in range(0,K+2-a):
                    new[a+b]+=ser[a]*p; p*=c
            ser=new
        for a in range(K+2):
            kz=a+n-1  # z^{-1} * z^n * z^a
            if kz<=K and ser[a]!=0:
                f[(kz,n)]=f.get((kz,n),0)+r[n-1]*ser[a]/E**n
    def add(A,B,s=1):
        C=dict(A)
        for k,v in B.items(): C[k]=C.get(k,0)+s*v
        return C
    def scal(A,c): return {k:v*c for k,v in A.items()}
    def zmul(A): return {(k+1,n):v for (k,n),v in A.items() if k+1<=K}
    def theta(A):  # (1-V/3) d/dV
        C={}
        for (k,n),v in A.items():
            if n==0: continue
            C[(k,n-1)]=C.get((k,n-1),0)+n*v
            C[(k,n)]=C.get((k,n),0)-F(n,3)*v
        return C
    def T(A):
        C=zmul(theta(A))
        for (k,n),v in A.items(): C[(k,n-1)]=C.get((k,n-1),0)+E*v
        return C
    alpha=E/3
    def Tsh(A,c): # T - E/3 + alpha... (z theta + Y + jz/3) = T - alpha + j z/3 ; here c=j
        C=add(T(A),scal(A,-alpha))
        return add(C,scal(zmul(A),F(c,3)))
    # y = (V-3)^{-3} = -1/27 (1-V/3)^{-3}
    yser=[]
    for m in range(NV+12):
        # coefficient of V^m in (1-x)^{-3}, x=V/3: C(m+2,2)/3^m
        yser.append(F(-1,27)*F((m+2)*(m+1)//2,3**m))
    def ymul(A):
        C={}
        for (k,n),v in A.items():
            for m in range(0,min(NV+4-n,len(yser))):
                C[(k,n+m)]=C.get((k,n+m),0)+v*yser[m]
        return C
    t1=T(T(T(f)))
    g=Tsh(Tsh(Tsh(f,2),1),0)
    t2=scal(ymul(g),27)
    res=add(t1,t2)
    bad=[(k,n,v) for (k,n),v in res.items() if v!=0 and n<=NV-4 and k<=K]
    return bad, len(res)
for E in [F(1),F(5,2),F(-7,3)]:
    bad,l=run(E); print(E, "nonzero:",bad[:5], "terms",l)
