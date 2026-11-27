from fractions import Fraction as Fr
N=10
# colour-0 T-system. T_k F = -F + (t/5)(k F - 4 t F')
def T(k,F):
    G=[Fr(0)]*(N+1)
    for n,c in enumerate(F):
        G[n]+= -c
        if n+1<=N: G[n+1]+= Fr(k-4*n,5)*c
    return G
def chain(F0):
    F1=T(-1,F0); F2=[-x for x in T(-2,F1)]; F3=T(2,F2); F4=[-x for x in T(1,F3)]; C=[-x for x in T(0,F4)]
    return F1,F2,F3,F4,C
# solve closing relation order by order
F0=[Fr(1)]+[Fr(0)]*N
for n in range(1,N):
    # choose F0[n] so that coefficient of t^{n+1} in C-F0 vanishes
    F0[n]=Fr(0)
    C=chain(F0)[4]
    d0=C[n+1]-F0[n+1]
    F0[n]=Fr(1)
    C=chain(F0)[4]
    d1=C[n+1]-F0[n+1]
    piv=d1-d0
    assert piv==4*n, (n,piv)
    F0[n]=-d0/piv
F1,F2,F3,F4,C=chain(F0)
print("F0",F0[:5]); print("closing defect",[C[i]-F0[i] for i in range(N)])
Fs=[F0,F1,F2,F3,F4]
print("F_b(0)",[F[0] for F in Fs])
# cyclotomic arithmetic: elements as length-5 lists of Fractions (coeffs of zeta^k), reduced via equality test modulo 1+z+..+z^4
def cz(k,c=Fr(1)):
    v=[Fr(0)]*5; v[k%5]=c; return v
def cadd(a,b): return [x+y for x,y in zip(a,b)]
def cmul(a,b):
    v=[Fr(0)]*5
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            v[(i+j)%5]+=x*y
    return v
def cnorm(a):
    m=a[4]; return [x-m for x in a[:4]]
def ceq(a,b): return cnorm(a)==cnorm(b)
# Rhat^{(i)}_b(t) = zeta^{i(b+1)} F_b(zeta^{-i} t); coefficient of t^n: zeta^{i(b+1)-i n} F_b[n]
M=N-1
def Rh(b,i,n,sign=1):
    return cz(i*(b+1)-i*n, Fs[b][n]*(sign**n))
# N_ij(t) = sum_{a} (1/5) Rhat_a^{(i)}(-t) Rhat_{3-a}^{(j)}(t)
ok=True
for i in range(5):
    for j in range(5):
        for n in range(M):
            s=[Fr(0)]*5
            for a in range(5):
                b=(3-a)%5
                for k in range(n+1):
                    s=cadd(s,cmul(Rh(a,i,k,-1),Rh(b,j,n-k)))
            s=[x/5 for x in s]
            target=cz(0,Fr(1)) if (i==j and n==0) else [Fr(0)]*5
            if not ceq(s,target): ok=False; print("fail",i,j,n,cnorm(s))
print("Rhat(-t)^T eta^-1 Rhat(t) = 1 through order",M-1,":",ok)
# tau check
tau_inv=[Fr(4,5*n+4) for n in range(8)]
def inv(a):
    b=[1/a[0]]
    for n in range(1,len(a)):
        b.append(-sum(a[k]*b[n-k] for k in range(1,n+1))/a[0])
    return b
print("tau",inv(tau_inv)[:4])
# c_g
B={2:Fr(-1,30),3:Fr(1,42),4:Fr(-1,30),5:Fr(5,66)}
for g in range(2,6):
    print(g,(-1)**(g-1)*B[g]*Fr(1,5**(g-1))/(2*g*(2*g-2)))
