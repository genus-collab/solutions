from fractions import Fraction as F
import cmath
# 1. c_g values
def bern(n):
    B=[F(1)]
    from math import comb
    for m in range(1,n+1):
        B.append(-sum(comb(m+1,k)*B[k] for k in range(m))/F(m+1))
    return B
B=bern(12)
for g in range(2,6):
    cg = (-1)**(g-1)*B[2*g]*F(1,5**(g-1))/(2*g*(2*g-2))
    print("c_%d ="%g, cg, " |B|form:", abs(B[2*g])/(2*g*(2*g-2)*5**(g-1)))
# 2. mu(b,k) table
nu=[0,-1,3,2,1]
for k in range(0,7):
    row=[]
    for b in range(5):
        e=4*k+nu[b]-nu[(b+k)%5]
        assert e%5==0
        row.append(e//5)
    print("k=%d mu="%k,row, " bound ok:", all(m>= max(0,(4*k-4)/5) for m in row))
# 3. Psi_Q^T eta^{-1} Psi_Q = 1 numerically with random m_a satisfying m_a m_{3-a}=1
import random
z=cmath.exp(2j*cmath.pi/5)
C0,C1,L=1.3,0.7,1.9
m=[L/C0, L**2/(C0*C1), C0*C1/L**2, C0/L, 1.0]
Psi=[[z**(j*(a+1))*m[a] for j in range(5)] for a in range(5)]
etainv=[[ (0.2 if (a+b)%5==3 else 0) for b in range(5)] for a in range(5)]
def mm(A,Bm): return [[sum(A[i][k]*Bm[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
def tr(A): return [[A[j][i] for j in range(5)] for i in range(5)]
P=mm(mm(tr(Psi),etainv),Psi)
print("max dev from I:", max(abs(P[i][j]-(1 if i==j else 0)) for i in range(5) for j in range(5)))
# unit: sum_j Delta_j^{-1/2} Psi_{aj} = 5 delta_{a3}; Delta_j^{-1/2} = zeta^j L/C0
print("unit image:", [round(abs(sum(z**j*L/C0*Psi[a][j] for j in range(5))),10) for a in range(5)])
# 4. iterate forward relation at z^0 to confirm m_a: X0 = L/C0, X_{a+1} = L X_a/(C1 kappa_a)
ck=[C1, L**5/(C0**2*C1**2), C1, C0, C0]
X=[L/C0]
for a in range(4): X.append(L*X[a]/ck[a])
print("X(0) vs m:", [round(X[a]-m[a],12) for a in range(5)], " closing:", round(L*X[4]-C0*X[0],12))
# 5. c(z)c(-z) coefficient of z^10
m5,m10=F(3,7),F(2,11)
print("z^10 coeff:", 2*m10-m5**2, " == -(m5^2-2m10):", -(m5**2-2*m10))
