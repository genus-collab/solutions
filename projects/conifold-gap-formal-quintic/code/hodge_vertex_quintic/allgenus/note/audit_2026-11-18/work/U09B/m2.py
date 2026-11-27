import cmath
z=cmath.exp(2j*cmath.pi/5)
eps=[1,-1,-1,1,1]
F=[[1,-3/20,9/800],[-1,-1/20,111/800],[-1,-9/20,63/800],[1,1/20,81/800],[1,-3/20,21/160]]
def M(n):
    # [t^n] M0: (Psi0^{-1} Rhat)_{ki}, Rhat^{(i)}_b = z^{i(b+1)} Fs_b(z^{-i} t)
    return [[sum(eps[b]*z**(-k*(b+1))*z**(i*(b+1))*z**(-i*n)*F[b][n] for b in range(5))/5 for i in range(5)] for k in range(5)]
for n in range(3):
    A=M(n)
    print(n,"max|M-M^T|",max(abs(A[k][i]-A[i][k]) for i in range(5) for k in range(5)))
M1=M(1); M2=M(2)
i,j=0,1
print("M1_ij",M1[i][j])
print("correct (upper=row): psi1 coeff -M2_ji =",-M2[j][i]," psi2 coeff -M2_ij =",-M2[i][j])
print("transposed: psi1 coeff -M2_ij =",-M2[i][j]," psi2 coeff -M2_ji =",-M2[j][i])
