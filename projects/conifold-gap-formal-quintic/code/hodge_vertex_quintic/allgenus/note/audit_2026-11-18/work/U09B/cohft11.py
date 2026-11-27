import cmath
z=cmath.exp(2j*cmath.pi/5)
eps=[1,-1,-1,1,1]
F1=[-3/20,-1/20,-9/20,1/20,-3/20]   # [t^1] Fs_b, colour 0, from the T-system with a=-3/20
# Psi0_{bi}=eps_b z^{i(b+1)}; Psi0^{-1}_{kb}=eps_b z^{-k(b+1)}/5
R1=[[sum(eps[b]*z**(-k*(b+1))*z**(i*b)*F1[b] for b in range(5))/5 for i in range(5)] for k in range(5)]
sym=max(abs(R1[k][i]-R1[i][k]) for i in range(5) for k in range(5))
print("R1 symmetric defect",sym)
sig=[z**(-i) for i in range(5)]       # varsigma_i at v=1
e=[-z**(2*i) for i in range(5)]       # (Psi^m)^{-1} H = -sum zeta^{2i} e_i at v=1
one=[z**i for i in range(5)]
def mat(M,x): return [sum(M[k][i]*x[i] for i in range(5)) for k in range(5)]
R1e=mat(R1,e); R11=mat(R1,one)
leg=-sum(sig[i]*R1e[i] for i in range(5))/24
dil=sum(sig[i]**2*e[i]*R11[i] for i in range(5))/24
loop=sum(R1[i][i]*sig[i]*e[i] for i in range(5))/2
print("leg",leg,"dil",dil,"loop",loop)
print("<H>_1 with R=M0(z):",leg+dil+loop)
print("<psi 1>_1 =", sum(sig[i]**-1*sig[i] for i in range(5))/24)
