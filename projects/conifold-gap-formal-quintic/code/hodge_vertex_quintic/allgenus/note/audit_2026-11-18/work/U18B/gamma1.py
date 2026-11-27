import cmath
zeta=cmath.exp(2j*cmath.pi/5)
def p(k,L): return sum((zeta**j-L)**(-k) for j in range(5))
def c1(L):
    C=-L*(1-L**5)/5
    f3=(-5*L**-2+p(2,L))/6
    f4=(5*L**-3+p(3,L))/12
    s1=(L**-1+p(1,L))/2
    s2=(-L**-2+p(2,L))/4
    b1=1/12
    return (s2+s1**2/2)*C+3*C**2*(s1*f3+f4)+7.5*f3**2*C**3+b1*(L**-1/5+p(1,L))
for eps in [1e-2,1e-3,1e-4,1e-5]:
    print(eps, c1(1+eps), c1(1+1j*eps))
print("zeta*L test c1(zeta L)=zeta^-1 c1(L):", c1(zeta*1.3), c1(1.3)/zeta)
print("L->0:", c1(1e-4), c1(1e-4j))
