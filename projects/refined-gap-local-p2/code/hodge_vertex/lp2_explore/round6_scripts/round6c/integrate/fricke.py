# Check of the s.9 erratum on BS 5163-5167 (floating, evidence only for the wording).
import cmath, math
def eta(tau,N=400):
    q=cmath.exp(2j*math.pi*tau); p=cmath.exp(2j*math.pi*tau/24)
    for n in range(1,N): p*= (1-q**n)
    return p
def y_of(tau):
    a=eta(3*tau)**12; b=eta(tau)**12
    return -a/(b+27*a)
def F21(u,N=20000):
    s=0.0; c=1.0
    for n in range(N):
        s+=c*u**n
        c*= (n+1/3)*(n+2/3)/(n+1)**2
    return s
for s in (0.45,0.6,0.8):
    tau=1j*s
    tf=-1/(3*tau)
    y=y_of(tau); yf=y_of(tf)
    u=-27*y; uf=-27*yf
    X=1/(1+27*y); Xf=1/(1+27*yf)
    A=F21(u.real); wc=F21(1-u.real)
    print('s',s,'u',u.real,'u_F',uf.real,'1-u',1-u.real,'X_F',Xf.real,'1/u',1/u.real,'X',X.real)
    print('  A(tau)',A,'A(Fricke)',F21(uf.real),'wc',wc,'A_cf=i sqrt3 tau A',(1j*math.sqrt(3)*tau*A).real,'-wc',-wc)
