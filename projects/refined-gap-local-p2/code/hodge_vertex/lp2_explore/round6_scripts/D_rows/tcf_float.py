from math import log, pi, sqrt
hs=[0.0]
for k in range(1,120002): hs.append(hs[-1]+1.0/k)
def tcf(y,s,nmax=40000):
    L=complex(log(abs(y)),s*pi)
    I01=L; I12=L*L/2
    ady=1.0
    for d in range(1,nmax):
        ady=ady*3*(3*d-1)*(3*d-2)/d**2*(-y)
        term=ady/(3*d)
        I01+=3*term
        I12+=3*L*term+9*term*(hs[3*d-1]-hs[d])
    return I01,I12
M=400
wc=[1.0]
for n in range(M): wc.append((n*(n-1)+2*n+2/9)*wc[n]/((n+1)**2))
ov=[];s=0.0
for n in range(M+1):
    s+=wc[n]; ov.append(s)
def tc(x): return sum(ov[n]*x**(n+1)/(n+1) for n in range(M))
for y in [-1/54,-1/40,-1/100,-1/33,-1/29]:
    x=1+27*y
    for s_ in (1,-1):
        a,b=tcf(y,s_)
        val=sqrt(3)/(2*pi)*(b-1j*pi*a-pi**2)
        print(f"y={y:.5f} s={s_:+d} BS t_cf={val:.10f}  t_c(x)={tc(x):.10f}")
