# Independent numeric check (evidence only) of Lemma t_cf: floating point, stdlib only.
import math
from fractions import Fraction as Fr
def a(d): return math.factorial(3*d)//math.factorial(d)**3
H=[Fr(0)]
for n in range(1,4000): H.append(H[-1]+Fr(1,n))
def series(y,N=400):
    # log-free parts: g=I0^[1]-log y; h=theta I1^[2]-A log y; g2 from theta g2 = h-g
    x=-y
    A=1.0;g=0.0;g2=0.0;h=0.0
    term=1.0  # a_d x^d
    for d in range(1,N):
        term=term*x*(3*d)*(3*d-1)*(3*d-2)/d**3
        A+=term
        g_d=term/d                # theta g = A-1 => g_d = a_d x^d/d
        h_d=3*term*float(H[3*d]-H[d])
        g+=g_d; h+=h_d
        g2+=(h_d-g_d)/d           # theta g2 = h - g
    return A,g,g2,h
def tcf(y):
    A,g,g2,h=series(y)
    l=math.log(-y)
    return math.sqrt(3)/(2*math.pi)*(l*l/2+g*l+g2-math.pi**2/2),math.sqrt(3)/(2*math.pi)*(A*l+h)
# t_c from Taylor at xi=0: (xi-1) t' = -w, w=2F1(1/3,2/3;1;xi)
def tc(xi,N=2000):
    w=1.0;c=[1.0]
    for n in range(N):
        c.append(c[-1]*(n+1/3)*(n+2/3)/(n+1)**2)
    # t'(xi)=w/(1-xi)= sum_n (sum_{k<=n} c_k) xi^n
    s=0.0;t=0.0;xp=xi
    for n in range(N):
        s+=c[n]
        t+=s*xp/(n+1)
        xp*=xi
    w=sum(c[k]*xi**k for k in range(N))
    return t,w
for y in [-1/54,-1/81,-1/270,-1/40]:
    xi=1+27*y
    a1,a2=tcf(y); b1,b2=tc(xi)
    print(y,a1,b1,a1-b1,'theta:',a2,-b2,a2+b2)
