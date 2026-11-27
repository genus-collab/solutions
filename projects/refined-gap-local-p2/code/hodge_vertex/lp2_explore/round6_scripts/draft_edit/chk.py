from fractions import Fraction as Fr
import math, cmath
# r_k recurrence, signs
r=[Fr(1)]; rm=Fr(0)
for m in range(0,210):
    nxt=(3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>0 else 0))/(27*(m+1))
    r.append(nxt)
print([str(x) for x in r[:5]])
print(''.join('+' if x>0 else ('-' if x<0 else '0') for x in r[:60]))
print('zeros up to 200', sum(1 for x in r[:201] if x==0))
# Lemma 1.2(c) float: t_cf with log y = log|y| +- i pi
# I_1^[2](y) and I_0^[1](y): theta I_1^[2] = A log y + h, I_0^[1]: per Prop 2.12 Step 0
def a(d): return math.comb(3*d,d)*math.comb(2*d,d)
def H(n): return sum(1.0/k for k in range(1,n+1))
def series(y,N=4000):
    # A=sum a_d (-y)^d
    A=0;h=0;ser_log=0;ser2=0
    p=1.0
    for d in range(0,N):
        ad=math.exp(math.lgamma(3*d+1)-3*math.lgamma(d+1)) * ((-y)**d)
        A+=ad
        if d>=1:
            h+=3*ad*(H(3*d)-H(d))
    return A,h
print('r ok')
