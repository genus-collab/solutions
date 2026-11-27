from fractions import Fraction as Fr
from math import factorial as f
# 1. F^T_2 at y^0, erratum 5.1 item 1: Sum_i W_i^0 (p3 - p1^3)/17280
def g2d0(al,ep):
    tot=Fr(0)
    for i in range(3):
        j,k=[x for x in range(3) if x!=i]
        wt=[3*al[i],al[j]-al[i],al[k]-al[i],-ep[0],-ep[1]]
        W=Fr(1)
        for w in wt: W*=w
        p1=sum(1/w for w in wt); p3=sum(1/w**3 for w in wt)
        tot+=W*(p3-p1**3)/17280
    return tot
print('g2d0',g2d0([Fr(1),Fr(2),Fr(5)],[Fr(3),Fr(7)]),'expected 19189/100800')
print('g2d0',g2d0([Fr(2),Fr(-3),Fr(7)],[Fr(1),Fr(4)]),'expected 158713/290304')
print('g2d0',g2d0([Fr(1),Fr(3),Fr(8)],[Fr(5,2),Fr(-3,7)]),'expected -10457963/29030400')
# 2. N(W)|_{y=0} at alpha=(1,2,5), expected -1440
al=[Fr(1),Fr(2),Fr(5)]
e1=sum(al);e2=al[0]*al[1]+al[0]*al[2]+al[1]*al[2];e3=al[0]*al[1]*al[2]
W=lambda Y:e1*Y**2+2*e2*Y+3*e3
print('NW',W(-al[0])*W(-al[1])*W(-al[2]), -e3*((al[0]-al[1])*(al[0]-al[2])*(al[1]-al[2]))**2)
# 3. GW_{g,0} = 9 c_g E^{2g-2} candidates, c_g = Int lambda_g lambda_{g-1} lambda_{g-2}
B={2:Fr(1,6),4:Fr(-1,30),6:Fr(1,42),8:Fr(-1,30)}
def c(g): return Fr(1,2*f(2*g-2))*abs(B[2*g-2])/(2*g-2)*abs(B[2*g])/(2*g)
print('9c2',9*c(2),'9c3',9*c(3),'9c4',9*c(4))
# 4. Taylor expansion of t_c, wc, F_cf at xi=0 (exact), check t_c = xi + 11 xi^2/18, F_cf(0)=-2/9, rho
N=8
a,b=Fr(1,3),Fr(2,3)
wc=[Fr(1)]
for n in range(N): wc.append(wc[-1]*(a+n)*(b+n)/((n+1)**2))
# d t_c/d xi = wc/(1-xi)
dt=[sum(wc[:n+1]) for n in range(N)]
tc=[Fr(0)]+[dt[n]/(n+1) for n in range(N)]
print('tc',tc[:4])
# F_cf = theta log wc = (xi-1) wc'/wc
dw=[(n+1)*wc[n+1] for n in range(N)]
def mul(p,q,M=N):
    r=[Fr(0)]*M
    for i,x in enumerate(p):
        for j,y in enumerate(q):
            if i+j<M:r[i+j]+=x*y
    return r
def inv(p,M=N):
    r=[Fr(0)]*M;r[0]=1/p[0]
    for n in range(1,M):
        r[n]=-sum(p[k]*r[n-k] for k in range(1,n+1) if k<len(p))/p[0]
    return r
F=mul(mul([Fr(-1),Fr(1)],dw),inv(wc))
print('Fcf',F[:3])
# 5. beta values and limits
def beta2(e1,e2): return (7*e1**4+10*e1**2*e2**2+7*e2**4)/(5760*e1*e2)
print('beta2',beta2(Fr(7,3),Fr(5,2)))
