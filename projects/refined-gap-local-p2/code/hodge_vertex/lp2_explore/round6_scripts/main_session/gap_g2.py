# Main-session check, written from BS eq:F2ref alone (printed in main.tex 4205-4220) and independent of agents' code.
# GW_2 = Sum_{k+h=2} E^(2k) (e1 e2)^(h-1) F_{k,h}. Substitute F -> F_cf, X -> 1/xi, expand in xi.
from fractions import Fraction as Fr
N=14
def mul(a,b):
    c=[Fr(0)]*N
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<N: c[i+j]+=x*y
    return c
def inv(a):
    b=[Fr(0)]*N; b[0]=1/a[0]
    for n in range(1,N):
        b[n]=-sum(a[k]*b[n-k] for k in range(1,n+1))/a[0]
    return b
# wc = 2F1(1/3,2/3;1;xi); A_cf=-wc (theta acts by y d/dy=(xi-1) d/dxi)
w=[Fr(1)]
for n in range(N-1): w.append(w[-1]*(Fr(1,3)+n)*(Fr(2,3)+n)/(n+1)**2)
dw=[(k+1)*w[k+1] for k in range(N-1)]+[Fr(0)]
xim1=[Fr(-1),Fr(1)]+[Fr(0)]*(N-2)
Fcf=mul(mul(xim1,dw),inv(w))   # F_cf = (xi-1) wc'/wc
# t_c: theta t_c = -wc?  (xi-1) t' = -wc => t' = wc/(1-xi)
s=inv([Fr(1),Fr(-1)]+[Fr(0)]*(N-2))
dt=mul(w,s)
tc=[Fr(0)]+[dt[k]/(k+1) for k in range(N-1)]  # tc=xi*(...)
print('F_cf',[str(x) for x in Fcf[:4]],'t_c',[str(x) for x in tc[:4]])
one=[Fr(1)]+[Fr(0)]*(N-1)
xi=[Fr(0),Fr(1)]+[Fr(0)]*(N-2)
def p(a,n):
    r=one
    for _ in range(n): r=mul(r,a)
    return r
F=Fcf; F2=p(F,2); F3=p(F,3)
# represent Laurent series as xi^-2 * series: multiply everything by xi^2
def add(*l):
    c=[Fr(0)]*N
    for a in l:
        for i,x in enumerate(a): c[i]+=x
    return c
def sc(a,k): return [x*k for x in a]
def sh(a,n): return ([Fr(0)]*n+a)[:N]  # times xi^n
# X=1/xi; xi^2*X^m = xi^(2-m)
F02=add(sc(sh(F3,3),Fr(5,8)),sc(sh(F2,3),Fr(5,8)),sc(sh(F2,2),Fr(-1,2)),
        sc(sh(F,1),Fr(13,96)),sc(sh(F,3),Fr(20,96)),sc(sh(F,2),Fr(-32,96)),
        sc(one,Fr(-27,2160)),sc(sh(one,1),Fr(98,2160)),sc(sh(one,2),Fr(-121,2160)),sc(sh(one,3),Fr(50,2160)))
F11=add(sc(sh(F2,2),Fr(1,16)),sc(sh(F,2),Fr(8,96)),sc(sh(F,1),Fr(-7,96)),
        sc(one,Fr(84,5760)),sc(sh(one,1),Fr(-196,5760)),sc(sh(one,2),Fr(121,5760)))
F20=add(sc(sh(F,1),Fr(1,384)),sc(one,Fr(-21,5760)),sc(sh(one,1),Fr(29,5760)),sc(sh(one,2),Fr(-8,5760)))
# tc^-2 *xi^2:
tcn=inv(tc[1:]+[Fr(0)])  # xi/t_c
t2=mul(tcn,tcn)
# (-3) beta_2 = -3 (7e1^4+10e1^2e2^2+7e2^4)/(5760 e1e2).  In E^2, e1e2 basis: e1^4+e2^4=E^4-4E^2 P+2P^2; e1^2e2^2=P^2
# => 7E^4 -28 E^2 P + 24 P^2. So GW_2 polar = -3/5760/P*(7E^4-28E^2P+24P^2)/t^2
# Coefficient of E^4/P: F20; E^2: F11; P: F02.
tg={'F20':Fr(-3*7,5760),'F11':Fr(3*28,5760),'F02':Fr(-3*24,5760)}
for nm,S in (('F20',F20),('F11',F11),('F02',F02)):
    diff=add(S,sc(t2,-tg[nm]))
    print(nm,'polar xi^-2',diff[0],'xi^-1',diff[1],'  const term',diff[2],'  target',tg[nm])
