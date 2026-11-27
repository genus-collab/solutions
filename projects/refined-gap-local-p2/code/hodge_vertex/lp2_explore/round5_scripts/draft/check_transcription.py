# Evidence only. Re-evaluates the formulas as TRANSCRIBED in s.7.8 of the draft, in the notation of the draft (a_1=-tau, a_2=1+tau, kappa_0, S, nabla).
from fractions import Fraction as Fr
from math import factorial, comb
class D:
    __slots__=('x','d')
    def __init__(s,x,d=0): s.x=Fr(x); s.d=Fr(d)
    def __add__(s,o): o=o if isinstance(o,D) else D(o); return D(s.x+o.x,s.d+o.d)
    __radd__=__add__
    def __neg__(s): return D(-s.x,-s.d)
    def __sub__(s,o): o=o if isinstance(o,D) else D(o); return D(s.x-o.x,s.d-o.d)
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o): o=o if isinstance(o,D) else D(o); return D(s.x*o.x,s.x*o.d+s.d*o.x)
    __rmul__=__mul__
    def inv(s): return D(1/s.x,-s.d/s.x**2)
    def __truediv__(s,o): o=o if isinstance(o,D) else D(o); return s*o.inv()
    def __rtruediv__(s,o): return D(o)*s.inv()
    def __pow__(s,n):
        r=D(1)
        for _ in range(n): r=r*s
        return r
    def iszero(s): return s.x==0 and s.d==0
def Scoeffs(K):  # (y/2)/sinh(y/2) = sum s_k y^(2k)
    sh=[Fr(1,(2*k+1)*4**k*factorial(2*k+1)//(2*k+1)) for k in range(K+1)]  # sinh(y/2)/(y/2)= sum (y/2)^{2k}/(2k+1)!
    sh=[Fr(1,4**k*factorial(2*k+1)) for k in range(K+1)]
    s=[Fr(0)]*(K+1); s[0]=Fr(1)
    for k in range(1,K+1): s[k]=-sum(sh[j]*s[k-j] for j in range(1,k+1))
    return s
def beta(a1,a2,n,S):
    # beta_n=(n-1)! [t^n] 1/(4 sinh sinh) = (n-1)! [t^(n+2)] S(a1 t)S(a2 t)/(a1 a2)
    m=(n+2)//2; tot=D(0)
    for k in range(m+1): tot=tot+S[k]*a1**(2*k)*S[m-k]*a2**(2*(m-k))
    return tot/(a1*a2)*factorial(n-1)
def check_E1(tau0,Hmax=16,T=30):
    tau=D(tau0,1); a1=-tau; a2=1+tau; sig=tau*a2; sigp=1+2*tau
    S=Scoeffs(Hmax+2)
    # Laurent series in 1/N: dict power p of N^p -> D ; truncated at N^-T
    B={}
    for h in range(2,Hmax+1): B[2-2*h]=beta(a1,a2,2*h-2,S)
    def shift(f):  # f(N+tau) expansion in 1/N  (only negative powers and small positive)
        out={}
        for p,c in f.items():
            if p>=0:
                for j in range(p+1):
                    q=p-j
                    if q>=-T: out[q]=out.get(q,D(0))+c*comb(p,j)*tau**j
            else:
                k=-p; j=0
                while p-j>=-T:
                    co=(-1)**j*comb(k+j-1,j)
                    out[p-j]=out.get(p-j,D(0))+c*co*tau**j; j+=1
        return out
    def add(f,g,c=1):
        o=dict(f)
        for p,v in g.items(): o[p]=o.get(p,D(0))+v*c
        return o
    def der(f): return {p-1:c*p for p,c in f.items() if p!=0}
    def mulMN(f,fac): # multiply by polynomial dict
        o={}
        for p,c in f.items():
            for q,e in fac.items():
                if p+q>=-T: o[p+q]=o.get(p+q,D(0))+c*e
        return o
    nabla=lambda f: add(shift(f),f,-1)
    # B_tau: tau derivative at fixed N: use the dual part
    Btau={p:D(c.d,0) for p,c in B.items()}   # only value needed; note: value part suffices for the comparison of value parts
    Bp=der(B)
    LHS=add(add({p:c*(-2) for p,c in nabla(Btau).items()}, {p:-c for p,c in add(Bp,shift(Bp)).items()}), mulMN(nabla(Bp),{1:D(2)/a2*(-1),0:tau/a2*(-1)}))
    # RHS
    ell={}
    for j in range(1,T+2): ell[-j]=(-1)**(j+1)*tau**j/D(j)
    MN={2:D(1),1:tau}
    invN={-1:D(1)}; invM=shift(invN)
    r1=mulMN(ell,add(MN,{0:D(0)-sigp/12}))
    RHS={p:-c for p,c in r1.items()}
    RHS=add(RHS,{1:tau,0:tau*tau/2})   # tau(M+N)/2 = tau N + tau^2/2
    RHS=add(RHS,add(invN,invM),0)
    co=tau*(1+2*sig)/24
    RHS=add(RHS,{p:-c*co for p,c in add(invN,invM).items()})
    RHS={p:c/(sig*sig) for p,c in RHS.items()}
    bad=[]
    # validity range: B known through N^{2-2Hmax}; everything is exact down to about N^{2-2Hmax}
    for p in range(1,2-2*Hmax+1,-1):
        l=LHS.get(p,D(0)); r=RHS.get(p,D(0))
        if l.x!=r.x: bad.append(p)
    return bad
for t in [Fr(3,11),Fr(-7,4),Fr(13,5),Fr(2,7)]:
    print('(E1) tau',t,'mismatches',check_E1(t))
# control: perturb coefficient
def Lop_fstar(tau0,K=12):
    # L f = 2P(f_tau/t - f_t/a2) + (2-P) f/a2 on f_* and f_sing formula; series in t with powers -2..K
    tau=D(tau0,1); a1=-tau;a2=1+tau;sig=tau*a2;sigp=1+2*tau
    S=Scoeffs(K+4)
    f={-2:D(1)/(a1*a2)}
    for m in range(1,K//2+2):
        tot=D(0)
        for k in range(m+1): tot=tot+S[k]*a1**(2*k)*S[m-k]*a2**(2*(m-k))
        f[2*m-2]=tot/(a1*a2)
    def P():
        o={}
        for j in range(1,K+6): o[j]=D(0)-(-tau)**j/D(factorial(j))
        return o
    Pser=P()
    ftau_over_t={p-1:D(c.d) for p,c in f.items()}
    ft={p-1:c*p for p,c in f.items() if p!=0}
    # value-only computation
    res={}
    def addto(p,v):
        res[p]=res.get(p,Fr(0))+v
    for p,c in Pser.items():
        for q,v in ftau_over_t.items(): addto(p+q,2*c.x*v.x)
        for q,v in ft.items(): addto(p+q,-2*c.x*v.x/a2.x)
    for q,v in f.items(): addto(q,2*v.x/a2.x)
    for p,c in Pser.items():
        for q,v in f.items(): addto(p+q,-c.x*v.x/a2.x)
    # also 𝔏 f_sing formula
    d0=(1+2*sig)/(24*sig)
    fs={-2:D(0)-1/sig,0:d0}
    res2={}
    def add2(p,v): res2[p]=res2.get(p,Fr(0))+v
    ftt={p-1:D(c.d) for p,c in fs.items()}; fts={p-1:c*p for p,c in fs.items() if p!=0}
    for p,c in Pser.items():
        for q,v in ftt.items(): add2(p+q,2*c.x*v.x)
        for q,v in fts.items(): add2(p+q,-2*c.x*v.x/a2.x)
    for q,v in fs.items(): add2(q,2*v.x/a2.x)
    for p,c in Pser.items():
        for q,v in fs.items(): add2(p+q,-c.x*v.x/a2.x)
    # claimed: (1/s^2)[2P/t^3 - tau(2-P)/t^2 - (s'/12)P/t + (tau(1+2s)/24)(2-P)]
    cl={}
    def add3(p,v): cl[p]=cl.get(p,Fr(0))+v
    s2=sig.x**2
    for p,c in Pser.items():
        add3(p-3,2*c.x/s2); add3(p-2,tau.x*c.x/s2); add3(p-1,-sigp.x/12*c.x/s2); add3(p,-(tau.x*(1+2*sig.x)/24)*c.x/s2)
    add3(-2,-2*tau.x/s2); add3(0,2*(tau.x*(1+2*sig.x)/24)/s2)
    b1=[p for p in sorted(res) if p<=K-2 and res[p]!=0]
    b2=[p for p in sorted(set(res2)|set(cl)) if p<=K-2 and res2.get(p,0)!=cl.get(p,0)]
    return b1,b2
for t in [Fr(3,11),Fr(-7,4)]:
    print('L f_*=0, L f_sing formula: tau',t,'nonzero/mismatch',Lop_fstar(t))
# (H0) Step 1 at genus 2, s=0, Omega=1 (the dilaton equation is all that is used), c=1: 2-2h = -2
W={(4,):Fr(1,1152),(2,3):Fr(29,5760),(2,2,2):Fr(7,240)}
Z=-W[(4,)]+W[(2,3)]-W[(2,2,2)]/6
t14=3*W[(4,)]; t123=4*W[(2,3)]; t1222=5*W[(2,2,2)]
d1q0=-W[(4,)]+(t14+2*W[(2,3)])-(2*t123+W[(2,2,2)])/2+t1222/6
print('(H0) step1, genus 2, trivial class: d1(q0) =',d1q0,'  -2 Z =',-2*Z, d1q0==-2*Z)
# Z_2 formula in s.7.7(f)
tau=Fr(3,11); sig=tau*(1+tau); s=1/sig
S=Scoeffs(4); a1=D(-tau);a2=D(1+tau)
print('Z_2 check',beta(a1,a2,2,S).x/sig==-(24+28*s+7*s*s)/5760)
