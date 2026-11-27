# Evidence only. Independent of rank2_check.py. Tests the sign of z in R against genus-1 localisation:
# theta F_1 = Int_{M_{1,1}} (R.omega)_{1,1}(H) computed from the graph sum, with R from the Laplace rows,
# versus Q^1 coefficient from virtual localisation, by hand: F_{1,1} = (1/24)(2 - sig (e1+e2)/(e1 e2)).
from fractions import Fraction as Fr
N=8
def S(c=0): return [Fr(c)]+[Fr(0)]*N
def add(x,y): return [a+b for a,b in zip(x,y)]
def sub(x,y): return [a-b for a,b in zip(x,y)]
def sc(c,x): return [c*a for a in x]
def mul(x,y):
    r=[Fr(0)]*(N+1)
    for i,a in enumerate(x):
        if a:
            for j in range(N+1-i): r[i+j]+=a*y[j]
    return r
def inv(x):
    r=[Fr(0)]*(N+1); r[0]=1/x[0]
    for n in range(1,N+1): r[n]=-sum(x[k]*r[n-k] for k in range(1,n+1))/x[0]
    return r
def th(x): return [k*a for k,a in enumerate(x)]
def pw(x,n):
    r=S(1)
    for _ in range(n): r=mul(r,x)
    return r
def run(a0,a1,a2,e1,e2):
    a0,a1,a2,e1,e2=map(Fr,(a0,a1,a2,e1,e2))
    geo=[Fr(1)]*(N+1); Qs=[Fr(0),Fr(1)]+[Fr(0)]*(N-1)
    s1=sc(Fr(-1),mul(add(S(a0),sc(a1+a2,Qs)),geo)); s2=sc(-a1*a2,mul(Qs,geo))
    def root(l0):
        L=S(l0)
        for _ in range(N+2):
            f=add(sub(mul(L,L),mul(s1,L)),s2); L=sub(L,mul(f,inv(sub(sc(Fr(2),L),s1))))
        return L
    Lp,Lm=root(Fr(0)),root(-a0)
    sig=a0+a1+a2; ee=e1*e2
    Del=lambda L: sc(ee,sub(add(S(a0*a1*a2),sc(2*a1*a2,L)),sc(sig,mul(L,L))))
    hz1=lambda v: -Fr(1,12)*sum(1/x for x in v)   # [z^1] h_v
    nw=(a1,a0+a2,-e1,-e2)
    tot_log=S(0); tot_R=S(0)
    for Li,Lj in ((Lp,Lm),(Lm,Lp)):
        D=Del(Li); A=sc(1/(a0+a2),sub(Li,S(a2))); B=sc(-1/a1,Li)
        # W_k via dA=A(1+A), dB=B(1+B)
        f1=lambda X: add(X,mul(X,X))
        A1,B1=f1(A),f1(B)
        W2=sub(sc(-(a0+a2),A1),sc(a1,B1))
        # W3 = -(a0+a2) (1+2A) A(1+A) - a1 (1+2B) B(1+B)
        d2=lambda X,X1: mul(add(S(1),sc(Fr(2),X)),X1)
        A2,B2=d2(A,A1),d2(B,B1)
        W3=sub(sc(-(a0+a2),A2),sc(a1,B2))
        # W4: theta^3 of X: theta[(1+2X)X1] = 2 X1^2 + (1+2X)^2 X1
        d3=lambda X,X1: add(sc(Fr(2),mul(X1,X1)),mul(mul(add(S(1),sc(Fr(2),X)),add(S(1),sc(Fr(2),X))),X1))
        W4=sub(sc(-(a0+a2),d3(A,A1)),sc(a1,d3(B,B1)))
        i2=inv(W2)
        ln1=sub(sc(Fr(1,8),mul(W4,mul(i2,i2))),sc(Fr(5,24),mul(mul(W3,W3),mul(i2,mul(i2,i2)))))
        c1=add(S(hz1(nw)),ln1)
        g=sc(Fr(1,2),mul(th(D),inv(D)))
        Rii=sub(c1,mul(g,inv(sub(Li,Lj))))
        tot_log=add(tot_log,mul(th(D),inv(D)))
        tot_R=add(tot_R,mul(Li,Rii))
    good=add(sc(Fr(1,48),tot_log),sc(Fr(1,2),tot_R))
    flip=sub(sc(Fr(1,48),tot_log),sc(Fr(1,2),tot_R))
    loc=Fr(1,24)*(2-sig*(e1+e2)/(e1*e2))
    q0=(a0/24)*sum(1/x for x in (-a0,a0+a1,a0+a2,-e1,-e2))
    print('weights',[str(x) for x in (a0,a1,a2,e1,e2)])
    print('  Q^0: formula',good[0],'localisation',q0,'flip',flip[0])
    print('  Q^1: formula',good[1],'localisation',loc,'flip',flip[1])
    print('  formula Q^d, d=1..8:',[str(x) for x in good[1:]])
    # (13.3) at a rational point
for w in ((2,3,5,7,11),(Fr(1,3),Fr(-2,5),7,Fr(3,2),-1),(5,1,Fr(1,2),2,3)): run(*w)
# (13.3): PF e^{W/z} = theta_U(km e^{W/z}), multiply out at rational points; z-coefficients separately
def pf(a0,a1,a2,Q,U):
    w=-a1*Q*U/(1-Q*U); thw=-a1*Q*U/(1-Q*U)**2
    W1=w-a2-(a0+a2)*U/(1-U)
    # PF e^{W/z} = [z (Q-1) thw + (Q-1) w^2 - (a0+(a1+a2)Q) w + a1 a2 Q] e^{W/z}
    lz1=(Q-1)*thw; lz0=(Q-1)*w*w-(a0+(a1+a2)*Q)*w+a1*a2*Q
    # km = -z a1 Q (1-U)/(1-QU); theta_U km = z * (-a1 Q) * U * d/dU[(1-U)/(1-QU)]
    dd=(-(1-Q*U)+Q*(1-U))/(1-Q*U)**2
    rz1=-a1*Q*U*dd; rz0=-a1*Q*(1-U)/(1-Q*U)*W1
    return lz1==rz1 and lz0==rz0
print('(13.3)',[pf(Fr(2),Fr(3),Fr(5),Fr(1,7),Fr(3,11)),pf(Fr(1,3),Fr(-2,5),Fr(7),Fr(-2,9),Fr(5,4))])
