from fractions import Fraction as Fr
N=12
a,b=Fr(1,3),Fr(2,3)
wc=[Fr(1)]
for n in range(N+2): wc.append(wc[-1]*(a+n)*(b+n)/((n+1)**2))
def mul(p,q,M=N):
    r=[Fr(0)]*M
    for i,x in enumerate(p[:M]):
        for j,y in enumerate(q[:M]):
            if i+j<M:r[i+j]+=x*y
    return r
def inv(p,M=N):
    r=[Fr(0)]*M;r[0]=1/p[0]
    for n in range(1,M):
        r[n]=-sum(p[k]*r[n-k] for k in range(1,n+1) if k<len(p))/p[0]
    return r
dw=[(n+1)*wc[n+1] for n in range(N+1)]
F=mul(mul([Fr(-1),Fr(1)],dw),inv(wc))
dt=[sum(wc[:n+1]) for n in range(N)]
tc=[Fr(0)]+[dt[n]/(n+1) for n in range(N-1)]
tt=tc[1:]            # t_c = xi * tt
# Laurent series in xi: dict power->coef
def L(*pairs): return dict(pairs)
def lmul(A,B,mx=4):
    R={}
    for i,x in A.items():
        for j,y in B.items():
            if i+j<=mx: R[i+j]=R.get(i+j,0)+x*y
    return R
def ladd(*Ls):
    R={}
    for A in Ls:
        for i,x in A.items(): R[i]=R.get(i,0)+x
    return R
def sc(c,A): return {i:c*x for i,x in A.items()}
Fs={i:F[i] for i in range(N)}
X={-1:Fr(1)}
Xinv={1:Fr(1)}
one={0:Fr(1)}
def p(A,n):
    R=one
    for _ in range(n): R=lmul(R,A,8)
    return R
F2=p(Fs,2);F3=p(Fs,3)
# BS eq:F2ref
F02=ladd(sc(Fr(5,8),lmul(F3,Xinv,8)), lmul(ladd(sc(Fr(5,8),Xinv),{0:Fr(-1,2)}),F2,8),
         sc(Fr(1,96),lmul(ladd(sc(13,X),sc(20,Xinv),{0:Fr(-32)}),Fs,8)),
         sc(Fr(-1,2160),ladd(sc(27,p(X,2)),sc(-98,X),{0:Fr(121)},sc(-50,Xinv))))
F11=ladd(sc(Fr(1,16),F2),sc(Fr(1,96),lmul(ladd({0:Fr(8)},sc(-7,X)),Fs,8)),sc(Fr(1,5760),ladd(sc(84,p(X,2)),sc(-196,X),{0:Fr(121)})))
F20=ladd(sc(Fr(1,384),lmul(X,Fs,8)),sc(Fr(-1,5760),ladd(sc(21,p(X,2)),sc(-29,X),{0:Fr(8)})))
# t_c^{-2}
ti=inv(tt); 
ti2=mul(ti,ti)
tm2={i-2:ti2[i] for i in range(N)}
bk={(2,0):Fr(-7,1920),(1,1):Fr(7,480),(0,2):Fr(-1,80)}
for name,S,key in (('F20',F20,(2,0)),('F11',F11,(1,1)),('F02',F02,(0,2))):
    D=ladd(S,sc(-bk[key],tm2))
    print(name,'t_c^-2 coef',bk[key],'diff coefs', [ (i,D.get(i,0)) for i in (-2,-1,0)])
# total at (3,7)
e1,e2=Fr(3),Fr(7);E=e1+e2;P=e1*e2
beta2=(7*e1**4+10*e1**2*e2**2+7*e2**4)/(5760*e1*e2)
tot=ladd(sc(E**4/P,F20),sc(E**2,F11),sc(P,F02),sc(3*beta2,tm2))
print('total diff',[(i,tot.get(i,0)) for i in (-2,-1,0)],'quoted -988747/4898880')
