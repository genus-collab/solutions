# Evidence only. Second-order test separating R(z), R(-z) and R^*: Int_{M_{0,5}} (R.omega)_{0,5}(H^5), by PPZ Def 2.2, 2.13,
# with R rebuilt from rows v_i = e_i^* o R^-1(-z) solving the row ODE with constants h_{wt_i}.
# True value: theta^2 <H,H,H> = Sum_{d>=1} d^2 Q^d/(e1 e2).
from fractions import Fraction as Fr
N=6
def S(c=0): return [Fr(c)]+[Fr(0)]*N
def add(*xs): return [sum(t) for t in zip(*xs)]
def sub(x,y): return [a-b for a,b in zip(x,y)]
def sc(c,x): return [c*a for a in x]
def mul(*xs):
    r=xs[0]
    for y in xs[1:]:
        o=[Fr(0)]*(N+1)
        for i,a in enumerate(r):
            if a:
                for j in range(N+1-i): o[i+j]+=a*y[j]
        r=o
    return r
def inv(x):
    r=[Fr(0)]*(N+1); r[0]=1/x[0]
    for n in range(1,N+1): r[n]=-sum(x[k]*r[n-k] for k in range(1,n+1))/x[0]
    return r
def th(x): return [k*a for k,a in enumerate(x)]
def thinv(x,c0): return [c0]+[x[k]/k for k in range(1,N+1)]
def mm(A,B): return [[add(*[mul(A[i][k],B[k][j]) for k in range(2)]) for j in range(2)] for i in range(2)]
def tr(A): return [[A[j][i] for j in range(2)] for i in range(2)]
def msc(c,A): return [[sc(c,A[i][j]) for j in range(2)] for i in range(2)]
def madd(A,B): return [[add(A[i][j],B[i][j]) for j in range(2)] for i in range(2)]
def mv(A,x): return [add(mul(A[0][0],x[0]),mul(A[0][1],x[1])),add(mul(A[1][0],x[0]),mul(A[1][1],x[1]))]
def run(a0,a1,a2,e1,e2):
    a0,a1,a2,e1,e2=map(Fr,(a0,a1,a2,e1,e2))
    geo=[Fr(1)]*(N+1); Qs=[Fr(0),Fr(1)]+[Fr(0)]*(N-1)
    s1=sc(Fr(-1),mul(add(S(a0),sc(a1+a2,Qs)),geo)); s2=sc(-a1*a2,mul(Qs,geo))
    def root(l0):
        L=S(l0)
        for _ in range(N+2):
            f=add(sub(mul(L,L),mul(s1,L)),s2); L=sub(L,mul(f,inv(sub(sc(Fr(2),L),s1))))
        return L
    Ls=[root(Fr(0)),root(-a0)]
    sig=a0+a1+a2; ee=e1*e2
    den=(a0+a1)*(a0+a2)*ee
    eta=[[sig/(a1*a2*den),1/den],[1/den,-a0/den]]
    det=eta[0][0]*eta[1][1]-eta[0][1]**2
    ei=[[S(eta[1][1]/det),S(-eta[0][1]/det)],[S(-eta[0][1]/det),S(eta[0][0]/det)]]
    Ds=[sc(ee,sub(add(S(a0*a1*a2),sc(2*a1*a2,L)),sc(sig,mul(L,L)))) for L in Ls]
    wts=[(a0,a1,a2,-e1,-e2),(-a0,a0+a1,a0+a2,-e1,-e2)]
    def hc(v):   # [z^1],[z^2] of h_v
        b=-Fr(1,12)*sum(1/x for x in v); return [Fr(1),b,b*b/2]
    cols=[[None,None],[None,None]]
    Rm=[None,None,None]
    colsn=[]
    for i in (0,1):
        L=Ls[i];Lj=Ls[1-i];D=Ds[i]
        g=sc(Fr(1,2),mul(th(D),inv(D)))
        Dop=lambda x: sub(th(x),mul(g,x))
        iv=inv(sub(L,Lj)); h=hc(wts[i])
        c=[S(1)]
        for n in (1,2):
            c.append(thinv(sc(Fr(-1),mul(Dop(Dop(c[n-1])),iv)),h[n]))
        r=[mul(L,c[0])]+[add(mul(L,c[n]),Dop(c[n-1])) for n in (1,2)]
        iD=inv(D)
        colsn.append([mv(ei,[mul(iD,c[n]),mul(iD,r[n])]) for n in range(3)])  # R_n e_i in {1,H}
    # R_n = Sum_i (R_n e_i) (x) e_i^*: matrix entries R[a][b] = Sum_i col_i[a] e_i^*(basis b), e_i^*(1)=1, e_i^*(H)=L_i
    for n in range(3):
        Rm[n]=[[add(colsn[0][n][a],colsn[1][n][a]) if b==0 else add(mul(colsn[0][n][a],Ls[0]),mul(colsn[1][n][a],Ls[1])) for b in (0,1)] for a in (0,1)]
    R1,R2=Rm[1],Rm[2]
    assert Rm[0]==[[S(1),S(0)],[S(0),S(1)]]
    etai=ei
    def total(B1,B2):
        # legs B(psi)=1+B1 psi+B2 psi^2, T=psi(1-B(psi)1), edge via the full formula, not assuming symplecticity
        one=[S(1),S(0)];H=[S(0),S(1)]
        est=lambda i,x: add(x[0],mul(Ls[i],x[1]))
        di=[inv(Ds[0]),inv(Ds[1])]
        h=[Ls[0],Ls[1]]
        a=[est(i,mv(B1,H)) for i in (0,1)]; b=[est(i,mv(B2,H)) for i in (0,1)]
        t2=[sc(Fr(-1),est(i,mv(B1,one))) for i in (0,1)]; t3=[sc(Fr(-1),est(i,mv(B2,one))) for i in (0,1)]
        # edge: E(x,y)= -B1 ei - x B2 ei - y ei B2^T holds under symplecticity, which is asserted separately for true R
        E0m=msc(Fr(-1),mm(B1,etai)); E1y=msc(Fr(-1),mm(etai,tr(B2)))
        ee2=lambda M,i,j: add(M[0][0],mul(M[0][1],Ls[j]),mul(Ls[i],M[1][0]),mul(Ls[i],Ls[j],M[1][1]))
        tot=S(0)
        for i in (0,1):
            hi=h[i]
            v=add(sc(Fr(5),mul(b[i],hi,hi,hi,hi)),sc(Fr(20),mul(a[i],a[i],hi,hi,hi)),mul(t3[i],hi,hi,hi,hi,hi),
                  sc(Fr(15),mul(a[i],t2[i],hi,hi,hi,hi)),sc(Fr(3),mul(t2[i],t2[i],hi,hi,hi,hi,hi)))
            tot=add(tot,mul(di[i],v))
        for i in (0,1):
            for j in (0,1):
                hi,hj=h[i],h[j]
                w=add(mul(ee2(E1y,i,j),hj,hj,hj),mul(ee2(E0m,i,j),add(sc(Fr(3),mul(a[j],hj,hj)),mul(t2[j],hj,hj,hj))))
                tot=add(tot,sc(Fr(10),mul(di[i],hi,hi,di[j],w)))
                for k in (0,1):
                    tot=add(tot,sc(Fr(15),mul(di[i],hi,hi,ee2(E0m,i,j),di[j],hj,ee2(E0m,j,k),di[k],h[k],h[k])))
        return tot
    expect=[Fr(d*d)/ee for d in range(N+1)]
    B1=msc(Fr(-1),R1); B2=madd(mm(R1,R1),msc(Fr(-1),R2))
    # symplecticity of B wrt eta: B2 ei + ei B2^T = B1 ei B1^T
    lhs=madd(mm(B2,etai),mm(etai,tr(B2))); rhs=mm(mm(B1,etai),tr(B1)); assert lhs==rhs
    print('weights',[str(x) for x in (a0,a1,a2,e1,e2)])
    print('  legs R^-1(psi)  [PPZ, claim]:', total(B1,B2)==expect)
    # z -> -z: legs R^-1(-psi)
    print('  legs R^-1(-psi) [z->-z]    :', total(msc(Fr(-1),B1),B2)==expect)
    # R^* in place of R: legs (R^*)^-1(psi)=R(-psi), B1'=-R1, B2'=R2
    print('  legs R(-psi)    [R->R^*]   :', total(B1,R2)==expect)
    print('  legs R(psi)     [R->R^-1]  :', total(R1,R2)==expect)
    print('  value claim, Q^0..Q^3      :',[str(x) for x in total(B1,B2)[:4]],' expected',[str(x) for x in expect[:4]])
    print('  value R->R^*,  Q^0..Q^2    :',[str(x) for x in total(B1,R2)[:3]])
for w in ((2,3,5,7,11),(Fr(1,3),Fr(-2,5),7,Fr(3,2),-1),(5,1,Fr(1,2),2,3)): run(*w)
