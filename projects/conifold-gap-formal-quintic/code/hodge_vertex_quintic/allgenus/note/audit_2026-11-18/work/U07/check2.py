# Check Lemma E, Corollary E-prime formula, Lemma T (a),(b),(d) on a random graded symplectic dressing with constant coefficients.
from fractions import Fraction as Fr
from math import factorial
import random
random.seed(11)
nu=[0,-1,3,2,1]; eps=[1,-1,-1,1,1]
NV=6       # modulo v^NV
K=16       # z-truncation for matrix series: powers 0..K-1

def M0(): return [[Fr(0)]*5 for _ in range(5)]
def I5():
    m=M0()
    for i in range(5): m[i][i]=Fr(1)
    return m
def mmul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
def madd(A,B): return [[A[i][j]+B[i][j] for j in range(5)] for i in range(5)]
def msc(c,A): return [[c*A[i][j] for j in range(5)] for i in range(5)]
def mT(A): return [[A[j][i] for j in range(5)] for i in range(5)]
eta=M0()
for a in range(5): eta[a][(3-a)%5]=Fr(5)
etainv=msc(Fr(1,25),eta)
assert mmul(eta,etainv)==I5()

# matrix series: list of K matrices
def smul(A,B):
    C=[M0() for _ in range(K)]
    for i in range(K):
        for j in range(K-i):
            C[i+j]=madd(C[i+j],mmul(A[i],B[j]))
    return C
def sone():
    S=[M0() for _ in range(K)]; S[0]=I5(); return S
def sexp(r):   # r has no constant term
    S=sone(); term=sone()
    for n in range(1,K):
        term=smul(term,r); term=[msc(Fr(1,n),t) for t in term]
        S=[madd(S[i],term[i]) for i in range(K)]
    return S
def sinv(A):   # A[0] invertible diagonal
    A0inv=M0()
    for i in range(5): A0inv[i][i]=1/A[0][i][i]
    B=[M0() for _ in range(K)]; B[0]=A0inv
    for n in range(1,K):
        acc=M0()
        for j in range(1,n+1): acc=madd(acc,mmul(A[j],B[n-j]))
        B[n]=msc(Fr(-1),mmul(A0inv,acc))
    return B
def sneg(A): return [msc(Fr((-1)**k),A[k]) for k in range(K)]   # z -> -z
def sT(A): return [mT(a) for a in A]

# random graded r_k = sum_b rho_{k,b} E_{b+k,b} with eta r_k symmetric (k odd) / antisymmetric (k even)
def rand_r():
    r=[M0() for _ in range(K)]
    for k in range(1,K):
        rho=[None]*5
        for b in range(5):
            if rho[b] is not None: continue
            a=(3-k-b)%5
            val=Fr(random.randint(-3,3),random.randint(1,4))
            if k%2==1:
                rho[b]=val; rho[a]=val
            else:
                if a==b: rho[b]=Fr(0)
                else: rho[b]=val; rho[a]=-val
        for b in range(5): r[k][(b+k)%5][b]=rho[b]
        er=mmul(eta,r[k])
        assert er==msc(Fr((-1)**(k+1)),mT(er))
    return r
R=sexp(rand_r())          # lowered-component form, R = 1+O(z), graded
# symplectic check
S=smul(sT(sneg(R)),smul([eta]+[M0()]*(K-1),R))
assert S[0]==eta and all(S[k]==M0() for k in range(1,K))
# graded check
for k in range(K):
    for a in range(5):
        for b in range(5):
            if (a-b-k)%5!=0: assert R[k][a][b]==0
# P = P0 R with P0 diagonal, P0^T eta P0 = eta: p_a p_{3-a}=1
p0=[Fr(3,2),Fr(-2),Fr(-1,2),Fr(2,3),Fr(1)]
for a in range(5): assert p0[a]*p0[(3-a)%5]==1
P0=M0()
for a in range(5): P0[a][a]=p0[a]
P=[mmul(P0,R[k]) for k in range(K)]
d=p0[3]
def vec(A): return [mmul(etainv,mmul(a,eta)) for a in A]
Pvec=vec(P); Rvec=vec(R)
P0vec=mmul(etainv,mmul(P0,eta)); assert P0vec[0][0]==d
Rvi=sinv(Rvec); Pvi=sinv(Pvec)

# edge E_{Rvec}: numerator N(psi,psi2) = etainv - Rvi(psi) etainv Rvi(psi2)^T ; E = N/(psi+psi2)
def Ncoef(i,j):
    m=mmul(Rvi[i],mmul(etainv,mT(Rvi[j])))
    m=msc(Fr(-1),m)
    if i==0 and j==0: m=madd(m,etainv)
    return m
D=K-1
E={}
for tot in range(0,D):
    for j in range(0,tot+1):
        i=tot-j
        # e_{i,j} = n_{i+1,j} - e_{i+1,j-1}
        m=Ncoef(i+1,j)
        if j>0: m=madd(m,msc(Fr(-1),E[(i+1,j-1)]))
        E[(i,j)]=m
    # divisibility check: n_{0,tot+1} = e_{0,tot}
    assert Ncoef(0,tot+1)==E[(0,tot)],tot
# symmetry E(psi,psi2)^T = E(psi2,psi)
for (i,j),m in E.items(): assert mT(m)==E[(j,i)]
print("edge is a power series and symmetric: ok (symplecticity)")

# same edge from Pvec instead of Rvec
def NcoefP(i,j):
    m=mmul(Pvi[i],mmul(etainv,mT(Pvi[j])))
    m=msc(Fr(-1),m)
    if i==0 and j==0: m=madd(m,etainv)
    return m
for i in range(0,6):
    for j in range(0,6):
        assert NcoefP(i,j)==Ncoef(i,j)
print("numerator with Pvec equals numerator with Rvec: ok")

# ---- Lemma E with a random f in H_- (no u): f = sum_{m=1..4} z^{-m} f_m
f={m:[Fr(random.randint(-3,3)) for _ in range(5)] for m in range(1,5)}
def mv(A,x): return [sum(A[i][k]*x[k] for k in range(5)) for i in range(5)]
def etapair(x,y): return sum(x[a]*sum(eta[a][b]*y[b] for b in range(5)) for a in range(5))
LMAX=8
# contraction: pairs psi^m e with eta([z^{-m-1}] f(-z), e)
C=[[Fr(0)]*5 for _ in range(LMAX)]
for l in range(LMAX):
    for m in range(0,4):
        fm=[(-1)**(m+1)*x for x in f[m+1]]       # [z^{-m-1}] f(-z)
        for b in range(5):
            C[l][b]+=sum(etapair(fm,[Fr(1) if c==a else Fr(0) for c in range(5)])*E[(m,l)][a][b] for a in range(5))
# right side: -R^{-1}(psi2) (Rf)_{>=0}(psi2)
Rf={}
for k in range(K):
    for m in range(1,5):
        Rf[k-m]=[x+y for x,y in zip(Rf.get(k-m,[Fr(0)]*5),mv(Rvec[k],f[m]))]
rhs=[[Fr(0)]*5 for _ in range(LMAX)]
for l in range(LMAX):
    for j in range(0,l+1):
        rhs[l]=[x-y for x,y in zip(rhs[l],mv(Rvi[l-j],Rf.get(j,[Fr(0)]*5)))]
assert C==rhs
print("Lemma E identity ok to order psi2^%d"%(LMAX-1))

# ---- vanishing solution mod v^NV, Pi_i
def pz(): return [Fr(0)]*NV
def padd(a,b): return [x+y for x,y in zip(a,b)]
def pscal(c,a): return [c*x for x in a]
YC=[{} for _ in range(5)]
for b in range(5):
    m=nu[b]
    for i in range(0,30):
        j=m+4*i+1
        if j<0 or j>=NV: continue
        c=eps[b]*Fr((-1)**i,4**i*factorial(i))*Fr((-1)**j,factorial(j))
        n=i+j
        t=pz(); t[j]=c
        YC[b][n]=padd(YC[b].get(n,pz()),t)
nmax=max(max(Y.keys()) for Y in YC)
assert nmax+2<K
def Pi(i,q):   # q: dict (k,b)-> Fraction (v-independent)
    tot=pz()
    for (k,b),c in q.items():
        n=k-i
        if n<0 or n not in YC[b]: continue
        tot=padd(tot,pscal(-c*(-1)**n,YC[b][n]))
    return tot
# Yp = P YC : dict (power,component) -> poly in v
Yp={}
for k in range(K):
    for b in range(5):
        for n,c in YC[b].items():
            for a in range(5):
                if P[k][a][b]!=0:
                    key=(k-n,a)
                    Yp[key]=padd(Yp.get(key,pz()),pscal(P[k][a][b],c))
# leaf point q_T = Pvec^{-1}(-z 1)
qT={}
for k in range(K-1):
    for a in range(5):
        qT[(k+1,a)]=-Pvi[k][a][0]
# Lemma T (b): Pi_i(q_T) = (-1)^{i+1} [z^{i-1}] Yp_0   (no (A) needed for this formula)
for i in range(0,7):
    lhs=Pi(i,qT)
    rhs_=pscal(Fr((-1)**(i+1)),Yp.get((i-1,0),pz()))
    assert lhs==rhs_,i
print("Lemma T(b) formula Pi_i(q_T)=(-1)^{i+1}[z^{i-1}]Yp_0 ok for i<=6 (mod v^%d)"%NV)
# [z^0] Yp_0 = 0 by the grading
assert all(x==0 for x in Yp.get((0,0),pz()))
# Lemma T (d): T_P - p_{0,0} v in v^2
TP=pscal(Fr(-1),Yp.get((-1,0),pz()))
assert TP[0]==0 and TP[1]==p0[0] and d*p0[0]==1
print("Lemma T(d) ok: T_P =",TP[:4])
# leaf: q_T = -(z/d) 1 + O(z^2)
assert qT[(1,0)]==-1/d and all(qT[(1,a)]==0 for a in range(1,5))

# ---- Corollary E-prime: (Pi_{k-1} (x) 1) E_P  vs  (-1)^{k+1} * ( - Rvec^{-1}(psi2) z^{-k} P0vec^{-1} etainv (Yp)_{>=k} )
for kk in range(1,4):
    i=kk-1
    for l in range(0,5):
        for b in range(5):
            col={(m,a):E[(m,l)][a][b] for m in range(0,K-2-l) for a in range(5)}
            lhs=Pi(i,col)
            tot=pz()
            for j in range(0,l+1):
                for a in range(5):
                    for c in range(5):
                        coef=Rvi[l-j][b][a]*(1/P0vec[a][a])*etainv[a][c]
                        if coef!=0:
                            tot=padd(tot,pscal(coef,Yp.get((j+kk,c),pz())))
            rhs_=pscal(Fr(-1)*(-1)**(kk+1),tot)
            assert lhs==rhs_,(kk,l,b)
print("Corollary E-prime transport formula ok (k=1..3, psi2^0..4, mod v^%d)"%NV)
