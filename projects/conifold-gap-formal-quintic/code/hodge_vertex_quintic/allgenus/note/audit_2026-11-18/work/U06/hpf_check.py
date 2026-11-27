# Independent sanity check for unit U06 (lem:gau-poly, lem:gau-wkb, thm:gau-R0 (ii)-(iii)):
#  1. solve HPF recursively in Q[u] (free coefficients 0), check polynomiality and resonance at n=5,10
#  2. check directly that A_n = L^{4n+1} ph_n(L^{-5}) satisfies PF5^{(0)} A = 0 in Q[L,L^-1][[z]]
#  3. check T_1(L g) = 5 Dg and the z^0 coefficient of PF5^{(0)}
#  4. check lem:gau-poly (ii),(iii): Phi*k(z) solves HPF iff ...; free coefficient formula
from fractions import Fraction as Fr
N = 11  # orders of zhat

# ---------- polynomials in u: dict {k: Fr}
def padd(a,b,c=1):
    r=dict(a)
    for k,v in b.items():
        r[k]=r.get(k,0)+c*v
        if r[k]==0: del r[k]
    return r
def pmul(a,b):
    r={}
    for k,v in a.items():
        for l,w in b.items():
            r[k+l]=r.get(k+l,0)+v*w
    return {k:v for k,v in r.items() if v!=0}
def pscal(a,c):
    return {k:v*c for k,v in a.items() if v*c!=0}
def theta(a):  # u d/du
    return {k:v*k for k,v in a.items() if k!=0}
ONE={0:Fr(1)}
OMU={0:Fr(1),1:Fr(-1)}   # 1-u
def upow(n): return {n:Fr(1)}

def O(e,Phi):
    # (O_e Phi)_m = Phi_m + (1-u)((4(m-1)-e)/5 - theta)Phi_{m-1}
    out=[]
    for m in range(len(Phi)):
        t=dict(Phi[m])
        if m>=1:
            inner=padd(pscal(Phi[m-1],Fr(4*(m-1)-e,5)),theta(Phi[m-1]),-1)
            t=padd(t,pmul(OMU,inner))
        out.append(t)
    return out
def shift(Phi,j,pre):  # multiply by pre(u) * zhat^j
    out=[{} for _ in Phi]
    for m in range(len(Phi)-j):
        out[m+j]=pmul(pre,Phi[m])
    return out
def HPF(Phi):
    P=[Phi]
    for j in range(1,6):
        P.append(O(-j,P[-1]))
    res=[padd(P[5][m],Phi[m],-1) for m in range(len(Phi))]
    sig=[Fr(3),Fr(17,5),Fr(9,5),Fr(274,625),Fr(24,625)]
    for k in range(1,6):
        term=shift(P[5-k],k,pmul(OMU,pscal(upow(k-1),sig[k-1])))
        res=[padd(res[m],term[m],-1) for m in range(len(Phi))]
    return res
def divOMU(p):
    # exact division of polynomial by (1-u)
    if not p: return {}
    d=max(p); q={}; carry=Fr(0)
    # p = (1-u) q ; q_k = p_k + q_{k-1}
    for k in range(0,d+1):
        carry=p.get(k,Fr(0))+carry
        if carry!=0: q[k]=carry
    assert carry==0, "not divisible by 1-u"
    return q
def solve(free):
    Phi=[{} for _ in range(N+2)]
    Phi[0]=dict(ONE)
    rs={}
    for n in range(1,N+1):
        Phi[n]={}
        res=HPF(Phi)
        rn=divOMU(res[n+1])   # coefficient of zhat^{n+1} with ph_n=0 is (1-u) r_n
        rs[n]=rn
        ph={}
        for k,v in rn.items():
            if 4*n==5*k:
                raise AssertionError("resonance violated at n=%d: %s"%(n,v))
            ph[k]=-v/(4*n-5*k)
        if n%5==0 and free.get(n//5,0)!=0:
            ph[4*n//5]=Fr(free[n//5])
        Phi[n]=ph
    # verify
    res=HPF(Phi)
    assert all(not res[m] for m in range(N+1)), "HPF not solved"
    return Phi,rs

Phi,rs=solve({})
print("ph_1..ph_6 (free coefficients 0):")
for n in range(1,7):
    print(n, {k:str(v) for k,v in sorted(Phi[n].items())}, "deg r_n =",max(rs[n]) if rs[n] else None)
print("resonance [u^4] r_5 =",rs[5].get(4,0),"  [u^8] r_10 =",rs[10].get(8,0))
print("degrees of ph_n:",[max(Phi[n]) if Phi[n] else 0 for n in range(N+1)])

# ---------- Laurent polynomials in L, PF5^{(0)} directly
def D(a):  # q d/dq with DL=(L^6-L)/5
    r={}
    for k,v in a.items():
        if k!=0:
            r[k+5]=r.get(k+5,0)+v*Fr(k,5)
            r[k]=r.get(k,0)-v*Fr(k,5)
    return {k:v for k,v in r.items() if v!=0}
Lp={1:Fr(1)}
qL={0:Fr(1,5**5),-5:Fr(-1,5**5)}
def M(A,c=1,kz=0):
    # (c*(zD+L) + kz*z) applied to series A (list of Laurent polys), truncated
    out=[{} for _ in A]
    for n in range(len(A)):
        out[n]=padd(out[n],pscal(pmul(Lp,A[n]),Fr(c)))
        if n+1<len(A):
            out[n+1]=padd(out[n+1],pscal(D(A[n]),Fr(c)))
            if kz: out[n+1]=padd(out[n+1],pscal(A[n],Fr(kz)))
    return out
def PF5_0(A):
    X=A
    for _ in range(5): X=M(X)
    Y=A
    for k in range(1,6): Y=M(Y,5,k)
    return [padd(padd(X[n],A[n],-1),pmul(qL,Y[n]),-1) for n in range(len(A))]
def flat(Phi):
    return [{(4*n+1)-5*k:v for k,v in Phi[n].items()} for n in range(len(Phi))]
A=flat(Phi[:N+1])
res=PF5_0(A)
print("PF5^(0) A^can = 0 through z^%d:"%N, all(not res[n] for n in range(N+1)))
# T_1 check: coefficient of z^1 of PF5^(0) applied to (L*g) with g=L^3 (arbitrary), expected 5 Dg
g={3:Fr(1),-2:Fr(7)}
test=[pmul(Lp,g)]+[{} for _ in range(3)]
r=PF5_0(test)
print("z^0 coefficient vanishes:", not r[0], "  T_1(Lg)=5Dg:", r[1]==pscal(D(g),Fr(5)))

# ---------- lem:gau-poly (ii),(iii): multiply by k(z), z^5 = u^4 zhat^5
def mulk(Phi,k5):  # k = 1 + sum k5[m] z^{5m}
    out=[dict(p) for p in Phi]
    for m,c in k5.items():
        for n in range(len(Phi)-5*m):
            out[n+5*m]=padd(out[n+5*m],pmul(upow(4*m),pscal(Phi[n],Fr(c))))
    return out
Phik=mulk(Phi,{1:Fr(7,3),2:Fr(-2)})
res=HPF(Phik)
print("Phi*k(z) with k in Q[[z^5]] solves HPF:", all(not res[m] for m in range(N+1)))
f1=Phik[5].get(4,0); f2=Phik[10].get(8,0)
print("free coefficients of Phi*k:",f1,f2," expected k5, k10 + k5*f1(Phi)+f2(Phi) =",Fr(7,3),Fr(-2)+Fr(7,3)*Phi[5].get(4,0)+Phi[10].get(8,0))
Phi2,_=solve({1:Fr(7,3),2:f2})
print("solution with the same free coefficients equals Phi*k:",all(Phi2[n]==Phik[n] for n in range(N+1)))

# ---------- h(z) = h_{w_0}: sum_{a in w_0} a^{1-2k} = -(T_k), T_k = 5^{1-2k} + sum_{j=1..4}(zeta^j-1)^{1-2k}
# power sums of 1/(1-zeta^j): roots x_j=1/(1-zeta^j) satisfy ((x-1)/x)^5=1 -> (x-1)^5 - x^5 = 0 -> -5x^4+10x^3-10x^2+5x-1=0
# i.e. x^4 - 2x^3 + 2x^2 - x + 1/5 = 0 ; Newton sums
e=[Fr(1),Fr(2),Fr(2),Fr(1),Fr(1,5)]
p=[Fr(4)]
for m in range(1,2*N):
    s=Fr(0)
    for i in range(1,min(m,4)+1):
        s+=(-1)**(i-1)*e[i]*(p[m-i] if m-i>0 else (Fr(m) if m==i else 0)) if m-i>0 else 0
    if m<=4: s+=(-1)**(m-1)*m*e[m]
    p.append(s)
# Bernoulli numbers
from math import comb
B=[Fr(1)]
for m in range(1,2*N+2):
    B.append(-sum(comb(m+1,j)*B[j] for j in range(m))/(m+1))
# log h = - sum_k b_k z^{2k-1} * sum_{a in w0} a^{1-2k};  w0 = {1-zeta^j} U {-5}
logh=[Fr(0)]*(N+1)
for k in range(1,(N+1)//2+1):
    if 2*k-1<=N:
        bk=B[2*k]/(2*k*(2*k-1))
        S=p[2*k-1]+Fr(-5)**(1-2*k)
        logh[2*k-1]=-bk*S
# exp
h=[Fr(1)]+[Fr(0)]*N
for n in range(1,N+1):
    h[n]=sum(Fr(m)*logh[m]*h[n-m] for m in range(1,n+1))/n
phat1=[sum(Phi[n].values(),Fr(0)) for n in range(N+1)]
print("h_1..h_4      :",[str(x) for x in h[1:5]])
print("ph_n(1), n<5  :",[str(x) for x in phat1[1:5]])
