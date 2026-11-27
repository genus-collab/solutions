# Sanity checks for unit U05 (lens A): Lemma gau-HPF (i)-(iii) to low order.
from fractions import Fraction as Fr
import random, cmath

# ---------- Laurent polynomials in s (dict exp->Fraction) ----------
def ladd(a,b,c=1):
    r=dict(a)
    for e,v in b.items():
        r[e]=r.get(e,0)+c*v
        if r[e]==0: del r[e]
    return r
def lmul(a,b):
    r={}
    for e1,v1 in a.items():
        for e2,v2 in b.items():
            r[e1+e2]=r.get(e1+e2,0)+v1*v2
    return {e:v for e,v in r.items() if v!=0}
def lscal(a,c): return {e:v*c for e,v in a.items() if v*c!=0}
def D(a):  # D s^c = (c/5)(s^c - s^(c-5))
    r={}
    for e,v in a.items():
        if e==0: continue
        r[e]=r.get(e,0)+Fr(e,5)*v
        r[e-5]=r.get(e-5,0)-Fr(e,5)*v
    return {e:v for e,v in r.items() if v!=0}
def theta_u(a): return {e:Fr(e,5)*v for e,v in a.items() if e!=0}
ONE_MINUS_U={0:Fr(1),5:Fr(-1)}
def spow(c): return {c:Fr(1)}

N=7
random.seed(1)
def randlaur():
    return {random.randint(-7,9):Fr(random.randint(-5,5),random.randint(1,4)) for _ in range(3)}
ph=[randlaur() for _ in range(N+1)]

# ---------- direct: PF5^(0)(L Phi) as series in z ----------
def Mop(F):  # (zD + L)F
    G=[]
    for m in range(N+1):
        t=lmul(spow(-1),F[m])
        if m>=1: t=ladd(t,D(F[m-1]))
        G.append(t)
    return G
def zshift(F,k): return [{} for _ in range(k)]+F[:N+1-k]
def sadd(F,G,c=1): return [ladd(F[m],G[m],c) for m in range(N+1)]
def sscal(F,p): return [lmul(p,F[m]) for m in range(N+1)]
LPhi=[lmul(spow(-4*n-1),ph[n]) for n in range(N+1)]
M5=LPhi
for _ in range(5): M5=Mop(M5)
prod=LPhi
for k in range(1,6):
    t=Mop(prod)
    prod=sadd([lscal(x,Fr(5)) for x in t], zshift(prod,1), Fr(k))
q={0:Fr(1,5**5),5:Fr(-1,5**5)}
direct=sadd(sadd(M5,LPhi,-1), sscal(prod,q), -1)

# ---------- HPF left side as series in zh ----------
def O(e,F):
    G=[]
    for m in range(N+1):
        t=dict(F[m])
        if m>=1:
            inner=ladd(lscal(F[m-1],Fr(4*(m-1)-e,5)), theta_u(F[m-1]), -1)
            t=ladd(t,lmul(ONE_MINUS_U,inner))
        G.append(t)
    return G
Phis=[ph]
for j in range(1,6): Phis.append(O(-j,Phis[-1]))
sig=[None,Fr(3),Fr(17,5),Fr(9,5),Fr(274,625),Fr(24,625)]
br=[{} for _ in range(N+1)]
for k in range(1,6):
    term=zshift(Phis[5-k],k)
    term=sscal(term,{5*(k-1):sig[k]})
    br=sadd(br,term)
lhs=sadd(sadd(Phis[5],ph,-1), sscal(br,ONE_MINUS_U), -1)
hpf_as_z=[lmul(spow(-4*n-1),lhs[n]) for n in range(N+1)]
print("(i) PF5^(0)(L Phi) == s^-1 * HPF-lhs to z^%d:"%N, all(direct[m]==hpf_as_z[m] for m in range(N+1)))

# ---------- (ii): coefficient of zh^(n+1) is (1-u)[(4n-5theta)ph_n + r_n] ----------
ok=True
for n in range(0,N):
    # Phi-flat: ph_n replaced by 0
    ph0=list(ph); ph0[n]={}
    P0=[ph0]
    for j in range(1,6): P0.append(O(-j,P0[-1]))
    def g(F,m): return F[m] if m>=0 else {}
    rn={}
    for j in range(1,6):
        x=P0[j-1][n]
        rn=ladd(rn, ladd(lscal(x,Fr(4*n+j,5)), theta_u(x), -1))
    rn=ladd(rn, P0[4][n], -3)
    rn=ladd(rn, lmul({5:Fr(17,5)}, g(Phis[3],n-1)), -1)
    rn=ladd(rn, lmul({10:Fr(9,5)}, g(Phis[2],n-2)), -1)
    rn=ladd(rn, lmul({15:Fr(274,625)}, g(Phis[1],n-3)), -1)
    rn=ladd(rn, lmul({20:Fr(24,625)}, g(ph,n-4)), -1)
    pred=lmul(ONE_MINUS_U, ladd(ladd(lscal(ph[n],Fr(4*n)), theta_u(ph[n]), -5), rn))
    ok = ok and (pred==lhs[n+1])
print("(ii) coefficient formula with r_n, n=0..%d:"%(N-1), ok, "; coefficient of zh^0 vanishes:", lhs[0]=={})

# ---------- (iii): solve recursion with polynomial ph_n, ph_0=1, compare ph_n(1) with h_n ----------
def solve(NN, free=lambda m: Fr(0)):
    phs=[{0:Fr(1)}]
    rs=[]
    for n in range(1,NN+1):
        cur=phs+[{}]+[{} for _ in range(NN+1-len(phs)-1)]
        F=[cur]
        for j in range(1,6): F.append(O(-j,F[-1]) if True else None)
        # recompute with truncation NN
        def g(Fm,m): return Fm[m] if m>=0 else {}
        rn={}
        for j in range(1,6):
            x=F[j-1][n]
            rn=ladd(rn, ladd(lscal(x,Fr(4*n+j,5)), theta_u(x), -1))
        rn=ladd(rn, F[4][n], -3)
        rn=ladd(rn, lmul({5:Fr(17,5)}, g(F[3],n-1)), -1)
        rn=ladd(rn, lmul({10:Fr(9,5)}, g(F[2],n-2)), -1)
        rn=ladd(rn, lmul({15:Fr(274,625)}, g(F[1],n-3)), -1)
        rn=ladd(rn, lmul({20:Fr(24,625)}, g(cur,n-4)), -1)
        rs.append(rn)
        new={}
        for e,v in rn.items():
            assert e%5==0 and e>=0
            k=e//5
            if 4*n-5*k==0:
                raise Exception("resonance violated n=%d: %s"%(n,v))
            new[e]=-v/(4*n-5*k)
        if n%5==0:
            fm=free(n//5)
            if fm!=0: new[4*n]=fm   # coefficient of u^(4n/5): exponent in s is 5*(4n/5)=4n
        phs.append(new)
    return phs,rs
global_N=N
N=6
phs,rs=solve(6)
print("r_1 =",rs[0]," ph_1 =",phs[1])
def at1(p): return sum(p.values())
z5=cmath.exp(2j*cmath.pi/5)
S3=sum((1-z5**j)**-3 for j in range(1,5)); S5=sum((1-z5**j)**-5 for j in range(1,5))
print("sum (1-zeta^j)^-3 =",S3)
a1=Fr(-3,20); a3=-(Fr(-1,30)/12)*(Fr(-1)-Fr(1,125))
h=[Fr(1),a1,a1**2/2,a1**3/6+a3,a1**4/24+a1*a3]
for n in range(1,5):
    print("n=%d  ph_n(u)=%s   ph_n(1)=%s   h_n=%s   equal:%s"%(n,{e//5:str(v) for e,v in sorted(phs[n].items())},at1(phs[n]),h[n],at1(phs[n])==h[n]))
print("deg check: resonance at n=5 passed (no exception); ph_5 has", len(phs[5]),"terms")
