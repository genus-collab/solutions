from fractions import Fraction as Fr
N=10
# series in u as dict exponent->coef (allow negative exponents), truncated at N
def D(f):
    # D = -(1-u) d/du
    out={}
    for e,c in f.items():
        if e==0: continue
        # d/du u^e = e u^{e-1}; times -(1-u) = -u^{e-1} + u^e
        out[e-1]=out.get(e-1,0)-e*c
        out[e]=out.get(e,0)+e*c
    return {e:c for e,c in out.items() if c!=0 and e<=N}
def add(*fs):
    out={}
    for f in fs:
        for e,c in f.items(): out[e]=out.get(e,0)+c
    return {e:c for e,c in out.items() if c!=0}
def scal(a,f): return {e:a*c for e,c in f.items()}
def mulu(f,k=1): return {e+k:c for e,c in f.items() if e+k<=N}
def PF4(f):
    D1=D(f);D2=D(D1);D3=D(D2);D4=D(D3)
    inner=add(scal(2,D3),scal(Fr(7,5),D2),scal(Fr(2,5),D1),scal(Fr(24,625),f))
    return add(mulu(D4),scal(-1,inner),mulu(inner))
# solve power series solution with given b0,b1,b2
def solve(b):
    f={i:Fr(b[i]) for i in range(3)}
    for n in range(3,N+1):
        f[n]=Fr(0)
        r=PF4(f).get(n-3,0)
        # coefficient of u^{n-3} has P_{-3}(n) b_n with P_{-3}(n)=n(n-1)^2(n-2)
        f[n]=-r/(n*(n-1)**2*(n-2))
    return f
f0=solve([1,0,0]);f1=solve([0,1,0]);f2=solve([0,0,1])
for f in (f0,f1,f2):
    r=PF4(f); print("PF4 residual low orders:",{e:c for e,c in r.items() if e<=N-4})
w1={e:f1.get(e,0)+Fr(7,10)*f2.get(e,0) for e in range(N+1)}
print("varpi1:",[w1[e] for e in range(5)])
print("f0:",[f0[e] for e in range(5)])
# tc = varpi1/f0 in frame (0,0)
inv={0:Fr(1)}
# compute 1/f0
g=[f0.get(e,0) for e in range(N+1)]
h=[Fr(0)]*(N+1);h[0]=Fr(1)/g[0]
for n in range(1,N+1):
    h[n]=-sum(g[k]*h[n-k] for k in range(1,n+1))/g[0]
tc=[sum(w1[k]*h[n-k] for k in range(n+1)) for n in range(N+1)]
print("tc frame(0,0):",tc[:5])
# P_{-1}'(0) = [u^{-1}] PF4(ell); compute PF4(ell) negative part: use D ell = 1-u^{-1}
# represent ell-derivatives: D ell = {0:1,-1:-1}
Dl={0:Fr(1),-1:Fr(-1)}
D2=D(Dl);D3=D(D2);D4=D(D3)
inner=add(scal(2,D3),scal(Fr(7,5),D2),scal(Fr(2,5),Dl))
r=add(mulu(D4),scal(-1,inner),mulu(inner))
print("r(1) negative part:",{e:c for e,c in sorted(r.items()) if e<0})
# obstruction matrix: [u^m] r(g) for g in f0,f1,f2 ; r(g) = PF4(g ell) - ell PF4 g ; compute via product rule numerically:
def r_of(g):
    # PF4(g*ell) non-ell part: D^k(g ell) = sum binom(k,j) D^{k-j}g D^j ell, j>=1 terms are non-ell part
    from math import comb
    Dg=[g]
    for k in range(4): Dg.append(D(Dg[-1]))
    Dle=[None,Dl,D2,D3,D4]
    def mul(a,b):
        out={}
        for e1,c1 in a.items():
            for e2,c2 in b.items():
                if e1+e2<=N: out[e1+e2]=out.get(e1+e2,0)+c1*c2
        return out
    def nonell(k):
        acc={}
        for j in range(1,k+1):
            acc=add(acc,scal(comb(k,j),mul(Dg[k-j],Dle[j])))
        return acc
    inner=add(scal(2,nonell(3)),scal(Fr(7,5),nonell(2)),scal(Fr(2,5),nonell(1)))
    return add(mulu(nonell(4)),scal(-1,inner),mulu(inner))
for name,g in (("f0",f0),("f1",f1),("f2",f2)):
    rr=r_of(g); print(name,[rr.get(m,0) for m in (-3,-2,-1)])
# tau, v
t=[Fr(4,5*n+4) for n in range(N+1)]
tau=[Fr(0)]*(N+1);tau[0]=Fr(1)
for n in range(1,N+1): tau[n]=-sum(t[k]*tau[n-k] for k in range(1,n+1))
print("tau:",tau[:4])
