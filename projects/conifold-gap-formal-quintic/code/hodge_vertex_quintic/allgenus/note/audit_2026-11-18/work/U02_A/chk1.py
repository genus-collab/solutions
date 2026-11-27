from fractions import Fraction as F
N=12
# series in u as lists of Fractions, index = power, allow negative via offset
def Dq(a):  # a: dict power->coef ; D = -(1-u) d/du
    r={}
    for n,c in a.items():
        if n!=0:
            r[n-1]=r.get(n-1,0)-n*c
            r[n]=r.get(n,0)+n*c
    return r
def mul(a,b):
    r={}
    for n,c in a.items():
        for m,d in b.items():
            if n+m<=N+2: r[n+m]=r.get(n+m,0)+c*d
    return r
def add(*xs):
    r={}
    for a in xs:
        for n,c in a.items(): r[n]=r.get(n,0)+c
    return r
def scal(k,a): return {n:k*c for n,c in a.items()}
u={1:F(1)}; one={0:F(1)}; omu={0:F(1),1:F(-1)}
def PF4(g):
    d1=Dq(g); d2=Dq(d1); d3=Dq(d2); d4=Dq(d3)
    inner=add(scal(2,d3),scal(F(7,5),d2),scal(F(2,5),d1),scal(F(24,625),g))
    return add(mul(u,d4),scal(-1,mul(omu,inner)))
# solve PF4 g=0 with given b0,b1,b2 via linear recursion: determine b_{n} successively
def solve(b0,b1,b2,rhs=None):
    g={0:F(b0),1:F(b1),2:F(b2)}
    for n in range(3,N+1):
        # coefficient of u^{n-3} in PF4 g is linear in b_n: compute with b_n=0 and b_n=1
        g0=dict(g); g0[n]=F(0)
        c0=PF4(g0).get(n-3,0)
        g1=dict(g); g1[n]=F(1)
        c1=PF4(g1).get(n-3,0)
        target = 0 if rhs is None else rhs.get(n-3,0)
        g[n]=(target-c0)/(c1-c0)
        assert c1-c0==(n)*(n-1)**2*(n-2), (n,c1-c0)
    return g
f0=solve(1,0,0); f1=solve(0,1,0); f2=solve(0,0,1)
for nm,f in (('f0',f0),('f1',f1),('f2',f2)):
    r=PF4(f); print(nm,[str(f.get(k,0)) for k in range(6)], 'residual low:',[str(r.get(k,0)) for k in range(-3,N-4)])
w1=add(f1,scal(F(7,10),f2))
print('varpi1',[str(w1.get(k,0)) for k in range(6)])
# tc in frame (0,0): w1/f0
def inv(a):
    r={0:1/a[0]}
    for n in range(1,N+1):
        s=sum(a.get(k,0)*r.get(n-k,0) for k in range(1,n+1))
        r[n]=-s/a[0]
    return r
tc=mul(w1,inv(f0))
print('tc',[str(tc.get(k,0)) for k in range(6)])
# P_j'(n) check: obstruction for log solution. PF4(g l) = l PF4 g + r(g). compute r via D acting on g*l: represent pair (A,B) = A*l + B
def Dpair(A,B):
    # D(A l + B) = (DA) l + A*(1-u^{-1}) + DB
    return Dq(A), add(Dq(B), A, scal(-1,{n-1:c for n,c in A.items()}))
def PF4pair(A,B):
    d=[(A,B)]
    for k in range(4): d.append(Dpair(*d[-1]))
    def comb(idx):
        return add(mul(u,d[4][idx]),scal(-1,mul(omu,add(scal(2,d[3][idx]),scal(F(7,5),d[2][idx]),scal(F(2,5),d[1][idx]),scal(F(24,625),d[0][idx])))))
    return comb(0),comb(1)
for nm,f in (('f0',f0),('f1',f1),('f2',f2),('w1',w1)):
    a,b=PF4pair(f,{})
    print(nm,'obstruction u^-3,u^-2,u^-1:',[str(b.get(k,0)) for k in (-3,-2,-1)])
# Jet matrix and Conc
def Conc(f,g):
    A={4:u,3:scal(-2,omu),2:scal(F(-7,5),omu),1:scal(F(-2,5),omu)}
    df=[f]; 
    for k in range(4): df.append(Dq(df[-1]))
    tot={}
    for k in range(1,5):
        w=mul(A[k],g)
        dw=[w]
        for j in range(4): dw.append(Dq(dw[-1]))
        for j in range(k):
            tot=add(tot,scal((-1)**j,mul(dw[j],df[k-1-j])))
    return tot
for (a,fa) in (('f0',f0),('f1',f1),('f2',f2)):
    for (b,fb) in (('f0',f0),('f1',f1),('f2',f2),('w1',w1)):
        c=Conc(fa,fb)
        print('Conc',a,b,[str(c.get(k,0)) for k in range(0,5)])
