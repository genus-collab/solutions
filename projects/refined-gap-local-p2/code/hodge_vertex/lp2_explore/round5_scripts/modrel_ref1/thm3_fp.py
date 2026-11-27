# Referee script: Prop 2.10 / Theorem 3 literally from Definition 1.3, with a RANDOM section c1, over F_p (p=2^61-1, p=1 mod 3),
# Taylor series in s around t=t0+s. Also Prop 3.1 for random c1. Negative controls.
import random,sys
p=2**61-1
inv=lambda x:pow(x%p,p-2,p)
# primitive cube root
def cuberoot():
    g=3
    while True:
        w=pow(g,(p-1)//3,p)
        if w!=1: return w
        g+=1
omega=cuberoot()
assert (omega*omega+omega+1)%p==0
def recur(n):
    r=[1]
    rm=0
    for m in range(0,n):
        rp=(3*(3*m*m+3*m+1)*r[m]-m**3*(r[m-1] if m>=1 else 0))%p*inv(27*(m+1))%p
        r.append(rp)
    return r
class Ser:
    N=30
def add(a,b): return [(x+y)%p for x,y in zip(a,b)]
def sub(a,b): return [(x-y)%p for x,y in zip(a,b)]
def scal(c,a): return [c*x%p for x in a]
def mul(a,b):
    N=len(a); r=[0]*N
    for i,x in enumerate(a):
        if x:
            for j in range(N-i):
                r[i+j]=(r[i+j]+x*b[j])%p
    return r
def sinv(a):
    N=len(a); r=[0]*N; i0=inv(a[0]); r[0]=i0
    for n in range(1,N):
        s=0
        for k in range(1,n+1): s=(s+a[k]*r[n-k])%p
        r[n]=(-s*i0)%p
    return r
def run(e1,e2,zmax,Ns,t0,seed,ctrl=None):
    rnd=random.Random(seed)
    E=(e1+e2)%p; eps=e1*e2%p
    def ser_const(c): return [c%p]+[0]*(Ns-1)
    tser=[t0%p,1]+[0]*(Ns-2)
    def d(a): return [ (k+1)*a[k+1]%p for k in range(Ns-1)]+[0]
    third=inv(3)
    def th(a): return scal(third,mul(tser,d(a)))  # theta=(t/3)d/dt
    y=mul(tser,mul(tser,tser))
    Dc=add(ser_const(1),scal(27,y))
    Dci=sinv(Dc)
    # random section c1
    c1=[rnd.randrange(1,p) for _ in range(Ns)]
    C1i=sinv(c1)
    XX=mul(th(c1),C1i)
    d2=mul(scal(E*inv(2)%p,sub(mul(C1i,C1i),ser_const(1))),Dci)
    if ctrl=='d2': d2=scal(2,d2)
    r=recur(zmax+2)
    if ctrl=='r3': r[3]=(r[3]+1)%p
    # coefficient [z^n] of prod_{j=1}^{L}(1+jz/E)^-1
    Ei=inv(E)
    def prodinv(L):
        c=[1]+[0]*zmax
        for j in range(1,L+1):
            f=j*Ei%p
            # divide by (1+f z)
            for n in range(1,zmax+1): c[n]=(c[n]-f*c[n-1])%p
        return c
    PI=[prodinv(l+1) for l in range(zmax+1)]
    rows={}
    for i in range(3):
        wi=pow(omega,i,p)
        ti=scal(wi,tser)
        den=sinv(scal(E,ti))
        v=mul(add(ser_const(1),scal(3,ti)),den)
        Y=sub(sinv(scal(E,v)),ser_const(0)) # placeholder, 1/(v) not E v
        Y=sub(sinv(v),ser_const(E*third))    # Y=1/v - E/3
        Ydirect=scal(-(E*third)%p,sinv(add(ser_const(1),scal(3,ti))))
        assert Y==Ydirect, "Y mismatch"
        # P0 = v S, z array
        vp=[v]
        for l in range(zmax+1): vp.append(mul(vp[-1],v))
        P0=[[0]*Ns for _ in range(zmax+1)]
        for l in range(zmax+1):
            for k in range(l,zmax+1):
                co=r[l]*PI[l][k-l]%p
                P0[k]=add(P0[k],scal(co,vp[l]))  # v^{l+1}
        def Th(f):
            g=[]
            for k in range(len(f)):
                x=mul(Y,f[k])
                if k>=1: x=add(x,th(f[k-1]))
                g.append(x)
            return g
        P1=[mul(C1i,x) for x in Th(P0)]
        ThP1=Th(P1)
        pref=mul(Dc,mul(c1,c1))
        P2=[mul(pref,sub(ThP1[k],mul(d2,P1[k]))) for k in range(zmax+1)]
        rows[i]=(P0,P1,P2,Y,Th)
    et00=eps*E%p*E%p*E%p*inv(9)%p
    if ctrl=='eta00': et00=2*et00%p
    et11=(-3*eps*E)%p; et12=(-3*eps)%p
    c0=3*eps*E%p
    margin=Ns-zmax-4
    bad=[]
    def sgn(f): return [ (x if k%2==0 else scal(p-1,x)) for k,x in enumerate(f)]
    def conv(f,g):
        out=[[0]*Ns for _ in range(zmax+1)]
        for a in range(zmax+1):
            for b in range(zmax+1-a):
                out[a+b]=add(out[a+b],mul(f[a],g[b]))
        return out
    for i in range(3):
        for j in range(i,3):
            Pi=rows[i];Pj=rows[j]
            A0=sgn(Pj[0]);A1=sgn(Pj[1]);A2=sgn(Pj[2])
            t00=conv(Pi[0],A0);t11=conv(Pi[1],A1);t12=conv(Pi[1],A2);t21=conv(Pi[2],A1)
            for k in range(zmax+1):
                val=add(add(scal(et00,t00[k]),scal(et11,t11[k])),scal(et12,add(t12[k],t21[k])))
                exp=[0]*Ns
                if i==j and k==0: exp[0]=c0
                dif=[(val[n]-exp[n])%p for n in range(margin)]
                if any(dif): bad.append((i,j,k))
    # Prop 3.1 at colour 0..2, general c1, up to zmax-1
    bad31=[]
    Rd=add(add(mul(Dc,th(XX)),mul(Dc,mul(XX,XX))),add(scal(27,mul(y,XX)),scal(6,y)))
    a22=mul(scal(p-1,C1i),add(scal(E,Dci),d2))
    a21=mul(scal(p-1,c1),add(add(scal(E,d2),mul(Dc,mul(d2,d2))),ser_const(E*E%p*third%p)))
    C3=mul(C1i,mul(C1i,C1i))
    for i in range(3):
        P0,P1,P2,Y,Th=rows[i]
        for lab,lhsrow,rhsf in [
          (0,P0,lambda k:P1[k]),
          (1,P1,lambda k:add(mul(mul(C1i,d2),P1[k]),mul(mul(C3,Dci),P2[k]))),
          (2,P2,lambda k:add(add(scal((-E*E*E%p)*inv(27)%p,P0[k]),mul(a21,P1[k])),add(mul(a22,P2[k]), scal(1,[0]*Ns) if k<2 else scal(p-1,mul(c1,mul(Rd,P1[k-2]))))))]:
            L=Th(lhsrow)
            for k in range(zmax):
                lhs=mul(C1i,L[k]); rhs=rhsf(k)
                if any((lhs[n]-rhs[n])%p for n in range(margin)): bad31.append((i,lab,k))
    return bad,bad31
if __name__=="__main__":
    rnd=random.Random(20261102)
    weights=[(1,1),(7,p-3),(inv(2),3*inv(5)%p),(5,2),(2,5),(1,p-2)]+[(rnd.randrange(2,p),rnd.randrange(2,p)) for _ in range(3)]
    zmax=int(sys.argv[1]) if len(sys.argv)>1 else 14
    for (e1,e2) in weights:
        Ns=zmax+8
        t0=rnd.randrange(2,p)
        b,b31=run(e1,e2,zmax,Ns,t0,rnd.randrange(10**9))
        print("weights",e1 if e1<1000 else 'big',e2 if e2<1000 else 'big','Pi fails',len(b),'Prop3.1 fails',len(b31),flush=True)
    print("controls")
    for ctrl in ['d2','r3','eta00']:
        b,b31=run(5,2,8,16,12345,7,ctrl)
        print(ctrl,'Pi fails',b[:6],'count',len(b),'P3.1 fails',len(b31))
