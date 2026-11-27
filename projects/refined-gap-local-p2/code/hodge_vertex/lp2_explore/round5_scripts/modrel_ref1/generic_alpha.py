# Referee script: at GENERIC alpha over F_p, power series in y.
# (1) eta_ab = Sum Phi_a Phi_b / Delta_i with modrel-free Phi_2 = Y(Y-e1 E)/(A D), Delta^-1 = 1/(3 eps W).
# (2) WKB solution of eq:pf at generic alpha (no closed form), rows from eq:RHb with BS's D,E; test
#     B_ij = Sum R_a^i(z) eta^ab R_b^j(-z) =0 (i != j) and theta B_ii=0.  Controls perturb E, D.
import random,sys
p=2**61-1
inv=lambda x:pow(x%p,p-2,p)
N=int(sys.argv[1]) if len(sys.argv)>1 else 24
K=int(sys.argv[2]) if len(sys.argv)>2 else 5
def add(a,b): return [(x+y)%p for x,y in zip(a,b)]
def sub(a,b): return [(x-y)%p for x,y in zip(a,b)]
def scal(c,a): return [c*x%p for x in a]
def mul(a,b):
    r=[0]*N
    for i,x in enumerate(a):
        if x:
            for j in range(N-i): r[i+j]=(r[i+j]+x*b[j])%p
    return r
def sinv(a):
    r=[0]*N; i0=inv(a[0]); r[0]=i0
    for n in range(1,N):
        s=0
        for k in range(1,n+1): s=(s+a[k]*r[n-k])%p
        r[n]=(-s*i0)%p
    return r
def th(a): return [n*a[n]%p for n in range(N)]
def const(c): return [c%p]+[0]*(N-1)
yy=[0,1]+[0]*(N-2)
# series from a_d,H_n
a=[1]; 
for d in range(1,N): a.append(a[-1]*3*(3*d-1)*(3*d-2)%p*inv(d*d)%p)
H=[0]
for n in range(1,3*N+1): H.append((H[-1]+inv(n))%p)
sg=lambda d:1 if d%2==0 else p-1
A=[a[d]*sg(d)%p for d in range(N)]
G=[(-a[d]*H[d]*sg(d))%p for d in range(N)]
h=[3*a[d]*(H[3*d]-H[d])%p*sg(d)%p for d in range(N)]
Ai=sinv(A)
Ey=th(mul(G,Ai))
D=add(const(1),th(mul(h,Ai)))
Dc=add(const(1),scal(27,yy)); X=sinv(Dc)
def roots(al):
    e1=sum(al)%p; e2=(al[0]*al[1]+al[0]*al[2]+al[1]*al[2])%p; e3=al[0]*al[1]*al[2]%p
    out=[]
    for i in range(3):
        Y=const(-al[i])
        for it in range(8):
            Y2=mul(Y,Y);Y3=mul(Y2,Y)
            L=add(add(add(mul(Dc,Y3),scal(e1,Y2)),scal(e2,Y)),const(e3))
            Lp=add(add(scal(3,mul(Dc,Y2)),scal(2*e1,Y)),const(e2))
            Y=sub(Y,mul(L,sinv(Lp)))
        out.append(Y)
    return e1,e2,e3,out
def test(al,e_1,e_2,ctrl=None):
    eps=e_1*e_2%p
    e1,e2,e3,Ys=roots(al)
    Ee=Ey; DD=D
    if ctrl=='E': Ee=add(Ey,scal(5,[0,0,0,1]+[0]*(N-4)))
    if ctrl=='E2': Ee=sub(D,X)         # missing the 1/2
    if ctrl=='D': DD=add(D,scal(5,[0,0,0,1]+[0]*(N-4)))
    Ws=[add(add(scal(e1,mul(Y,Y)),scal(2*e2,Y)),const(3*e3)) for Y in Ys]
    # eta_ab
    res={}
    i3=inv(3*eps)
    target={(0,0):const(i3*inv(e3)),(0,1):const(0),(0,2):const(0),(1,1):const(0),(1,2):const(-i3),(2,2):const(i3*e1)}
    Phi=[[const(1) for _ in Ys],[mul(Y,Ai) for Y in Ys],[mul(mul(Y,sub(Y,scal(e1,Ee))),sinv(mul(A,DD))) for Y in Ys]]
    fails=[]
    for ab,tg in target.items():
        s=[0]*N
        for i in range(3):
            s=add(s,mul(mul(Phi[ab[0]][i],Phi[ab[1]][i]),scal(i3,sinv(Ws[i]))))
        if s!=tg: fails.append(ab)
    # WKB
    Rows=[]
    for i,Y in enumerate(Ys):
        W=Ws[i]
        # f0=(W0/W)^{1/2}: theta g=-(1/2)(theta W/W) g
        q=mul(th(W),sinv(W))
        g=[0]*N; g[0]=1
        for n in range(1,N):
            s=0
            for k in range(1,n+1): s=(s+q[k]*g[n-k])%p
            g[n]=(-s*inv(2)%p)*inv(n)%p
        f=[g]
        def Th(fa,KK):
            out=[]
            for k in range(KK+1):
                x=mul(Y,fa[k]) if k<len(fa) else [0]*N
                if k>=1 and k-1<len(fa): x=add(x,th(fa[k-1]))
                out.append(x)
            return out
        def Pop(fa,KK):
            # prod (Th+alpha_j) f + y (3Th)(3Th+z)(3Th+2z) f
            cur=[fa[k] if k<len(fa) else [0]*N for k in range(KK+1)]
            res=cur
            for j in range(3):
                t=Th(res,KK)
                res=[add(t[k],scal(al[j],res[k])) for k in range(KK+1)]
            # second
            r=[cur[k] for k in range(KK+1)]
            for jj in [2,1,0]:
                t=Th(r,KK)
                new=[]
                for k in range(KK+1):
                    x=scal(3,t[k])
                    if k>=1 and jj>0: x=add(x,scal(jj,r[k-1]))
                    new.append(x)
                r=new
            return [add(res[k],mul(yy,r[k])) for k in range(KK+1)]
        for k in range(1,K+1):
            Gk=Pop(f,k+1)[k+1]
            # check lower vanish
            rhs=mul(mul(mul(Y,sinv(W)),Gk),sinv(f[0]))
            assert rhs[0]==0,"log term at y^0"
            hk=[0]*N
            for n in range(1,N): hk[n]=rhs[n]*inv(n)%p
            f.append(mul(f[0],hk))
        chk=Pop(f,K+1)
        for k in range(K+1): assert not any(chk[k]),("PF fails",k)
        R0=f
        T=Th(R0,K); R1=[mul(x,Ai) for x in T]
        T1=Th(R1,K)
        R2=[mul(sub(T1[k],mul(scal(e1,Ee),R1[k])),sinv(DD)) for k in range(K+1)]
        Rows.append((R0,R1,R2))
    # inverse metric
    et00=3*eps*e3%p; et11=(-3*eps*e1)%p; et12=(-3*eps)%p
    Bf=[]
    def sgn(fa): return [ (x if k%2==0 else scal(p-1,x)) for k,x in enumerate(fa)]
    def conv(fa,gb):
        out=[[0]*N for _ in range(K+1)]
        for x in range(K+1):
            for y2 in range(K+1-x):
                out[x+y2]=add(out[x+y2],mul(fa[x],gb[y2]))
        return out
    vals={}
    for i in range(3):
        for j in range(i,3):
            Pi=Rows[i];Pj=Rows[j]
            A0=sgn(Pj[0]);A1=sgn(Pj[1]);A2=sgn(Pj[2])
            t00=conv(Pi[0],A0);t11=conv(Pi[1],A1);t12=conv(Pi[1],A2);t21=conv(Pi[2],A1)
            for k in range(K+1):
                v=add(add(scal(et00,t00[k]),scal(et11,t11[k])),scal(et12,add(t12[k],t21[k])))
                if i!=j:
                    if any(v): Bf.append((i,j,k))
                else:
                    if any(v[1:]): Bf.append((i,j,k,'nonconst'))
                    if k==0: vals[i]=v[0]
    w0=[Ws[i][0] for i in range(3)]
    z0ok=all(vals[i]==3*eps*w0[i]%p for i in range(3))
    return fails,Bf,z0ok
rnd=random.Random(3)
for trial in range(3):
    al=[rnd.randrange(2,p) for _ in range(3)]
    e_1,e_2=rnd.randrange(2,p),rnd.randrange(2,p)
    print('generic',test(al,e_1,e_2),flush=True)
al=[11,23,inv(7)]
print('small',test(al,5,2))
for c in ['E','E2','D']:
    f,b,z=test(al,5,2,ctrl=c)
    print('ctrl',c,'eta fails',f,'B fails',len(b),b[:4],'z0 ok',z)
