# Generic alpha: Givental graph sum with rows from the I-function (J=I restricted to P_j), times h_{wt_j} (BS eq:canR with c_j),
# as exact truncated y-series, against direct localisation. Tests Stage A (A2, A3, graph sum), the sign of c_j, and the role of the eps weights.
# Evidence only.
import sys, math, itertools, time
from fractions import Fraction as Fr
from collections import Counter
from hod import tau, B
from tf import theta_coeffs, mul_ser, stable
from loc import Loc

class Ctx:
    def __init__(s,ny): s.ny=ny
def ymul(a,b,ny):
    o=[Fr(0)]*ny
    for i,x in enumerate(a):
        if x==0: continue
        for j in range(ny-i):
            o[i+j]+=x*b[j]
    return o
def yinv(a,ny):
    o=[Fr(0)]*ny; o[0]=1/a[0]
    for n in range(1,ny):
        o[n]=-sum(a[k]*o[n-k] for k in range(1,n+1))/a[0]
    return o
def ysqrt(a,ny):
    # a[0]=1
    o=[Fr(0)]*ny; o[0]=Fr(1)
    for n in range(1,ny):
        o[n]=(a[n]-sum(o[k]*o[n-k] for k in range(1,n)))/2
    return o
def yadd(a,b): return [x+y for x,y in zip(a,b)]
def ytheta(a): return [k*x for k,x in enumerate(a)]
def yscal(a,c): return [x*c for x in a]

class Gen:
    def __init__(s,alpha,e1,e2,ny,Z,sign=1,eps_in_h=True,mode_h='wt'):
        s.a=[Fr(x) for x in alpha]; s.e1=Fr(e1); s.e2=Fr(e2); s.ny=ny; s.Z=Z
        ny_=ny
        a=s.a; kap=s.e1*s.e2
        e1a=sum(a); e2a=a[0]*a[1]+a[0]*a[2]+a[1]*a[2]; e3a=a[0]*a[1]*a[2]
        s.e1a=e1a
        # A,D,E(y)
        ad=[Fr(1)]+[Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d for d in range(1,ny)]
        Dc=[Fr(1),Fr(27)]+[Fr(0)]*(ny-2) if ny>=2 else [Fr(1)]
        X=yinv(Dc,ny)
        Ai=yinv(ad,ny)
        G=[Fr(0)]+[-Fr(math.factorial(3*d),math.factorial(d)**3)*(-1)**d*sum(Fr(1,k) for k in range(1,d+1)) for d in range(1,ny)]
        Ey=ytheta(ymul(G,Ai,ny))           # E(y) by definition
        Dy=ymul(X,ymul(Ai,Ai,ny),ny)       # D=X/A^2
        s.chk_modrel=[2*Ey[k]-(Dy[k]-X[k]) for k in range(ny)]
        s.A=ad; s.Ai=Ai; s.Di=yinv(Dy,ny); s.Ey=Ey
        s.Q=None
        # eta^{ab}
        W=[]
        for i in range(3):
            j,k=[x for x in range(3) if x!=i]
            W.append(3*a[i]*(a[j]-a[i])*(a[k]-a[i])*kap)
        s.wj=W
        Hv=[-a[i] for i in range(3)]
        eta=[[sum(Hv[k]**(p+q)/W[k] for k in range(3)) for q in range(3)] for p in range(3)]
        A_=eta
        det=(A_[0][0]*(A_[1][1]*A_[2][2]-A_[1][2]*A_[2][1])-A_[0][1]*(A_[1][0]*A_[2][2]-A_[1][2]*A_[2][0])+A_[0][2]*(A_[1][0]*A_[2][1]-A_[1][1]*A_[2][0]))
        adj=[[0]*3 for _ in range(3)]
        for p in range(3):
            for q in range(3):
                r=[x for x in range(3) if x!=p]; c=[x for x in range(3) if x!=q]
                adj[q][p]=(-1)**(p+q)*(A_[r[0]][c[0]]*A_[r[1]][c[1]]-A_[r[0]][c[1]]*A_[r[1]][c[0]])/det
        s.etai=adj
        s.rows=[]; s.Del=[]
        for j in range(3):
            aj=[a[l]-a[j] for l in range(3) if l!=j]
            # Y_j(y)
            Y=[Fr(0)]*ny; Y[0]=-a[j]
            # fixed point iteration via Newton (coefficient recursion): delta = -27 y Y^3 / prod_{l!=j}(Y+a_l)
            for it in range(ny):
                Y3=ymul(Y,ymul(Y,Y,ny),ny)
                den=[Fr(1)]+[Fr(0)]*(ny-1)
                for l in range(3):
                    if l!=j:
                        q=list(Y); q[0]+=a[l]
                        den=ymul(den,q,ny)
                rhs=ymul(Y3,yinv(den,ny),ny)
                new=[Fr(0)]*ny; new[0]=-a[j]
                for n in range(1,ny):
                    new[n]=-27*rhs[n-1]
                Y=new
            s_Y=Y
            # Delta_j=3 kap W(Y)
            Y2=ymul(Y,Y,ny)
            WY=yadd(yadd(yscal(Y2,e1a),yscal(Y,2*e2a)),[3*e3a]+[Fr(0)]*(ny-1))
            Del=yscal(WY,3*kap)
            s.Del.append(Del)
            ratio=yscal(Del,1/W[j])
            sq=ysqrt(ratio,ny)       # (Delta/w)^{1/2}
            # delta=Y+a_j, uhat=Sum delta_k y^k/k
            uh=[Fr(0)]*ny
            for k in range(1,ny):
                uh[k]=Y[k]/k
            # z-Laurent series: dict? represent a y-coefficient as list in z from -ny..Z+ny
            low=-(ny); top=Z+ny+2
            LEN=top-low+1
            def zero(): return [Fr(0)]*LEN
            def zmul(u,v):
                o=zero()
                for i,x in enumerate(u):
                    if x==0: continue
                    for k,w in enumerate(v):
                        if w==0: continue
                        e=i+k+2*low
                        if low<=e<=top: o[e-low]+=x*w
                return o
            # C_d(z)
            def serdiv_lin(ser,c0,c1):  # ser/(c0+c1 z) with z>=0 series as list by power
                o=[Fr(0)]*len(ser)
                for n in range(len(ser)):
                    o[n]=(ser[n]-(c1*o[n-1] if n>0 else 0))/c0
                return o
            I=[zero() for _ in range(ny)]
            for d in range(ny):
                ser=[Fr(0)]*(top+d+1); ser[0]=Fr(1)
                # numerator prod_{m=0}^{3d-1}(3a_j-m z)
                for m in range(3*d):
                    o=[Fr(0)]*len(ser)
                    for n in range(len(ser)):
                        o[n]=3*a[j]*ser[n]-(m*ser[n-1] if n>0 else 0)
                    ser=o
                for al in aj:
                    for m in range(1,d+1):
                        ser=serdiv_lin(ser,al,Fr(m))
                # 1/(d! z^d)
                f=Fr(1,math.factorial(d))
                for n in range(len(ser)):
                    e=n-d
                    if low<=e<=top: I[d][e-low]=ser[n]*f
            # exp(-uh/z)
            ex=[zero() for _ in range(ny)]
            ex[0][0-low]=Fr(1)
            # exp series by recursion in y: n ex_n = Sum_k k (-uh_k/z) ex_{n-k}
            for n in range(1,ny):
                t=zero()
                for k in range(1,n+1):
                    c=-k*uh[k]
                    if c==0: continue
                    src=ex[n-k]
                    for i,x in enumerate(src):
                        if x!=0 and i-1>=0: t[i-1]+=c*x
                ex[n]=[x/n for x in t]
            # product
            R=[zero() for _ in range(ny)]
            for p in range(ny):
                for q in range(ny-p):
                    w=zmul(ex[p],I[q])
                    R[p+q]=[x+y for x,y in zip(R[p+q],w)]
            # poles in z must cancel
            pol=[R[n][i] for n in range(ny) for i in range(0,-low)]
            s.polecheck=all(x==0 for x in pol)
            # now coefficient arrays: Zser[k] = y-series of [z^k] R
            Zs=[[R[n][k-low] for n in range(ny)] for k in range(Z+3)]
            # multiply by sq : P0 = sq * R * h
            Zs=[ymul(c,sq,ny) for c in Zs]
            # h_{wt}
            wt=[3*a[j]]+aj
            if eps_in_h: wt=wt+[-s.e1,-s.e2]
            if mode_h=='wt':
                hs=[Fr(1)]+[Fr(0)]*(Z+2)
                for w in wt:
                    hs=mul_ser(hs,theta_coeffs(-Fr(sign)/w,Z+2),Z+2)
            else:
                hs=[Fr(1)]+[Fr(0)]*(Z+2)
            P0=[[Fr(0)]*ny for _ in range(Z+3)]
            for k in range(Z+3):
                for l in range(k+1):
                    if hs[l]!=0:
                        P0[k]=yadd(P0[k],yscal(Zs[k-l],hs[l]))
            # T = z theta + Y - (z/2) theta log Del
            dl=ymul(ytheta(Del),yinv(Del,ny),ny)
            def Top(ser,extra=None):
                o=[]
                for k in range(Z+3):
                    x=ymul(Y,ser[k],ny)
                    if extra is not None: x=yadd(x,yscal(ymul(extra,ser[k],ny),-1))
                    if k>=1:
                        x=yadd(x,ytheta(ser[k-1]))
                        x=yadd(x,yscal(ymul(dl,ser[k-1],ny),Fr(-1,2)))
                    o.append(x)
                return o
            P1=[ymul(x,Ai,ny) for x in Top(P0)]
            P2=[ymul(x,s.Di,ny) for x in Top(P1,extra=yscal(Ey,e1a))]
            s.rows.append((P0,P1,P2))
        s.em={}; s.vm={}; s.pm={}
    def pair(s,i,j,M,b):
        key=(i,j,M,b)
        if key in s.pm: return s.pm[key]
        o=[Fr(0)]*s.ny
        for p in range(3):
            for q in range(3):
                c=s.etai[p][q]
                if c==0: continue
                o=yadd(o,yscal(ymul(s.rows[i][p][M],s.rows[j][q][b],s.ny),c))
        s.pm[key]=o
        return o
    def sympl_check(s,kmax=6):
        bad=0
        for i in range(3):
            for j in range(3):
                for k in range(kmax+1):
                    t=[Fr(0)]*s.ny
                    for m in range(k+1):
                        t=yadd(t,yscal(s.pair(i,j,m,k-m),Fr((-1)**(k-m))))
                    if k==0 and i==j:
                        t=[x-y for x,y in zip(t,s.Del[i])]
                    if any(x!=0 for x in t): bad+=1
        return bad
    def edge(s,i,j,m,n):
        key=(i,j,m,n)
        if key in s.em: return s.em[key]
        o=[Fr(0)]*s.ny
        for k in range(n+1):
            M=m+k+1; b=n-k
            o=yadd(o,yscal(s.pair(i,j,M,b),Fr(-1)*(-1)**(M+b)*(-1)**k))
        s.em[key]=o
        return o
    def vert(s,i,g,hs):
        hs=tuple(sorted(hs)); key=(i,g,hs)
        if key in s.vm: return s.vm[key]
        n=len(hs); B_=3*g-3+n-sum(hs); tot=[Fr(0)]*s.ny
        if B_>=0:
            def ms(rem,mx):
                if rem==0: yield (); return
                for p in range(min(rem,mx),0,-1):
                    for r in ms(rem-p,p): yield (p+1,)+r
            for ds in ms(B_,B_):
                val=tau(g,tuple(sorted(hs+ds)))
                if val==0: continue
                sym=1
                for c in Counter(ds).values(): sym*=math.factorial(c)
                term=[Fr(0)]*s.ny; term[0]=val/sym
                for d in ds:
                    term=ymul(term,yscal(s.rows[i][0][d-1],Fr((-1)**d)),s.ny)
                tot=yadd(tot,term)
            if g>=2:
                for _ in range(g-1): tot=ymul(tot,s.Del[i],s.ny)
            elif g==0:
                tot=ymul(tot,yinv(s.Del[i],s.ny),s.ny)
        s.vm[key]=tot
        return tot
    def total(s,g):
        res=[Fr(0)]*s.ny
        for gen,ed,aut in stable(g):
            nv=len(gen); val=[0]*nv
            for u,v in ed: val[u]+=1; val[v]+=1
            bud=[3*gen[i]-3+val[i] for i in range(nv)]
            hs=[[] for _ in range(nv)]
            acc=[Fr(0)]*s.ny
            for col in itertools.product(range(3),repeat=nv):
                def go(k,pref):
                    nonlocal acc
                    if k==len(ed):
                        x=pref
                        for i in range(nv):
                            x=ymul(x,s.vert(col[i],gen[i],tuple(hs[i])),s.ny)
                            if all(c==0 for c in x): return
                        acc=yadd(acc,x); return
                    u,v=ed[k]
                    ru=bud[u]-sum(hs[u])
                    for a in range(ru+1):
                        rv=(bud[v]-sum(hs[v])) if u!=v else ru-a
                        for b in range(rv+1):
                            e=s.edge(col[u],col[v],a,b)
                            if all(c==0 for c in e): continue
                            hs[u].append(a); hs[v].append(b)
                            go(k+1,ymul(pref,e,s.ny))
                            hs[u].pop(); hs[v].pop()
                go(0,[Fr(1)]+[Fr(0)]*(s.ny-1))
            res=yadd(res,yscal(acc,Fr(1,aut)))
        return res

def Qseries(ny):
    g=[Fr(0)]+[3*Fr(math.factorial(3*d-1),math.factorial(d)**3)*(-1)**d for d in range(1,ny)]
    ex=[Fr(0)]*ny; ex[0]=Fr(1)
    for n in range(1,ny):
        ex[n]=sum(k*g[k]*ex[n-k] for k in range(1,n+1))/n
    return ex   # Q/y

if __name__=="__main__":
    g=int(sys.argv[1]); ny=int(sys.argv[2])
    Z=3*g
    pts=[((Fr(2,3),Fr(5,7),Fr(-3,11)),3,7),((Fr(1,2),Fr(2),Fr(-5,3)),Fr(1,3),Fr(4,5))]
    if len(sys.argv)>3: pts=pts[:int(sys.argv[3])]
    for al,e1,e2 in pts:
        t0=time.time()
        G=Gen(al,e1,e2,ny,Z)
        print(f"g={g} alpha={[str(x) for x in al]} eps=({e1},{e2}) z-poles cancel:",G.polecheck," modrel check:",all(x==0 for x in G.chk_modrel)," symplecticity violations:",G.sympl_check(),flush=True)
        F=G.total(g)
        # localisation
        L=Loc(al,e1,e2)
        gw=[L.gw(g,d) for d in range(ny)]
        Qy=Qseries(ny)
        # Sum_d Q^d GW = Sum_d y^d (Q/y)^d GW
        ser=[Fr(0)]*ny; pw=[Fr(1)]+[Fr(0)]*(ny-1)
        for d in range(ny):
            for k in range(ny-d):
                ser[d+k]+=gw[d]*pw[k]
            pw=ymul(pw,Qy,ny)
        print("   graph sum - localisation at y^k:",[str(a-b) for a,b in zip(F,ser)],f"[{time.time()-t0:.0f}s]",flush=True)
        print("   localisation y^k:",[float(x) for x in ser],flush=True)
        for lbl,kw in (("sign flip of c_j",dict(sign=-1)),("no eps weights in h",dict(eps_in_h=False)),("no h at all",dict(mode_h='none'))):
            Gb=Gen(al,e1,e2,ny,Z,**kw)
            Fb=Gb.total(g)
            print("   control",lbl,": diff:",[float(a-b) for a,b in zip(Fb,ser)],flush=True)
