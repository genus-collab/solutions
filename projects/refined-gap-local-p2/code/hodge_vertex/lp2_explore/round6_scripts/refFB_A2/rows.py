# Rows at generic alpha, derived from the I-function alone (independent of BS eq:RHb coefficients):
# C^c_{11} determined from z->infinity behaviour of S(H)|_{P_j} at three colours.
from fractions import Fraction as Fr
from math import factorial
from engine import bern

class Ctx: pass

# y-series: list length D+1 of Fractions
def ymul(a,b,D):
    r=[Fr(0)]*(D+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,w in enumerate(b):
            if i+j>D: break
            r[i+j]+=x*w
    return r
def yinv(a,D):
    r=[Fr(0)]*(D+1); r[0]=1/a[0]
    for n in range(1,D+1):
        r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1))/a[0]
    return r
def ysqrt1(a,D):
    # a[0]=1
    assert a[0]==1
    r=[Fr(0)]*(D+1); r[0]=Fr(1)
    for n in range(1,D+1):
        r[n]=(a[n]-sum(r[k]*r[n-k] for k in range(1,n)))/2
    return r
def ytheta(a): return [i*x for i,x in enumerate(a)]
def yadd(a,b): return [x+w for x,w in zip(a,b)]
def yscal(c,a): return [c*x for x in a]
def ypow(a,n,D):
    if n<0: return ypow(yinv(a,D),-n,D)
    r=[Fr(0)]*(D+1); r[0]=Fr(1)
    for _ in range(n): r=ymul(r,a,D)
    return r

# yz-series: list length D+1 of dict z-power -> Fraction
def yzmul(a,b,D,Z):
    r=[{} for _ in range(D+1)]
    for i in range(D+1):
        for j in range(D+1-i):
            for p,x in a[i].items():
                for q,w in b[j].items():
                    if p+q<=Z:
                        r[i+j][p+q]=r[i+j].get(p+q,Fr(0))+x*w
    return r
def yz_from_y(a): return [{0:x} if x!=0 else {} for x in a]
def yz_zscale(a,s): return [{p*1:x*(s**p) for p,x in c.items()} for c in a]

def yzop(r,Y,D,Z,extra=None):
    """(z theta + Y + extra) r"""
    out=[{} for _ in range(D+1)]
    for d in range(D+1):
        for p,x in r[d].items():
            if d!=0 and p+1<=Z:
                out[d][p+1]=out[d].get(p+1,Fr(0))+d*x
    Yz=yz_from_y(Y if extra is None else yadd(Y,extra))
    m=yzmul(Yz,r,D,Z)
    for d in range(D+1):
        for p,x in m[d].items():
            out[d][p]=out[d].get(p,Fr(0))+x
    return out

def theta_series(z_power_exps,Z):
    pass

def h_series(wts,Z,sign=1):
    """h_v(z)^{sign}=exp(-sign sum B_{2k}/(2k(2k-1)) p_{2k-1}(1/v) z^{2k-1}), list by power of z to Z"""
    e=[Fr(0)]*(Z+1)
    k=1
    while 2*k-1<=Z:
        e[2*k-1]=-sign*bern(2*k)/(2*k*(2*k-1))*sum(Fr(1)/w**(2*k-1) for w in wts)
        k+=1
    # exp
    r=[Fr(0)]*(Z+1); r[0]=Fr(1)
    for n in range(1,Z+1):
        r[n]=sum(k*e[k]*r[n-k] for k in range(1,n+1))/n
    return r

def setup(alpha,eps,D,Zmax,cflip=False):
    c=Ctx(); c.alpha=alpha; c.eps=eps; c.D=D; c.Zmax=Zmax
    Z=Zmax+D+3; c.Z=Z
    a0,a1,a2=alpha
    e1=a0+a1+a2; e2=a0*a1+a0*a2+a1*a2; e3=a0*a1*a2
    c.e1,c.e2,c.e3=e1,e2,e3
    # mirror map
    logQ=[Fr(0)]*(D+1)
    for d in range(1,D+1):
        logQ[d]=3*Fr(factorial(3*d-1),factorial(d)**3)*(-1)**d
    c.logQhat=logQ # log Q=log y+logQhat
    A=yadd([Fr(1)]+[Fr(0)]*D,ytheta(logQ)); c.A=A; Ainv=yinv(A,D); c.Ainv=Ainv
    # Q(y)
    Q=[Fr(0)]*(D+1); ex=[Fr(0)]*(D+1); ex[0]=Fr(1)
    for n in range(1,D+1):
        ex[n]=sum(k*logQ[k]*ex[n-k] for k in range(1,n+1))/n
    c.Qhat=ex  # Q=y*Qhat
    # large z: sigma1,sigma2 for each colour; determine C
    sig=[]
    for j in range(3):
        aj=alpha[j]
        s1=[Fr(0)]*(D+1); s2=[Fr(0)]*(D+1)
        for d in range(1,D+1):
            # shat_d = x*3a*prod_{m=1}^{3d-1}(3a x-m)/(d! prod_{l!=j} prod_{m=1}^d((al-a)x+m)), expand to x^2
            num=[3*aj]  # poly in x, trunc deg 1 after factoring x
            num=[3*aj,Fr(0)]
            for m in range(1,3*d):
                num=[num[0]*(-m), num[1]*(-m)+num[0]*3*aj]
            den=[Fr(factorial(d)),Fr(0)]
            for l in range(3):
                if l==j: continue
                for m in range(1,d+1):
                    den=[den[0]*m, den[1]*m+den[0]*(alpha[l]-aj)]
            # ratio
            r0=num[0]/den[0]; r1=(num[1]-r0*den[1])/den[0]
            s1[d]=r0; s2[d]=r1
        sig.append((s1,s2))
    c.sig=sig
    # check theta sigma1 = -a (A-1)
    g=[]
    for j in range(3):
        aj=alpha[j]; s1,s2=sig[j]
        chk=[x+aj*(w-(1 if i==0 else 0)) for i,(x,w) in enumerate(zip(ytheta(s1),A))]
        assert all(x==0 for x in chk), "mirror map check"
        inner=ymul(yadd(ytheta(s2),yscal(-aj,s1)),Ainv,D)
        t=ytheta(inner); t[0]+=aj*aj
        g.append(ymul(t,Ainv,D))
    # Solve C2 a^2 - C1 a + C0 = g_j for three j, order by order
    C2=[Fr(0)]*(D+1);C1=[Fr(0)]*(D+1);C0=[Fr(0)]*(D+1)
    M=[[alpha[j]**2,-alpha[j],Fr(1)] for j in range(3)]
    def solve3(M,b):
        A_=[row[:]+[b[i]] for i,row in enumerate(M)]
        for i in range(3):
            p=next(r for r in range(i,3) if A_[r][i]!=0); A_[i],A_[p]=A_[p],A_[i]
            A_[i]=[x/A_[i][i] for x in A_[i]]
            for r in range(3):
                if r!=i:
                    f=A_[r][i]; A_[r]=[x-f*w for x,w in zip(A_[r],A_[i])]
        return [A_[i][3] for i in range(3)]
    for d in range(D+1):
        s=solve3(M,[g[j][d] for j in range(3)])
        C2[d],C1[d],C0[d]=s
    c.C2,c.C1,c.C0=C2,C1,C0
    c.Dfun=ymul(A,C2,D); c.e1E=ymul(A,C1,D)
    # weights
    wts=[]
    for i in range(3):
        wts.append([3*alpha[i]]+[alpha[l]-alpha[i] for l in range(3) if l!=i]+[-eps[0],-eps[1]])
    c.wts=wts
    W0=[]
    for i in range(3):
        p=Fr(1)
        for w in wts[i]: p*=w
        W0.append(p)
    c.W0=W0
    # eta
    eta=[[sum((-alpha[i])**(a+b)/W0[i] for i in range(3)) for b in range(3)] for a in range(3)]
    c.eta=eta
    # inverse
    M=[row[:]+[Fr(int(i==k)) for k in range(3)] for i,row in enumerate(eta)]
    for i in range(3):
        p=next(r for r in range(i,3) if M[r][i]!=0); M[i],M[p]=M[p],M[i]
        M[i]=[x/M[i][i] for x in M[i]]
        for r in range(3):
            if r!=i:
                f=M[r][i]; M[r]=[x-f*w for x,w in zip(M[r],M[i])]
    c.etainv=[row[3:] for row in M]
    # Y_j
    c.Y=[];c.Wval=[];c.unit=[];c.Delta=[]
    c.r=[];c.Pt=[]
    for j in range(3):
        aj=alpha[j]
        Y=[Fr(0)]*(D+1); Y[0]=-aj
        # Newton / fixed point: prod(Y+al)+27 y Y^3=0, solve order by order
        for n in range(1,D+1):
            # residual with Y[n]=0
            def Lhat(Yc):
                p=ymul(ymul(yadd(Yc,[alpha[0]]+[Fr(0)]*D),yadd(Yc,[alpha[1]]+[Fr(0)]*D),D),yadd(Yc,[alpha[2]]+[Fr(0)]*D),D)
                Y3=ymul(ymul(Yc,Yc,D),Yc,D)
                ys=[Fr(0)]+Y3[:-1] if D>=1 else [Fr(0)]
                return yadd(p,yscal(Fr(27),ys))
            res=Lhat(Y)[n]
            # derivative at y=0 of prod(Y+al): prod_{l!=j}(al-aj)
            dd=Fr(1)
            for l in range(3):
                if l!=j: dd*=alpha[l]-aj
            Y[n]=-res/dd
        c.Y.append(Y)
        Wv=yadd(yadd(yscal(e1,ymul(Y,Y,D)),yscal(2*e2,Y)),[3*e3]+[Fr(0)]*D)
        c.Wval.append(Wv)
        assert Wv[0]==aj*(alpha[(j+1)%3]-aj)*(alpha[(j+2)%3]-aj)
        c.Delta.append(yscal(3*eps[0]*eps[1],Wv))
        unit=ysqrt1(yscal(1/Wv[0],Wv),D)
        c.unit.append(unit)
        # shat as y-list of dict z-power
        shat=[{} for _ in range(D+1)]
        shat[0]={0:Fr(1)}
        for d in range(1,D+1):
            # N_d(z)=prod_{m=0}^{3d-1}(3a-mz); D_d=d! z^d prod_{l!=j}prod_m(al-a+m z)
            poly=[Fr(1)]
            for m in range(0,3*d):
                new=[Fr(0)]*(len(poly)+1)
                for k,x in enumerate(poly):
                    new[k]+=x*3*aj; new[k+1]+=x*(-m)
                poly=new
            ser=[Fr(0)]*(Z+d+2)
            for k,x in enumerate(poly):
                if k<len(ser): ser[k]=x
            # divide by each (b+m z)
            for l in range(3):
                if l==j: continue
                b=alpha[l]-aj
                for m in range(1,d+1):
                    new=[Fr(0)]*len(ser)
                    for k in range(len(ser)):
                        new[k]=(ser[k]-(m*new[k-1] if k>0 else 0))/b
                    ser=new
            ser=[x/factorial(d) for x in ser]
            shat[d]={k-d:x for k,x in enumerate(ser) if x!=0 and k-d<=Z}
        # exp(-uhat/z), uhat=sum Y_d y^d/d
        uh=[Fr(0)]+[Y[d]/d for d in range(1,D+1)]
        ex=[{} for _ in range(D+1)]; ex[0]={0:Fr(1)}
        # exp via recursion theta E=-theta(uhat)/z E -> d E_d=-sum_k Y_k z^-1 E_{d-k}
        for d in range(1,D+1):
            acc={}
            for k in range(1,d+1):
                for p,x in ex[d-k].items():
                    acc[p-1]=acc.get(p-1,Fr(0))+(-Y[k])*x
            ex[d]={p:x/d for p,x in acc.items()}
        r0=yzmul(ex,shat,D,Z-D-1)
        r1=yzop(r0,Y,D,Z-D-1)
        r1=yzmul(yz_from_y(Ainv),r1,D,Z-D-1)
        r2=yzop(r1,Y,D,Z-D-1,extra=yscal(-1,c.e1E))
        r2=yzmul(yz_from_y(yinv(c.Dfun,D)),r2,D,Z-D-1)
        rr=[r0,r1,r2]
        # regularity: no negative z powers
        for a in range(3):
            for d in range(D+1):
                for p,x in rr[a][d].items():
                    if p<0 and x!=0:
                        raise Exception("pole in r_%d at y^%d z^%d: %s"%(a,d,p,x))
        c.r.append(rr)
        # P~ := unit * r_a (without h)
        Pt=[yzmul(yz_from_y(unit),rr[a],D,Zmax+1) for a in range(3)]
        c.Pt.append(Pt)
    # full P: P^{(i)}_a=Pt*h_{wt_i}
    c.P=[]
    for i in range(3):
        hh=h_series(c.wts[i],Zmax+1,sign=(-1 if cflip else 1))
        hz=[{k:x for k,x in enumerate(hh) if x!=0}]+[{} for _ in range(D)]
        c.P.append([yzmul(c.Pt[i][a],hz,D,Zmax+1) for a in range(3)])
    return c

def checks(c):
    D=c.D; out={}
    # C0=0, r0|z=0 unit^-1, Psi_{a j}: r_a|z=0 =  lambda^a-like
    out['C0=0']=all(x==0 for x in c.C0)
    ok=True
    for j in range(3):
        v=[c.Pt[j][0][d].get(0,Fr(0)) for d in range(D+1)]
        ok&= v==[Fr(1)]+[Fr(0)]*D
    out['Pt0|z=0==1 (Delta=3e1e2 W)']=ok
    ok=True
    for j in range(3):
        for a in range(3):
            if c.r[j][a][0]!={0:(-c.alpha[j])**a}: ok=False
    out['r|y=0 z-free']=ok
    # symplectic
    bad=0
    for i in range(3):
        for j in range(3):
            tot=[{} for _ in range(D+1)]
            for a in range(3):
                for b in range(3):
                    Pjm=[{p:x*(1 if p%2==0 else -1) for p,x in c.Pt[j][b][d].items()} for d in range(D+1)]
                    m=yzmul(c.Pt[i][a],Pjm,D,c.Zmax)
                    for d in range(D+1):
                        for p,x in m[d].items():
                            tot[d][p]=tot[d].get(p,Fr(0))+c.etainv[a][b]*x
            for d in range(D+1):
                for p,x in tot[d].items():
                    exp=Fr(0)
                    if p==0 and i==j: exp=c.Delta[i][d]
                    if x!=exp: bad+=1
    out['symplectic defects']=bad
    return out

if __name__=="__main__":
    c=setup((Fr(1),Fr(2),Fr(5)),(Fr(3),Fr(7)),3,8)
    print(checks(c))
    print('D',c.Dfun,'e1E',c.e1E)
    print('Delta0',c.Delta[0])
