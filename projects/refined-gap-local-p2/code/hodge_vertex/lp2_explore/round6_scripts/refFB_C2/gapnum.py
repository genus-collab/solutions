# Exact evaluation of G_g[F_cf] at rational t near -1/3, compared with (-3)^(g-1) beta_(2g-2)/t_c^(2g-2).
# Per-graph, per-colouring bounds as delta -> 0. Evidence only.
import sys, math, itertools, time
from fractions import Fraction as Fr
from tf import Closed, TF, stable
from hod import B

def cf_series(N):
    # wc = 2F1(1/3,2/3;1;xi)
    w=[Fr(1)]
    for n in range(N):
        w.append(w[n]*(Fr(1,3)+n)*(Fr(2,3)+n)/((n+1)**2))
    # t_c = Int wc/(1-xi)
    p=[Fr(0)]*(N+1); s=Fr(0)
    for n in range(N+1):
        s+=w[n]; p[n]=s
    tc=[Fr(0)]+[p[n]/(n+1) for n in range(N)]
    return w,tc
def ev(c,x):
    r=Fr(0)
    for a in reversed(c): r=r*x+a
    return r
def dev(c,x):
    r=Fr(0)
    for n in range(len(c)-1,0,-1): r=r*x+n*c[n]
    return r
def rnd(x,digits=120):
    d=10**digits
    return Fr(int(x*d),d)
def beta(n,e1,e2):
    # (n-1)! [t^n] 1/(4 sinh(e1 t/2) sinh(e2 t/2))
    # x/(2 sinh(x/2)) = Sum_k c_k x^(2k), c_k = B_{2k}(1/2)/(2k)!, B_n(1/2)=(2^{1-n}-1)B_n
    N=(n+2)//2+1
    c=[(Fr(2)**(1-2*k)-1)*B[2*k]/math.factorial(2*k) for k in range(N+1)]
    tot=Fr(0)
    # 1/(4 sinh sinh)=(1/(e1 e2 t^2)) Sum c_i c_j e1^{2i} e2^{2j} t^{2i+2j}
    for i in range(N+1):
        for j in range(N+1):
            if 2*i+2*j-2==n: tot+=c[i]*c[j]*Fr(e1)**(2*i)*Fr(e2)**(2*j)
    return tot/(Fr(e1)*Fr(e2))*math.factorial(n-1)

if __name__=="__main__":
    g=int(sys.argv[1])
    deltas=[Fr(1,10**k) for k in (2,3,4,5)]
    if len(sys.argv)>2: deltas=deltas[:int(sys.argv[2])]
    w,tc=cf_series(60)
    print("beta2(3,7)",beta(2,3,7),"check",Fr(7*81+10*9*49+7*2401,5760*21))
    for e1,e2 in ((3,7),(1,2),(Fr(5,2),Fr(-3,7))):
        E=Fr(e1)+Fr(e2); wts=[E,-Fr(e1),-Fr(e2)]
        cl=Closed(E,3*g+1)
        pc=(-3)**(g-1)*beta(2*g-2,e1,e2)
        print(f"g={g} eps=({e1},{e2}) polar coefficient {pc}")
        last=None
        graphs=stable(g)
        for d in deltas:
            t=Fr(-1,3)+d
            xi=1+27*t**3
            wv=rnd(ev(w,xi)); tcv=rnd(ev(tc,xi)); dw=rnd(dev(w,xi))
            Fcf=rnd((xi-1)*dw/wv); c1=-wv
            t0=time.time()
            hv=TF(cl,e1,e2,t,Fcf,c1,'hv',0,wts)
            tot,per=hv.total(g,per=True)
            diff=tot.a-pc/tcv**(2*g-2)
            # max colouring size per graph with an edge
            maxc=0.0
            if g<=3:
                for gen,ed,aut in graphs:
                    if not ed: continue
                    for col in itertools.product(range(3),repeat=len(gen)):
                        x=hv.graph(gen,ed,[col])
                        maxc=max(maxc,abs(float(x.a))+abs(float(x.b)))
            # vertex graph, colour by colour
            vg=[(gen,ed) for gen,ed,a in graphs if not ed][0]
            vc=[hv.graph(vg[0],vg[1],[(c,)]) for c in range(3)]
            vdiff0=vc[0].a-pc/tcv**(2*g-2)
            print(f"  delta={d}: xi={float(xi):.3e} tc={float(tcv):.6e} Im={float(tot.b):.1e} G-polar={float(diff):.10f}; colour0 vertex-polar={float(vdiff0):.10f}; col1 vertex={float(vc[1].a):.6f}{float(vc[1].b):+.6f}w; max|T(G,c)| over edge graphs={maxc:.6f} [{time.time()-t0:.1f}s]",flush=True)
            # controls: F_cf shifted by constant; c1-solution wrong; t_c wrong normalisation
            bad=TF(cl,e1,e2,t,Fcf+Fr(1,9),c1,'hv',0,wts).total(g)
            bad2=TF(cl,e1,e2,t,Fcf+xi/5,c1,'hv',0,wts).total(g)
            print(f"     controls: F+1/9 -> G-polar={float(bad.a-pc/tcv**(2*g-2)):.6f}; F+xi/5 -> {float(bad2.a-pc/tcv**(2*g-2)):.6f}; polar with xi instead of t_c: {float(tot.a-pc/xi**(2*g-2)):.6f}",flush=True)
