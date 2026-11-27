from hodge import *
tau=Fr(2,7); z=z_of_tau(tau)
print("Xi_4(7)",peval(Xi_poly(4,(7,)),z),"doc 193375/214990848")
for h in range(2,8):
    # sharpness: largest D with vanishing; single part and top coefficient
    print(h,"Xi(2h-1) z-poly leading", [str(c) for c in Xi_poly(h,(2*h-1,))][-1], "beta-like?")
# check e_k expansions directly at h=2: e_1(alpha+alpha^2)=lambda1+lambda1^2-2lambda2 ; e_2=(a1+a1^2)(a2+a2^2)=lambda2(1+lambda1+lambda2)=lambda2+lambda1lambda2(+lambda2^2=0)
# Lambda(-1)=1-lambda1+lambda2. lambda1=ch1, lambda2=ch1^2/2.
# Xi^(k)_2(R) direct via lambda polynomial
def int_class(g,classes,ds):
    # classes: dict chs tuple->coef
    return sum(c*chint(g,tuple(sorted(ch)),tuple(sorted(ds))) for ch,c in classes.items())
def mul(A,B):
    out={}
    for a,x in A.items():
        for b,y in B.items():
            k=tuple(sorted(a+b)); out[k]=out.get(k,Fr(0))+x*y
    return out
L0={():Fr(1)}; L1={(1,):Fr(1)}; L2={(1,1):Fr(1,2)}
def add(*As):
    out={}
    for A,c in As:
        for k,v in A.items(): out[k]=out.get(k,Fr(0))+c*v
    return out
Lam=add((L0,1),(L1,-1),(L2,1))
e1=add((L1,1),(mul(L1,L1),1),(L2,-2))
e2=add((L2,1),(mul(L1,L2),1))
E={0:L0,1:e1,2:e2}
for R in [(2,),(3,),(2,2),(4,)]:
    r=len(R);N=3*2-3+r
    res=[]
    for k in range(3):
        cls=mul(Lam,E[k])
        tot=Fr(0)
        for e,w in leaf_groups(R,N).items():
            Mleft=N-sum(e)
            for M1 in range(Mleft+1):
                for par in partitions_ge(M1,M1 if M1>0 else 1):
                    ds=tuple(sorted(e+tuple(p+1 for p in par)))
                    wt=Fr((-1)**len(par))
                    for v,c in Counter(par).items(): wt/=factorial(c)
                    for ch,cf in cls.items():
                        if sum(ch)==Mleft-M1:
                            tot+=w*wt*cf*chint(2,tuple(sorted(ch)),ds)
        res.append(tot)
    print("R",R,"direct Xi^(k)",[str(x) for x in res],"from z-poly",[str(x) for x in Xi_poly(2,R)])
