from hodge import *
def mul(A,B):
    out={}
    for a,x in A.items():
        for b,y in B.items():
            k=tuple(sorted(a+b)); out[k]=out.get(k,Fr(0))+x*y
    return out
def add(*As):
    out={}
    for A,c in As:
        for k,v in A.items(): out[k]=out.get(k,Fr(0))+c*v
    return out
L0={():Fr(1)}; L1={(1,):Fr(1)}; L2={(1,1):Fr(1,2)}; L3={(1,1,1):Fr(1,6),(3,):Fr(2)}
Lam=add((L0,1),(L1,-1),(L2,1),(L3,-1))
e1=add((L1,1),(mul(L1,L1),1),(L2,-2))
e2=add((L2,1),(mul(L1,L2),1),(L3,-3),(mul(L2,L2),1),(mul(L1,L3),-2))
e3=mul(L3,add((L0,1),(L1,1),(L2,1),(L3,1)))
E={0:L0,1:e1,2:e2,3:e3}
g=3
bad=0;n=0
for R in [(2,),(3,),(2,2),(4,),(3,2),(2,2,2),(5,),(4,2),(3,3),(3,2,2),(6,),(5,2),(4,3),(1,2),(1,1,3),(1,1)]:
    r=len(R);N=3*g-3+r
    res=[]
    for k in range(4):
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
                        if sum(ch)==Mleft-M1 and (not ch or max(ch)<=2*g-1):
                            tot+=w*wt*cf*chint(g,tuple(sorted(ch)),ds)
        res.append(tot)
    p=Xi_poly(g,R); p=p+[Fr(0)]*(4-len(p))
    n+=1
    if p!=res: bad+=1
    print(R,"D=",sum(R)-len(R),[str(x) for x in res],"match" if p==res else "MISMATCH")
print("genus3 direct e_k: cases",n,"bad",bad)
