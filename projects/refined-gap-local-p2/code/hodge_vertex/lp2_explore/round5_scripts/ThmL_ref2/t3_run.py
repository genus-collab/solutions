# Relation (4.1), G>=3, Hodge-side values, with string stripping of unit parts in the value function
# (string itself is tested directly on small cases first). Controls that must fail are run for each relation.
import sys,time
import t2_relation as T2
from t2_relation import *
import hodge
_cache={}
def Xstrip(g,U,tau):
    U=tuple(sorted(U)); key=(g,U,tau)
    if key in _cache: return _cache[key]
    Up=tuple(u for u in U if u>1); v=len(U)-len(Up)
    if v==0 or g==0:
        val=hodge.X(g,U,tau) if g>0 else X0(U)
    elif Up:
        J=sum(Up); r=Fr(1)
        for i in range(v): r*=(J+i)
        val=r*hodge.X(g,Up,tau)
    elif g==1:
        val=Fr(factorial(v-1))*hodge.X(1,(1,),tau)
    else:
        val=Fr(0)
    _cache[key]=val
    return val
T2.Xg=lambda g,U,tau: Xstrip(g,U,tau)
# string direct test
t=Fr(3,5)
for g,U in [(2,(3,1)),(2,(2,3,1)),(3,(5,1)),(3,(5,1,1)),(1,(1,1,1)),(1,(4,1)),(3,(1,)),(3,(1,1,1)),(2,(3,3,1,1))]:
    assert hodge.X(g,U,t)==Xstrip(g,U,t),(g,U)
print("string stripping agrees with direct on 9 cases",flush=True)
tau=Fr(sys.argv[1]); G=int(sys.argv[2]); kmax=int(sys.argv[3]); rcap=int(sys.argv[4])
cases=[(2,()),(3,()),(4,()),(5,()),(6,()),(7,()),(8,()),(9,()),
       (2,(2,)),(2,(3,)),(2,(9,)),(3,(5,)),(5,(5,)),(4,(3,)),(6,(17,)),(3,(21,)),(5,(13,)),
       (2,(3,3)),(2,(2,2)),(4,(2,7)),(3,(8,13)),(2,(2,2,2)),(3,(2,3,4)),(2,(5,5,5))]
t0=time.time(); nt=0;fails=0
ctl={k:[0,0] for k in ("pert","rising","beta","cwrong","even")}
for k0,R in cases:
    if k0>kmax or len(R)>rcap: continue
    for n in (2*G,2*G+1):
        v=relation(k0,n,R,tau); nt+=1
        if v!=0: fails+=1
        out=[]
        if n==2*G:
            T=tuple(sorted(R+(k0,)))
            o=Opts(); o.pert=(G,T); w=relation(k0,n,R,tau,o); ctl["pert"][1]+=1; ctl["pert"][0]+= (w!=0)
            # response equals -(k0-1) Pi sigma^G tau delta for k0 not in R, single target
            if k0 not in R:
                D=sum(T)-len(T); delta=D-2*G+2; s=tau*(1+tau)
                pred=-(k0-1)*Pi(k0,tau)*s**G*tau*delta
                out.append(("resp_ok",w==pred or delta==0))
            o=Opts(); o.fall=False; w=relation(k0,n,R,tau,o); ctl["rising"][1]+=1; ctl["rising"][0]+=(w!=0)
            o=Opts(); o.beta_sign=-1; w=relation(k0,n,R,tau,o); ctl["beta"][1]+=1; ctl["beta"][0]+=(w!=0)
            o=Opts(); o.cwrong=True; w=relation(k0,n,R,tau,o); ctl["cwrong"][1]+=1; ctl["cwrong"][0]+=(w!=0)
            o=Opts(); o.evenwrong=True; w=relation(k0,n,R,tau,o); ctl["even"][1]+=1; ctl["even"][0]+=(w!=0)
        print("k0",k0,"R",R,"n",n,"|mu|",k0+sum(R),"value",v,out,"t",round(time.time()-t0,1),flush=True)
print("tau",tau,"G",G,"relations",nt,"failures",fails,"controls nonzero/total",ctl)
