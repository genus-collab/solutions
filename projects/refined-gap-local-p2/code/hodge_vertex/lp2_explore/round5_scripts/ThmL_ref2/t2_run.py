import sys,time
from t2_relation import *
tau=Fr(sys.argv[1]); G=int(sys.argv[2])
sets=[]
# k0, R : large parts beyond the author's |mu|<=20 window
cases=[(2,()),(3,()),(4,()),(5,()),(7,()),(9,()),(11,()),(13,()),
       (2,(2,)),(2,(9,)),(3,(5,)),(5,(5,)),(7,(11,)),(9,(13,)),(6,(17,)),(12,(23,)),
       (2,(3,3)),(4,(2,7)),(3,(8,13)),(6,(6,11)),(5,(19,4))]
t0=time.time()
fails=0;n_tot=0;nontriv=0
for k0,R in cases:
    if k0>int(sys.argv[3]): continue
    for n in (2*G,2*G+1):
        v=relation(k0,n,R,tau)
        # nontrivial content: compare with perturbation
        n_tot+=1
        if v!=0: fails+=1
        print("k0",k0,"R",R,"n",n,"|mu|",k0+sum(R),"value",v,"t",round(time.time()-t0,1),flush=True)
print("tau",tau,"G",G,"relations",n_tot,"failures",fails)
