# Target coefficient of s.5, including k0 repeated in R (report tested only k0 not in R), and every choice of k0 in T.
from t2_relation import *
import t2_relation as T2
tests=0;bad=0
for tau in [Fr(3,5),Fr(-5,3)]:
    s=tau*(1+tau)
    for G,T in [(3,(2,2)),(3,(2,2,2)),(3,(3,3)),(3,(2,3)),(3,(4,)),(2,(2,)),(2,(2,2)),(4,(2,2)),(4,(3,3)),(4,(2,2,3)),(4,(4,4)),(3,(5,5)),(4,(2,5))]:
        D=sum(T)-len(T); delta=D-2*G+2
        for k0 in set(T):
            L=list(T); L.remove(k0); R=tuple(L)
            o=Opts(); o.pert=(G,tuple(sorted(T)))
            base=relation(k0,2*G,R,tau)
            w=relation(k0,2*G,R,tau,o)
            pred=-(k0-1)*Pi(k0,tau)*s**G*tau*delta
            tests+=1
            flag=(base==0 and w-base==pred)
            if not flag: bad+=1
            print("G",G,"T",T,"k0",k0,"delta",delta,"base",base,"resp==pred",w-base==pred)
print("tests",tests,"failures",bad)
