from fractions import Fraction as Fr
from check_pf_sympl import *
import time
t0=time.time()
print("PF z^100, E=29/21:", pf_check(Fr(29,21),100), round(time.time()-t0,1))
E=Fr(-7,10)
for (t,i,j) in [(Fr(-1,7),0,0),(Fr(-1,7),1,1),(Fr(-1,7),0,1),(Fr(3,11),1,2),(Fr(-1,3)+Fr(1,1000),0,0)]:
    out=pair_tilde(E,t,70,i,j)
    nz=[n for n,x in enumerate(out) if not Qw._c(x).iszero()]
    print("Pitilde z^70",t,i,j,nz[:4],out[0],round(time.time()-t0,1))
