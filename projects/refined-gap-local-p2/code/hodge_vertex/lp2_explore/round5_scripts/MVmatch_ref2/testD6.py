import hodge
hodge.BB[:] = hodge.bern(160)
from testD import *
for h,R in [(6,(10,)),(6,(5,5)),(5,(4,4))]:
    print(h,R,"D",D(R),Xi_schur(h,R,F(3,2)),flush=True)
