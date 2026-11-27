from schur_side import *
from hodge import *
K=Fld(None)
for tau in [Fr(2,7),Fr(-5,3),Fr(13,9)]:
    z=z_of_tau(tau); sig=tau*(1+tau)
    for h in (2,3,4):
        N=3*h-3;Pp=2*h-2+2*N+1
        F,items=compute_F((),N,Pp,tau,K)
        tot=Fr(0)
        for a in range(1,N+1):
            tot+=(-1)**a*(1/sig)**a*laguerre(N-a,2*h-2+a,-a)*coeff(F,items,(),a,2*h-2+a)
        tot=tot/sig**(h-1)+laguerre(N,2*h-2,0)*peval(int_omega_poly(h,()),z)
        print(tau,h,tot==peval(Xi_poly(h,()),z), tot)
