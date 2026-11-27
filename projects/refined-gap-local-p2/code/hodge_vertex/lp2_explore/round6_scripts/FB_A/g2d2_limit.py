# Evidence only. y^2 coefficient of the genus-2 graph sum at weights near the diagonal, versus the refined value at T_F.
# GW_{2,2} = F[2] + 6 F[1] because Q = y(1 - 6y + ...). Expected limit from refined BPS data of local P2 (d=1: (0,1); d=2: (0,5/2)).
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2d1 import *

def ser_coef_deg2(num_fn, e1, e2, E):
    return None

def poly_deg(fn_terms, N):
    pass

def refined_d2(e1,e2):
    # degree-2 part (in eps) of -chi_{5/2}(e^{E/2})/(4 sinh(e1/2) sinh(e2/2)) + (1/2) chi_1(e^{E})/(4 sinh(e1) sinh(e2))
    E=e1+e2
    def cosh_sum(vals,order):
        # sum over exponents a of exp(a) truncated to degree 'order': returns dict deg->sum a^deg/deg!
        return {d:sum(Fr(a)**d for a in vals)/math.factorial(d) for d in range(order+1)}
    # 1/(4 sinh(a/2) sinh(b/2)) = (1/ab)(1 - (a^2+b^2)/24 + 7(a^4+b^4)/5760 + a^2b^2/576 + ...)
    def invsinh(a,b):
        return {0:Fr(1),2:-(a*a+b*b)/24,4:Fr(7,5760)*(a**4+b**4)+a*a*b*b/576}
    # primitive part: exponents E*(5/2,3/2,1/2,-1/2,-3/2,-5/2)
    vals=[E*Fr(k,2) for k in (5,3,1,-1,-3,-5)]
    ch=cosh_sum(vals,4)
    ins=invsinh(e1,e2)
    prim=-sum(ch[d]*ins[4-d] for d in (0,2,4))/(e1*e2)
    # double cover: chi_1(v^2): exponents E*(2,0,-2) (v=e^{E/2}, v^2 = e^E: chi_1(v^2)= v^4+1+v^-4 = e^{2E}+1+e^{-2E})
    vals2=[2*E,Fr(0),-2*E]
    ch2=cosh_sum(vals2,4)
    ins2=invsinh(2*e1,2*e2)
    dc=Fr(1,2)*sum(ch2[d]*ins2[4-d] for d in (0,2,4))/(4*e1*e2)
    return prim+dc

if __name__=="__main__":
    for e1,e2 in [(Fr(3),Fr(7)),(Fr(1),Fr(1)),(Fr(5,2),Fr(-3,7))]:
        E=e1+e2
        print("eps",e1,e2,"refined d=2 value",float(refined_d2(e1,e2)),"; refined d=1 value",float(refined_d1(e1,e2)))
        for k in [2,3,4]:
            dl=Fr(1,10**k)
            al=(E/3+dl,E/3+2*dl-3*dl*dl,E/3-3*dl+dl*dl)
            T=Theory(al,e1,e2)
            F=T.F2()
            print("  delta 1e-%d: F[1] = %.12f ; F[2]+6F[1] = %.12f"%(k,float(F[1]),float(F[2]+6*F[1])))
