# fresh: genus-2 d=0 from W_i^0 * Int_{Mbar_2} Hodge_{wt_i} = W*(p3-p1^3)/17280 (lambda1^3=1/2880, lambda1 lambda2=1/5760)
from fractions import Fraction as F
def f(al,e1,e2):
    tot=F(0)
    for i in range(3):
        j,k=[x for x in range(3) if x!=i]
        w=[F(3*al[i]),F(al[j]-al[i]),F(al[k]-al[i]),F(-e1),F(-e2)]
        W=F(1)
        for a in w: W*=a
        p1=sum(1/a for a in w); p3=sum(1/a**3 for a in w)
        tot+=W*(p3-p1**3)/17280
    return tot
print(f((1,2,5),3,7), F(19189,100800))
print(f((2,-3,7),1,4), F(158713,290304))
print(f((1,3,8),F(5,2),F(-3,7)), F(-10457963,29030400))
