from fractions import Fraction as F
import cmath, math
z5 = cmath.exp(2j*math.pi/5)
# power sums over w_0 = {1-zeta^j}_{j=1..4} U {-5}
for k in (1,2,3):
    s = sum((1-z5**j)**(1-2*k) for j in range(1,5)) + (-5)**(1-2*k)
    print("k",k,"sum a^(1-2k) over w0:",s)
# Bernoulli numbers
def bern(n):
    B=[F(1)]
    for m in range(1,n+1):
        B.append(-sum(math.comb(m+1,k)*B[k] for k in range(m))/F(m+1))
    return B
B=bern(12)
b=[None]+[B[2*k]/(2*k*(2*k-1)) for k in range(1,6)]
print("b1,b2:",b[1],b[2])
# h1 with -5 and +5
s1m = F(2)+F(-1,5); s1p=F(2)+F(1,5)
print("h1(-5):",-b[1]*s1m," h1(+5):",-b[1]*s1p)
# z^3 coefficient of log h: -b2 * sum a^-3 ; sum_{j}(1-z^j)^-3 = -1 ; (-5)^-3
s3 = F(-1)+F(-1,125)
print("log h z^3 coeff:", -b[2]*s3, " GJR17:", F(1,360)*(F(1,(-5)**3)-1))
# Bernoulli polynomial B6(1/5)/6
def bernpoly(n,x):
    return sum(math.comb(n,k)*B[k]*x**(n-k) for k in range(n+1))
v = bernpoly(6,F(1,5))/6
print("B6(1/5)/6 =", v, " target 4537/3937500 =", F(4537,3937500))
# G_inf coefficient check from (37): m*(-1)^{mj+1}B_{mj+1}(k/m)/(mj(mj+1)) at m=5,j=1,k=1
m=5;j=1
c = m*(-1)**(m*j+1)*bernpoly(m*j+1,F(1,m))/(m*j*(m*j+1))
print("(37) coefficient z^5:", c)
