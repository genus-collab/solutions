# Evidence-only checks for the port of conifold s.13.1 to local P2 (c -> E).
from fractions import Fraction as F
from math import factorial
from itertools import product

# (1) Genus-0 test of the string / lam_0 identities with general leaf scale c (trivial CohFT, r >= 3 so no exceptional terms):
#     Xi(1,R) = J_R/c^2 Xi(R),  Xi(0,R) = (D_R - 2h + 2)/c^2 Xi(R)  with h = 0.
def coeff_leaf(j, a, c):   # [psi^a] 1/(c - j psi)
    return F(j)**a / c**(a+1)
def coeff_T(d, c):         # T = -psi^2/(c-psi) = -sum_{d>=2} psi^d / c^(d-1)
    return F(0) if d < 2 else -F(1)/c**(d-1)
def corr0(powers):         # <prod tau_a>_0 = (n-3)!/prod a!
    n = len(powers)
    if n < 3 or sum(powers) != n-3: return F(0)
    r = F(factorial(n-3))
    for a in powers: r /= factorial(a)
    return r
def compositions(total, k, mins):
    if k == 0:
        if total == 0: yield ()
        return
    for a in range(mins[0], total+1):
        for rest in compositions(total-a, k-1, mins[1:]):
            yield (a,)+rest
def Xi0(R, c):
    r = len(R); tot = F(0)
    for m in range(0, max(r-2,0)+1):
        n = r+m; dim = n-3
        if dim < 0: continue
        s = F(0)
        for pw in compositions(dim, n, [0]*r+[2]*m):
            term = corr0(list(pw))
            if term == 0: continue
            for i,j in enumerate(R): term *= coeff_leaf(j, pw[i], c)
            for k in range(m): term *= coeff_T(pw[r+k], c)
            s += term
        tot += s/factorial(m)
    return tot
c = F(7,3)
ok = True
for R in [(2,3,5),(2,2,4,3),(3,3,3),(2,5,7,1),(4,2,2,2)]:
    J = sum(R); D = J-len(R)
    x = Xi0(R,c)
    s1 = Xi0((1,)+R, c); s0 = Xi0((0,)+R, c)
    a = (s1 == F(J)/c**2 * x); b = (s0 == F(D+2)/c**2 * x)
    ok &= a and b
    print("R=",R,"Xi=",x," string:",a," lam0:",b)
# control: wrong scale power must fail
R=(2,3,5); x=Xi0(R,c); print("control (J/c instead of J/c^2):", Xi0((1,)+R,c) == F(10)/c*x)
print("genus-0 string/lam0 identities with c = 7/3:", ok)

# (2) gap bookkeeping: lp2 single-vertex pole predicted = (3 e1 e2 E^3)^(g-1) Z[g,0];
#     with e2 = 1, e1 = x, E = c this is (-3)^(g-1) (-x c^3)^(g-1) Z[g,0] = (-3)^(g-1) gap_g(conifold).
#     Compare with lp2 log value gap_3 = -169/27216 at x0 = 1/3, via (2g-3)! [t^(2g)] t^2/(4 sinh(e1 t/2) sinh(e2 t/2)).
def sinh_series(a, N):  # coefficients of sinh(a t) up to t^N
    return [F(0) if k%2==0 else a**k/factorial(k) for k in range(N+1)]
def mul(p,q,N):
    r=[F(0)]*(N+1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):
            if i+j<=N: r[i+j]+=x*y
    return r
def inv_shifted(p, N):  # p = t^2 * u(t), return series of 1/u up to t^N
    u = p[2:]+[F(0)]*2
    inv=[F(0)]*(N+1); inv[0]=1/u[0]
    for n in range(1,N+1):
        inv[n] = -sum(u[k]*inv[n-k] for k in range(1,n+1) if k < len(u))/u[0]
    return inv
for g, x0, logged in [(3, F(1,3), F(-169,27216)), (3, F(4,5), F(-1414951,224000000))]:
    N = 2*g+2
    prod4 = [4*v for v in mul(sinh_series(x0/2,N+2), sinh_series(F(1,2),N+2), N+2)]
    inv = inv_shifted(prod4, N)       # t^2/(4 sinh sinh) series
    coef = inv[2*g]
    cf = factorial(2*g-3)*coef        # conifold gap_g at e1 = x0, e2 = 1
    pred = F(-3)**(g-1)*cf
    print("g=",g,"x0=",x0," conifold gap =",cf," (-3)^(g-1) gap =",pred," lp2 logged:",logged, pred==logged)
# scaling of Delta prefactor: 3 e1 e2 E^3 / (-x c^3) at e2=1, e1=x, E=c
x0=F(5,7); print("3 e1 e2 E^3/(-x c^3) =", 3*x0*(1+x0)**3/(-(x0*(1+x0)**3)))
print("--- per-colouring logs, single-vertex poles ---")
for g, x0, logged in [(2,F(1,3),F(-83,6480)),(2,F(1,8),F(-9773,327680)),(3,F(1,8),F(-1666067,117440512))]:
    N=2*g+2
    prod4=[4*v for v in mul(sinh_series(x0/2,N+2), sinh_series(F(1,2),N+2), N+2)]
    inv=inv_shifted(prod4,N); pred=F(-3)**(g-1)*factorial(2*g-3)*inv[2*g]
    print("g=",g,"x0=",x0,"pred",pred,"logged",logged,pred==logged)
