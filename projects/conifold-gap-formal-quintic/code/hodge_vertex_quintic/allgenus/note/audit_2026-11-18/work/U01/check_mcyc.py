# Check the operator identity Mcyc^LR = kappa3 * PF4 on test functions (incl. l-dependent ones),
# and the idempotent/Delta formulas of cor:set-cyclicPF5(a) via the structure constants.
exec(open('check_kappa.py').read().split('# Birkhoff')[0])
C0inv = s_inv(C0)
T = p_trim([s_mul(I1t,C0inv), s_one()])
C1 = p_D(T)[0]; C1inv = s_inv(C1)
L = s_zero()
for n in range(N+1):
    c = Fr(1)
    for j in range(n): c *= (Fr(-1,5)-j)
    c /= factorial(n)
    L[n] = c*(-5**5)**n
k1 = s_mul(s_pow(L,5), s_inv(s_mul(s_pow(C0,2), s_pow(C1,3))))
k3 = s_mul(C0, C1inv)
def dT(p): return p_scal_series(C1inv, p_D(p))
def Mcyc(p):
    x = p_scal_series(C0inv, p)
    x = dT(dT(x))
    x = p_scal_series(s_inv(k1), x)
    return dT(dT(x))
import random
random.seed(1)
tests = []
for t in range(3):
    f = [Fr(random.randint(-9,9)) for _ in range(N+1)]
    g = [Fr(random.randint(-9,9)) for _ in range(N+1)]
    tests.append([f])            # l-free
    tests.append([f,g])          # degree 1 in l
    tests.append([f,g,f])        # degree 2 in l
ok = True
for p in tests:
    a = Mcyc(p); b = p_scal_series(k3, PF4(p))
    if p_trim(a) != p_trim(b): ok = False
print("Mcyc^LR == kappa3*PF4 on test functions:", ok)
