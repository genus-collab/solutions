# E6: check that R = exp(r), r = sum_k r_k z^k with r_k^T = (-1)^(k+1) r_k (r_odd symmetric, r_even antisymmetric)
# satisfies R(-z)^T R(z) = 1 mod z^(N+1), and conversely that log of a random symplectic R has that shape.
# Exact arithmetic over Q (fractions). Matrices 3x3, truncation N = 5.
from fractions import Fraction as F
import random
random.seed(7)
n, N = 3, 5

def zero(): return [[F(0)]*n for _ in range(n)]
def eye(): return [[F(int(i==j)) for j in range(n)] for i in range(n)]
def add(A,B): return [[A[i][j]+B[i][j] for j in range(n)] for i in range(n)]
def scal(c,A): return [[c*A[i][j] for j in range(n)] for i in range(n)]
def mul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def tr(A): return [[A[j][i] for j in range(n)] for i in range(n)]
# series: list of N+1 matrices
def smul(S,T): return [ (lambda k: __import__('functools').reduce(add,[mul(S[a],T[k-a]) for a in range(k+1)],zero()))(k) for k in range(N+1)]
def sadd(S,T): return [add(S[k],T[k]) for k in range(N+1)]
def sscal(c,S): return [scal(c,S[k]) for k in range(N+1)]
def sexp(X):  # X has zero constant term
    res=[eye()]+[zero() for _ in range(N)]; term=[eye()]+[zero() for _ in range(N)]
    for j in range(1,N+1):
        term=sscal(F(1,j),smul(term,X)); res=sadd(res,term)
    return res
def slog(R):  # R(0)=1
    Y=[zero()]+R[1:]; res=[zero() for _ in range(N+1)]; p=[eye()]+[zero() for _ in range(N)]
    for j in range(1,N+1):
        p=smul(p,Y); res=sadd(res,sscal(F((-1)**(j+1),j),p))
    return res
def sneg_t(S): return [scal(F((-1)**k),tr(S[k])) for k in range(N+1)]  # S(-z)^T
def rnd(): return [[F(random.randint(-4,4),random.randint(1,3)) for j in range(n)] for i in range(n)]
r=[zero()]
for k in range(1,N+1):
    A=rnd(); r.append(add(A,tr(A)) if k%2 else add(A,scal(F(-1),tr(A))))
R=sexp(r)
P=smul(sneg_t(R),R)
ok1 = all(P[k]==(eye() if k==0 else zero()) for k in range(N+1))
print("exp(r) symplectic mod z^%d:"%(N+1), ok1)
# converse: build symplectic R as product of two exps, take log, test shape
r2=[zero()]
for k in range(1,N+1):
    A=rnd(); r2.append(add(A,tr(A)) if k%2 else add(A,scal(F(-1),tr(A))))
R2=smul(R,sexp(r2))
P2=smul(sneg_t(R2),R2)
ok2 = all(P2[k]==(eye() if k==0 else zero()) for k in range(N+1))
L=slog(R2)
ok3 = all(tr(L[k])==(L[k] if k%2 else scal(F(-1),L[k])) for k in range(1,N+1))
print("product symplectic:", ok2, " log has r_k^T=(-1)^(k+1) r_k:", ok3)
