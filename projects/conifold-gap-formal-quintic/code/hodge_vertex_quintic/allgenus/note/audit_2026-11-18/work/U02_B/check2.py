# Check of Psi_0, M_0 symplecticity (lem:set-M0sympl), unit (lem:set-dilaton (ii)) to low order, in Q(zeta).
from fractions import Fraction as Fr
import json
M = 8
Fs = [[Fr(c) for c in F] for F in json.load(open("Fs.json"))["F"]]

# Q(zeta): vectors of length 4 in basis 1,z,z^2,z^3 with z^4 = -1-z-z^2-z^3
class Z:
    __slots__ = ("c",)
    def __init__(self, c):
        self.c = [Fr(x) for x in c]
    @staticmethod
    def zeta(k):
        k %= 5
        if k < 4:
            c = [0, 0, 0, 0]; c[k] = 1
            return Z(c)
        return Z([-1, -1, -1, -1])
    def __add__(self, o):
        o = o if isinstance(o, Z) else Z([o, 0, 0, 0])
        return Z([a + b for a, b in zip(self.c, o.c)])
    __radd__ = __add__
    def __neg__(self):
        return Z([-a for a in self.c])
    def __sub__(self, o):
        return self + (-(o if isinstance(o, Z) else Z([o, 0, 0, 0])))
    def __mul__(self, o):
        if not isinstance(o, Z):
            return Z([a * o for a in self.c])
        r = [Fr(0)] * 5
        for i, a in enumerate(self.c):
            if a == 0: continue
            for j, b in enumerate(o.c):
                r[(i + j) % 5] += a * b
        # reduce z^4 = -(1+z+z^2+z^3)
        return Z([r[k] - r[4] for k in range(4)])
    __rmul__ = __mul__
    def iszero(self):
        return all(a == 0 for a in self.c)
    def __eq__(self, o):
        o = o if isinstance(o, Z) else Z([o, 0, 0, 0])
        return (self - o).iszero()
    def __repr__(self):
        return str([str(a) for a in self.c])

eps = [1, -1, -1, 1, 1]
nu = [0, -1, 3, 2, 1]
# Rhat[b][i] = series in t: zeta^{i(b+1)} F_b(zeta^{-i} t)
Rhat = [[[Z.zeta(i * (b + 1) - i * n) * Fs[b][n] for n in range(M)] for i in range(5)] for b in range(5)]

def smul(a, b):
    r = [Z([0, 0, 0, 0]) for _ in range(M)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < M:
                r[i + j] = r[i + j] + x * y
    return r
def sneg(a):  # a(-t)
    return [x * ((-1) ** n) for n, x in enumerate(a)]

# N_ij(t) = (1/5) sum_a Rhat_a^{(i)}(-t) Rhat_{3-a}^{(j)}(t)
ok = True
for i in range(5):
    for j in range(5):
        acc = [Z([0, 0, 0, 0]) for _ in range(M)]
        for a in range(5):
            p = smul(sneg(Rhat[a][i]), Rhat[(3 - a) % 5][j])
            acc = [x + y for x, y in zip(acc, p)]
        acc = [x * Fr(1, 5) for x in acc]
        target = [Z([1 if (i == j and n == 0) else 0, 0, 0, 0]) for n in range(M)]
        if not all(x == y for x, y in zip(acc, target)):
            ok = False
            print("FAIL", i, j, acc)
print("Rhat(-t)^T eta^{-1} Rhat(t) = 1 mod t^%d:" % M, ok)

# M_0 = Psi_0^T eta^{-1} Rhat ; entries (M0)_{ji} = (1/5) sum_b eps_{3-b} zeta^{j(4-b)} Rhat_b^{(i)}
M0 = [[None] * 5 for _ in range(5)]
for j in range(5):
    for i in range(5):
        acc = [Z([0, 0, 0, 0]) for _ in range(M)]
        for b in range(5):
            f = Z.zeta(j * (4 - b)) * Fr(eps[(3 - b) % 5], 5)
            acc = [x + f * y for x, y in zip(acc, Rhat[b][i])]
        M0[j][i] = acc
print("M0(0) = 1:", all(M0[j][i][0] == (1 if i == j else 0) for i in range(5) for j in range(5)))
# M0(-t)^T M0(t) = 1
ok = True
for i in range(5):
    for j in range(5):
        acc = [Z([0, 0, 0, 0]) for _ in range(M)]
        for k in range(5):
            acc = [x + y for x, y in zip(acc, smul(sneg(M0[k][i]), M0[k][j]))]
        if not all(acc[n] == (1 if (i == j and n == 0) else 0) for n in range(M)):
            ok = False
print("M0(-t)^T M0(t) = 1 mod t^%d:" % M, ok)
print("M0 first-order coefficient matrix (rows j, cols i):")
for j in range(5):
    print([M0[j][i][1] for i in range(5)])
# cyclic covariance: (M0)_{j+1,i+1}(t) = (M0)_{ji}(t)? and Galois
cyc = all(M0[(j + 1) % 5][(i + 1) % 5][n] == M0[j][i][n] for i in range(5) for j in range(5) for n in range(M))
print("M0 commutes with the cyclic permutation of colours:", cyc)
# unit: M0^{-1}(z) applied to 1' = sum_i zeta^i e_i : translation t(z) = z(1' - M0^{-1}(z) 1'); M0^{-1}(t) = M0(-t)^T
unit = [Z.zeta(i) for i in range(5)]
tr = []
for j in range(5):
    acc = [Z([0, 0, 0, 0]) for _ in range(M)]
    for i in range(5):
        # (M0(-t)^T)_{ji} = M0_{ij}(-t)
        acc = [x + unit[i] * y for x, y in zip(acc, sneg(M0[i][j]))]
    tr.append(acc)
print("M0^{-1}(t) 1' , coefficients of t^0, t^1, t^2 (components j=0..4):")
for n in range(3):
    print(n, [tr[j][n] for j in range(5)])
