"""Fresh implementation of Definition 1.3 (draft s.1.3) using only fractions.

Representation: a "row" is a list over z-order n of Laurent polynomials in v:
    dict {exponent: Fraction}.  All arithmetic exact.
"""
from fractions import Fraction as Fr


def rec_r(n, bad=None):
    """(Rec): r_{-1}=0, r_0=1, 27(m+1)r_{m+1}=3(3m^2+3m+1)r_m-m^3 r_{m-1}."""
    r = [Fr(1)]
    rm1 = Fr(0)
    for m in range(0, n):
        nxt = (3 * (3 * m * m + 3 * m + 1) * r[m] - m ** 3 * rm1) / (27 * (m + 1))
        rm1 = r[m]
        r.append(nxt)
    if bad is not None:
        # control: perturb one r_l
        l, val = bad
        r[l] = r[l] + val
    return r


def lp_add(a, b, c=Fr(1)):
    out = dict(a)
    for k, x in b.items():
        out[k] = out.get(k, Fr(0)) + c * x
        if out[k] == 0:
            del out[k]
    return out


def lp_mul(a, b):
    out = {}
    for k, x in a.items():
        for l, y in b.items():
            out[k + l] = out.get(k + l, Fr(0)) + x * y
    return {k: x for k, x in out.items() if x != 0}


def lp_scale(a, c):
    return {k: x * c for k, x in a.items() if x * c != 0} if c != 0 else {}


def lp_shift(a, s):
    return {k + s: x for k, x in a.items()}


def lp_dv(a):
    return {k - 1: k * x for k, x in a.items() if k != 0}


def lp_eval(a, v):
    # v may be Fraction or anything supporting ** with ints and *,+
    tot = 0
    for k, x in a.items():
        tot = tot + (v ** k) * x
    return tot


def series_S(E, N, r=None, sign_denom=1):
    """S(v,z) as list [z^n] -> dict in v, n=0..N.
    S = Sum_l r_l (vz)^l / prod_{j=1}^{l+1}(1 + j z/E)
    Control: sign_denom=-1 replaces 1+jz/E by 1-jz/E.
    """
    if r is None:
        r = rec_r(N + 1)
    out = [{} for _ in range(N + 1)]
    for l in range(0, N + 1):
        # 1/prod(1+jz/E) series up to z^(N-l)
        ser = [Fr(0)] * (N - l + 1)
        ser[0] = Fr(1)
        for j in range(1, l + 2):
            c = sign_denom * Fr(j) / E
            # divide by (1+c z)
            for n in range(1, N - l + 1):
                ser[n] = ser[n] - c * ser[n - 1]
        for n in range(0, N - l + 1):
            if ser[n] != 0 and r[l] != 0:
                out[n + l][l] = out[n + l].get(l, Fr(0)) + r[l] * ser[n]
    return out


def z_mul_lp(row, lp):
    return [lp_mul(x, lp) for x in row]


def theta_v(row, E):
    """theta acting on a function of v: theta v = 1/E - v/3."""
    tv = {0: Fr(1) / E, 1: Fr(-1, 3)}
    return [lp_mul(lp_dv(x), tv) for x in row]


def z_shift(row, N):
    """multiply by z, truncate at z^N"""
    return [{}] + row[:N]


def row_add(a, b, c=Fr(1)):
    return [lp_add(x, y, c) for x, y in zip(a, b)]


def Th_plain(f, E, N):
    """Th = Y + z theta on f, Y=-E/3+1/v."""
    Y = {0: -E / 3, -1: Fr(1)}
    return row_add(z_mul_lp(f, Y), z_shift(theta_v(f, E), N))


def z_flip(row):
    """z -> -z"""
    return [lp_scale(x, Fr((-1) ** n)) for n, x in enumerate(row)]


def Th_flip(f, E, N):
    """Y - z theta on f (row at -z)."""
    Y = {0: -E / 3, -1: Fr(1)}
    return row_add(z_mul_lp(f, Y), z_shift(theta_v(f, E), N), Fr(-1))


def abc_vs(E, N, r=None, sign_denom=1):
    """Returns a,b as v-Laurent rows at +z, and the pieces of c:
    c = Delta_c * Th b, Th b kept as v-Laurent row (Thb)."""
    S = series_S(E, N, r, sign_denom)
    a = [lp_mul(x, {1: Fr(1)}) for x in S]
    b = Th_plain(a, E, N)
    Thb = Th_plain(b, E, N)
    return a, b, Thb


def abc_flip(E, N, r=None, sign_denom=1):
    """a,b,Thb at -z, each expressed as a power series in z"""
    S = series_S(E, N, r, sign_denom)
    a = z_flip([lp_mul(x, {1: Fr(1)}) for x in S])
    b = Th_flip(a, E, N)
    Thb = Th_flip(b, E, N)
    return a, b, Thb


class Qw:
    """Q(omega), omega^2=-1-omega."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fr(a)
        self.b = Fr(b)

    def _c(o):
        return o if isinstance(o, Qw) else Qw(o)

    def __add__(self, o):
        o = Qw._c(o)
        return Qw(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Qw(-self.a, -self.b)

    def __sub__(self, o):
        o = Qw._c(o)
        return Qw(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return Qw._c(o) - self

    def __mul__(self, o):
        o = Qw._c(o)
        # (a+bw)(c+dw)=ac+(ad+bc)w+bd w^2 ; w^2=-1-w
        ac = self.a * o.a
        bd = self.b * o.b
        return Qw(ac - bd, self.a * o.b + self.b * o.a - bd)

    __rmul__ = __mul__

    def inv(self):
        # norm = a^2 - ab + b^2 ; conj(w)=w^2=-1-w
        n = self.a ** 2 - self.a * self.b + self.b ** 2
        # conj = a + b w^2 = (a-b) - b w
        return Qw((self.a - self.b) / n, -self.b / n)

    def __truediv__(self, o):
        return self * Qw._c(o).inv()

    def __rtruediv__(self, o):
        return Qw._c(o) * self.inv()

    def __pow__(self, n):
        if n == 0:
            return Qw(1)
        if n < 0:
            return self.inv() ** (-n)
        out = Qw(1)
        for _ in range(n):
            out = out * self
        return out

    def __eq__(self, o):
        o = Qw._c(o)
        return self.a == o.a and self.b == o.b

    def __ne__(self, o):
        return not self == o

    def __repr__(self):
        return f"({self.a}+{self.b}w)"

    def iszero(self):
        return self.a == 0 and self.b == 0


OMEGA = Qw(0, 1)
