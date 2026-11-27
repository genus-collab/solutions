# -*- coding: utf-8 -*-
"""tr_field.py -- coefficient fields for the TR computations of the allgenus/tr task (new file, stdlib only).

Two fields with the same operator interface (+, -, *, /, ==, unary -, .iszero(), .zeta(k), .conj(m)):
  Qz   : exact Q(zeta), zeta a primitive fifth root of unity; element a0 + a1 z + a2 z^2 + a3 z^3 with Fractions.
  Fp(p): integers mod a prime p = 1 mod 5 with a fixed primitive fifth root of unity zeta (p = 2^89-1, 2^61-1 are used).
Mixed arithmetic with int and Fraction is supported (they are coerced).  F.one, F.zero, F.of(x), F.z(k) = zeta^k.
Rational reconstruction: Fp.ratrec(a) (lattice reduction, returns a Fraction or None).
Self-test: python tr_field.py
"""
from fractions import Fraction as Fr

P61 = (1 << 61) - 1
P89 = (1 << 89) - 1


# ------------------------------------------------------------------ Q(zeta5)
class Qz(object):
    __slots__ = ('c',)

    def __init__(self, c=(0, 0, 0, 0)):
        self.c = tuple(Fr(x) for x in c)

    @staticmethod
    def of(x):
        if isinstance(x, Qz): return x
        return Qz((x, 0, 0, 0))

    @staticmethod
    def z(k):
        k %= 5
        if k == 4: return Qz((-1, -1, -1, -1))
        c = [0, 0, 0, 0]; c[k] = 1
        return Qz(c)

    def _co(self, o):
        if isinstance(o, Qz): return o
        if isinstance(o, (int, Fr)): return Qz((o, 0, 0, 0))
        return NotImplemented

    def __add__(self, o):
        if isinstance(o, (int, Fr)): return Qz((self.c[0] + o, self.c[1], self.c[2], self.c[3]))
        o = self._co(o)
        return Qz(tuple(a + b for a, b in zip(self.c, o.c)))
    __radd__ = __add__

    def __neg__(self): return Qz(tuple(-a for a in self.c))

    def __sub__(self, o):
        if isinstance(o, (int, Fr)): return Qz((self.c[0] - o, self.c[1], self.c[2], self.c[3]))
        o = self._co(o)
        return Qz(tuple(a - b for a, b in zip(self.c, o.c)))

    def __rsub__(self, o): return (-self) + o

    def __mul__(self, o):
        if isinstance(o, (int, Fr)):
            if o == 0: return QZ0
            return Qz(tuple(a * o for a in self.c))
        a = self.c; b = o.c
        if b[1] == 0 and b[2] == 0 and b[3] == 0:
            r = b[0]
            return Qz((a[0] * r, a[1] * r, a[2] * r, a[3] * r)) if r else QZ0
        if a[1] == 0 and a[2] == 0 and a[3] == 0:
            r = a[0]
            return Qz((b[0] * r, b[1] * r, b[2] * r, b[3] * r)) if r else QZ0
        c = [0] * 7
        for i in range(4):
            ai = a[i]
            if ai == 0: continue
            for j in range(4):
                if b[j] != 0: c[i + j] += ai * b[j]
        c[0] += c[5]; c[1] += c[6]; r = c[4]
        return Qz((c[0] - r, c[1] - r, c[2] - r, c[3] - r))
    __rmul__ = __mul__

    def conj(self, m):
        """Galois sigma_m: zeta -> zeta^m (m = 1..4)"""
        v = [0] * 5
        for k in range(4): v[(m * k) % 5] += self.c[k]
        r = v[4]
        return Qz((v[0] - r, v[1] - r, v[2] - r, v[3] - r))

    def norm(self):
        n = self * self.conj(2) * self.conj(3) * self.conj(4)
        assert n.c[1] == n.c[2] == n.c[3] == 0
        return n.c[0]

    def inv(self):
        num = self.conj(2) * self.conj(3) * self.conj(4)
        n = (self * num).c
        assert n[1] == n[2] == n[3] == 0 and n[0] != 0, "inverse of 0 in Q(zeta)"
        return num * (1 / n[0])

    def __truediv__(self, o):
        if isinstance(o, (int, Fr)): return Qz(tuple(a / o for a in self.c))
        return self * o.inv()

    def __rtruediv__(self, o): return self.inv() * o

    def __eq__(self, o):
        if isinstance(o, (int, Fr)): o = Qz((o, 0, 0, 0))
        return self.c == o.c

    def __ne__(self, o): return not self.__eq__(o)
    def __hash__(self): return hash(self.c)
    def iszero(self): return not any(self.c)
    def is_rational(self): return self.c[1] == 0 and self.c[2] == 0 and self.c[3] == 0
    def rational(self):
        assert self.is_rational(), self
        return self.c[0]

    def __repr__(self):
        terms = []
        for k, a in enumerate(self.c):
            if a: terms.append(str(a) + ("" if k == 0 else (" z" if k == 1 else " z^%d" % k)))
        return "(" + (" + ".join(terms) if terms else "0") + ")"


QZ0 = Qz()


class QzField(object):
    """field object for Qz (uniform interface with FpField)"""
    name = "Q(zeta5)"
    exact = True

    def __init__(self):
        self.zero = Qz(); self.one = Qz((1, 0, 0, 0))
        self._z = [Qz.z(k) for k in range(5)]

    def of(self, x): return Qz.of(x)
    def z(self, k): return self._z[k % 5]
    def fmt(self, a): return repr(a)


# ------------------------------------------------------------------ F_p
def make_fp(p, zeta=None, order=5):
    """returns (FpField object); elements are instances of a class bound to p; zeta: optional fixed fifth root of unity"""
    assert p % order == 1

    class Fp(object):
        __slots__ = ('v',)
        P = p

        def __init__(self, v=0): self.v = v % p

        @staticmethod
        def of(x):
            if isinstance(x, Fp): return x
            if isinstance(x, Fr): return Fp(x.numerator * pow(x.denominator, p - 2, p))
            return Fp(x)

        def __add__(self, o):
            if isinstance(o, Fp): return Fp(self.v + o.v)
            return Fp(self.v + Fp.of(o).v)
        __radd__ = __add__

        def __sub__(self, o):
            if isinstance(o, Fp): return Fp(self.v - o.v)
            return Fp(self.v - Fp.of(o).v)

        def __rsub__(self, o): return Fp(Fp.of(o).v - self.v)
        def __neg__(self): return Fp(-self.v)

        def __mul__(self, o):
            if isinstance(o, Fp): return Fp(self.v * o.v)
            return Fp(self.v * Fp.of(o).v)
        __rmul__ = __mul__

        def inv(self):
            assert self.v, "inverse of 0 mod p"
            return Fp(pow(self.v, p - 2, p))

        def __truediv__(self, o): return self * Fp.of(o).inv()
        def __rtruediv__(self, o): return Fp.of(o) * self.inv()

        def __eq__(self, o):
            if not isinstance(o, Fp): o = Fp.of(o)
            return self.v == o.v

        def __ne__(self, o): return not self.__eq__(o)
        def __hash__(self): return hash(self.v)
        def iszero(self): return self.v == 0
        def __repr__(self): return "Fp(%d)" % self.v

    class FpField(object):
        name = "F_p, p = %d" % p
        exact = False

        def __init__(self):
            self.p = p; self.cls = Fp
            self.zero = Fp(0); self.one = Fp(1)
            self.order = order
            if zeta is not None:
                zz = zeta % p; assert pow(zz, order, p) == 1 and zz != 1
            else:
                g = 2
                while True:
                    zz = pow(g, (p - 1) // order, p)
                    if all(pow(zz, order // q, p) != 1 for q in range(2, order + 1) if order % q == 0 and all(q % r for r in range(2, q))):
                        break
                    g += 1
            self._z = [Fp(pow(zz, k, p)) for k in range(order)]

        def of(self, x): return Fp.of(x)
        def z(self, k): return self._z[k % self.order]

        def ratrec(self, a, bound=None):
            """rational reconstruction of a (Fp or int): Fraction r/s with |r|, s <= sqrt(p/2), or None"""
            a = a.v if isinstance(a, Fp) else a % p
            B = bound or int((p // 2) ** 0.5)
            r0, r1, s0, s1 = p, a, 0, 1
            while r1 > B:
                q = r0 // r1
                r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
            if s1 == 0 or abs(s1) > B: return None
            f = Fr(r1, s1)
            return f if (f.numerator - a * f.denominator) % p == 0 else None

        def fmt(self, a):
            f = self.ratrec(a)
            return str(f) if f is not None else "~%d" % (a.v if isinstance(a, Fp) else a)

    return FpField()


if __name__ == "__main__":
    import random
    rng = random.Random(1)
    Q = QzField()
    def rq(): return Qz(tuple(Fr(rng.randrange(-9, 10), rng.randrange(1, 7)) for _ in range(4)))
    z = Q.z(1)
    assert z * z * z * z * z == Q.one and not (z - Q.one).iszero()
    assert sum((Q.z(k) for k in range(5)), Q.zero).iszero()
    for _ in range(50):
        a, b, c = rq(), rq(), rq()
        assert (a * b) * c == a * (b * c) and a * (b + c) == a * b + a * c
        if not a.iszero(): assert a * a.inv() == Q.one
        assert (a * b).conj(2) == a.conj(2) * b.conj(2)
    assert (1 - z).norm() == 5
    for p in (P89, P61):
        F = make_fp(p); zz = F.z(1)
        assert (zz * zz * zz * zz * zz - F.one).iszero() and not (zz - F.one).iszero()
        x = F.of(Fr(-7, 12)); assert F.ratrec(x) == Fr(-7, 12)
        assert F.of(Fr(3, 5)) * 5 == F.of(3)
    print("tr_field ok: Q(zeta5) ring axioms, Galois, norm(1-zeta) = 5; F_p zeta and rational reconstruction at 2^89-1, 2^61-1")
