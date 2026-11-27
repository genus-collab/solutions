# -*- coding: utf-8 -*-
"""ref_kpoly.py -- univariate polynomials over a field ring K (KRing / IRing of ref_sys): lists of K-elements, low -> high."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sys import *


def ktrim(K, a):
    a = list(a)
    while a and K.iszero(a[-1]): a.pop()
    return a


def kdivmod(K, a, b):
    a = ktrim(K, a); b = ktrim(K, b)
    if not b: raise ZeroDivisionError
    if len(a) < len(b): return [], a
    inv = K.inv(b[-1]); q = [K.zero] * (len(a) - len(b) + 1); db = len(b) - 1
    for s in range(len(a) - len(b), -1, -1):
        f = K.mul(a[s + db], inv)
        if not K.iszero(f):
            q[s] = f
            for i in range(db + 1): a[s + i] = K.sub(a[s + i], K.mul(f, b[i]))
    return ktrim(K, q), ktrim(K, a[:db])


def kmonic(K, a):
    a = ktrim(K, a)
    if not a: return a
    inv = K.inv(a[-1]); return [K.mul(inv, x) for x in a]


def kgcd(K, a, b):
    a = ktrim(K, a); b = ktrim(K, b)
    while b: a, b = b, kdivmod(K, a, b)[1]
    return kmonic(K, a)


def kgcdl(K, fs):
    g = []
    for f in fs:
        f = ktrim(K, f)
        if f: g = kgcd(K, g, f) if g else kmonic(K, f)
    return g


def keval(K, a, x):
    r = K.zero
    for c in reversed(a): r = K.add(K.mul(r, x), c)
    return r


def keval_int(K, a, x):
    """evaluate at an integer x (cheap: scalar multiplications only)"""
    r = K.zero
    for c in reversed(a): r = K.add(K.scal(x, r), c)
    return r


def kpow(K, x, n):
    r = K.one
    for _ in range(n): r = K.mul(r, x)
    return r


def kres_true(K, f, g):
    """resultant at the true degrees (Euclid), Sylvester sign convention"""
    f = ktrim(K, f); g = ktrim(K, g)
    if not f or not g: return K.zero
    res = K.one
    while True:
        m = len(f) - 1; n = len(g) - 1
        if n == 0: return K.mul(res, kpow(K, g[0], m))
        if m == 0: return K.mul(res, kpow(K, f[0], n))
        if m < n:
            f, g = g, f
            if (m * n) & 1: res = K.sub(K.zero, res)
            continue
        r = kdivmod(K, f, g)[1]
        if not r: return K.zero
        k = len(r) - 1; res = K.mul(res, kpow(K, g[-1], m - k))
        if (m * n) & 1: res = K.sub(K.zero, res)
        f, g = g, r


def kres(K, f, g, df, dg):
    f = ktrim(K, f); g = ktrim(K, g); m = len(f) - 1; n = len(g) - 1
    assert m <= df and n <= dg
    if (m < df and n < dg) or not f or not g: return K.zero
    if m == df: return K.mul(kpow(K, f[-1], dg - n), kres_true(K, f, g))
    r = K.mul(kpow(K, g[-1], df - m), kres_true(K, f, g))
    return K.sub(K.zero, r) if (dg * (df - m)) & 1 else r


def knewton(K, x0, ys, p):
    """Newton interpolation over K on the integer nodes x0, x0 + 1, ..."""
    n = len(ys); c = list(ys); inv = [0] + [pow(j, p - 2, p) for j in range(1, n)]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1): c[i] = K.scal(inv[j], K.sub(c[i], c[i - 1]))
    out = [c[n - 1]]
    for i in range(n - 2, -1, -1):
        xi = x0 + i; new = [K.zero] + list(out)
        for k, v in enumerate(out): new[k] = K.sub(new[k], K.scal(xi, v))
        new[0] = K.add(new[0], c[i]); out = new
    return ktrim(K, out)


def kmult_root(K, a, r):
    """multiplicity of the root r (K-element) in a"""
    m = 0; a = ktrim(K, a)
    while a and K.iszero(keval(K, a, r)):
        a = kdivmod(K, a, [K.sub(K.zero, r), K.one])[0]; m += 1
    return m, a
