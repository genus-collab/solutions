# -*- coding: utf-8 -*-
"""g5l_lib.py -- shared helpers for the LOCAL audit of the genus-5 cluster run (hodge_vertex_quintic/natural_form/genus5_cluster/local).

Reads only: the package's parts/<tag>/*.json (per-inner-graph pieces) and results/<tag>.json (merged tables), both of
cluster_genus5_quintic/.  Writes nothing.  The package itself is imported for g5_checks (tags, read_parts, merge) and for
fq_fp (the two Mersenne primes); nothing in the package is modified.

Conventions (see g5_engine.py / fq4_split.py):
  inner graph t  = a leg-free stable graph of genus 5 (4555 of them, the index is the position in data/graphs_5_0_12.json);
  piece (t -> q) = the sum of the terms of t whose OUTER edge set S contracts t to the outer graph q; |S| = E(q);
  with --maxout m only the subsets with |S| <= m are kept, so exactly the outer graphs q with E(q) <= m get pieces,
  and those pieces are then COMPLETE (every S of size <= m is included).
  A row of the table = sum over all inner graphs t of the piece (t -> q); K = 8 residues of u^0..u^7 mod p.
"""
import os, sys, json

PKG = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))), "cluster_genus5_quintic")
assert os.path.isdir(PKG), PKG
if PKG not in sys.path:
    sys.path.insert(0, PKG)

K = 8           # genus 5: coefficients u^0 .. u^(2g-3)
NIN = 4555      # leg-free stable graphs of genus 5


def prime(name="89"):
    from fq_fp import P1, P2
    return P2 if name == "89" else P1


def result_path(tag):
    return os.path.join(PKG, "results", "%s.json" % tag)


def load_result(tag):
    """the merged table of a collect step: dict with 'table' {str(q): [K residues]}, 'info', 'failed', 'outer'"""
    with open(result_path(tag)) as f:
        return json.load(f)


def load_parts(tag):
    """(names, pieces, times) straight from the package's own reader; pieces: {str(t): {str(q): [K residues]}}"""
    import g5_checks as C
    return C.read_parts(tag)


def merge(pieces, p):
    """sum the pieces into rows {q: [K residues]} (the package's own merge)"""
    import g5_checks as C
    return C.merge(pieces, K, p)


def outer_meta(tag="g5_89_nat"):
    """per outer graph: dict(E, gv, ed, aut); index = outer index = inner index (same graph list)"""
    return load_result(tag)["outer"]


def row(table, q):
    """a row of a table, as a list of K ints ([0]*K if absent).  Accepts both key styles: the JSON tables of results/ are keyed
    by the string of the outer index, the tables returned by g5_checks.merge by the integer."""
    v = table.get(str(q))
    if v is None:
        v = table.get(int(q))
    return list(v) if v else [0] * K


def piece_sum(res, p):
    """sum of all pieces of one inner graph (res = {str(q): [K residues]}) -> [K residues]"""
    s = [0] * K
    for c in res.values():
        for k in range(K):
            s[k] = (s[k] + c[k]) % p
    return s


def fmt_rat(x, p, margin=10 ** 6):
    """rational reconstruction of a residue for printing only (None if it does not fit); evidence is the residue itself.
    Extended Euclid on (p, x) keeping r_i = t_i * x mod p, so x = r_i / t_i whenever both are small; verified mod p.
    Accepted only when |num| * |den| <= p / margin, i.e. with `margin` decimal digits of slack over the uniqueness bound
    |num| * |den| <= p / 2, so that a residue which one prime cannot determine prints as None instead of as a spurious
    14-digit-over-14-digit fraction."""
    if x % p == 0:
        return "0"
    from fractions import Fraction as Fr
    nmax = p // margin
    a, b = p, x % p
    ta, tb = 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        ta, tb = tb, ta - q * tb
        if ta and abs(a) * abs(ta) <= nmax:
            f = Fr(a, ta)
            if f.denominator % p and (f.numerator - f.denominator * x) % p == 0:
                return str(f)
    return None


def log_opener(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    f = open(path, "w", encoding="utf-8")

    def log(*x):
        s = " ".join(str(y) for y in x)
        print(s, flush=True)
        f.write(s + "\n")
        f.flush()
    return log, f
