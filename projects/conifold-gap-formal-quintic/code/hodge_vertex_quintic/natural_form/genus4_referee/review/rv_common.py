# -*- coding: utf-8 -*-
"""rv_common.py -- adversarial review of fq4_split.py / genus4_split.py: the mutant factory.

Nothing in the project folder is written: sys.dont_write_bytecode is set before any project import, the project modules are imported
read-only from PROJ, and every mutant is a COPY of fq4_split.py written to review/mutants/fq4_split_<name>.py with one textual
replacement (asserted to match exactly once in the original source).  Every copy (also 'baseline') carries one extra non-mutating
patch: a module global POPMAX (default 99) that restricts the outer-edge subsets to popcount <= POPMAX, used for the cheap genus-4
runs restricted to the outer graphs with few edges (with POPMAX = 99 the code is the original one).
MUTANTS: name -> (list of (old, new) replacements, one-line description).
"""
import os, sys, importlib.util
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
MUT = os.path.join(HERE, "mutants"); OUT = os.path.join(HERE, "out")
if PROJ not in sys.path: sys.path.insert(0, PROJ)

POP_OLD = 'self.outer_pos[self.contract(G, S)]) for S in range(1 << E) if '
POP_NEW = 'self.outer_pos[self.contract(G, S)]) for S in range(1 << E) if bin(S).count("1") <= POPMAX and '
IMP_OLD = 'import os, sys, itertools, math\n'
IMP_NEW = 'import os, sys, itertools, math\nPOPMAX = 99\n'

EM_OLD = 'self.scale(tp[a + b + 1], Eh[(a, b)][i][j]) if Eh[(a, b)][i][j] else None'
INV_OLD = 'inv = [Fr(4, 5 * k + 4) for k in range(n)]'

MUTANTS = {
    'baseline': ([], 'exact copy (control: every check must pass)'),
    'm01_EM_sign': ([(EM_OLD, 'self.scale(tp[a + b + 1], (-Eh[(a, b)][i][j]) % p) if Eh[(a, b)][i][j] else None')],
                    'E_M -> -E_M'),
    'm02_tau1': ([(INV_OLD, 'inv = [Fr(4, 5 * k + (5 if k == 1 else 4)) for k in range(n)]')],
                 '1/tau: coefficient of u^1 4/9 -> 4/10 (tau_1 and all later tau_k change)'),
    'm03_tau3': ([(INV_OLD, 'inv = [Fr(4, 5 * k + (5 if k == 3 else 4)) for k in range(n)]')],
                 '1/tau: coefficient of u^3 4/19 -> 4/20 (tau_3, tau_4, tau_5 change)'),
    'm04_tau4': ([(INV_OLD, 'inv = [Fr(4, 5 * k + (5 if k == 4 else 4)) for k in range(n)]')],
                 '1/tau: coefficient of u^4 4/24 -> 4/25 (tau_4, tau_5 change; invisible at genus 3 through u^3)'),
    'm05_tau5': ([(INV_OLD, 'inv = [Fr(4, 5 * k + (5 if k == 5 else 4)) for k in range(n)]')],
                 '1/tau: coefficient of u^5 4/29 -> 4/30 (only tau_5 changes)'),
    'm06_taupow': ([(EM_OLD, 'self.scale(tp[a + b], Eh[(a, b)][i][j]) if Eh[(a, b)][i][j] else None')],
                   'E_M(a,b) = tau^(a+b) Eh_0 instead of tau^(a+b+1)'),
    'm07_taupow_hi': ([(EM_OLD, 'self.scale(tp[a + b + 1 if a + b + 1 <= 6 else a + b], Eh[(a, b)][i][j]) if Eh[(a, b)][i][j] else None')],
                      'wrong tau power only for zhat-orders a+b+1 >= 7 (these orders do not exist at genus 3)'),
    'm08_EM_transpose': ([(EM_OLD, 'self.scale(tp[a + b + 1], Eh[(a, b)][j][i]) if Eh[(a, b)][j][i] else None')],
                         'E_M(a,b)[i][j] -> E_M(a,b)[j][i]'),
    'm09_no_mult': ([('c0 = R.inv(math.prod(math.factorial(x) for x in mult.values()))', 'c0 = 1')],
                    'drop the 1/prod(mult!) symmetry factor of equal leaves'),
    'm10_contract_loop': ([('else: ne[cid[find(u)]] += 1', 'else: ne[cid[find(u)]] += (0 if u == w else 1)')],
                          'contraction: an inner loop does not raise the genus of its component'),
    'm11_contract_inverted': ([('if not (S >> e) & 1: par[find(u)] = find(w)', 'if (S >> e) & 1: par[find(u)] = find(w)'),
                               ('if (S >> e) & 1: out.append((cid[find(u)], cid[find(w)]))', 'if not (S >> e) & 1: out.append((cid[find(u)], cid[find(w)]))')],
                              'contraction contracts the OUTER edges instead of the inner ones (keys stay valid)'),
    'm12_contract_gswap': ([('gv = [gv[c] + ne[c] - nv[c] + 1 for c in range(len(roots))]',
                             'gv = [gv[c] + ne[c] - nv[c] + 1 for c in range(len(roots))][::-1]')],
                           'contraction: the component genera are attached to the components in reversed order'),
    'm13_flip': ([('mats[e][c][key[pos]] if flip else mats[e][key[pos]][c]', 'mats[e][key[pos]][c] if flip else mats[e][c][key[pos]]')],
                 'colour_sum: edge matrix read transposed'),
    'm14_phase': ([('pow(self.zph[j], n, p)', 'pow(self.zph[j], n + 1, p)')], 'vertex phase zeta^(4j(n+1)) instead of zeta^(4jn)'),
    'm15_xpow': ([('self.scale(self.xpow[n], ', 'self.scale(self.xpow[n - 1], ')], 'xhat^(n-1) instead of xhat^n at every vertex'),
    'm16_popcount': ([('if bin(S).count("1") < K]', 'if bin(S).count("1") < K - 1]')], 'subset restriction popcount < K-1'),
    'm17_valprune': ([('if (S >> e) & 1) >= K: continue', 'if (S >> e) & 1) >= K - 1: continue')], 'valuation pruning one order too early'),
    'm18_leafmax': ([('partitions_min2(m, dim, self.KZ + 1)', 'partitions_min2(m, dim, self.KZ)')], 'largest leaf exponent 3g-2 dropped'),
    'm19_perms': ([('for perm in set(itertools.permutations(hp)): lst.append((perm, vec))', 'for perm in [hp]: lst.append((perm, vec))')],
                  'only the sorted assignment of psi-exponents to the half-edges of a vertex'),
    'm20_keys_swap': ([('keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]', 'keys = [(hp[2 * e + 1], hp[2 * e]) for e in range(E)]')],
                      'edge exponents (a, b) swapped relative to the colours'),
    'm21_no_aut': ([("ia = self.R.inv(G['aut'])", 'ia = 1')], 'drop 1/|Aut|'),
    'm22_EO_sign': ([('self.EO = {key: [[self.sub(self.ER[key][i][j], self.EM[key][i][j])', 'self.EO = {key: [[self.add(self.ER[key][i][j], self.EM[key][i][j])')],
                    'E_out = E_R + E_M'),
    'm23_Minv_sign': ([('(1 if k % 2 == 0 else p - 1) % p for j in range(5)]', '1 % p for j in range(5)]')],
                      'M_0^-1_k = C_k^T without the sign (-1)^k'),
    'm24_tau_reciprocal': ([('tau.append(-sum(inv[j] * tau[k - j] for j in range(1, k + 1)))\n    return tau', 'tau.append(-sum(inv[j] * tau[k - j] for j in range(1, k + 1)))\n    return inv')],
                           'tau = 2F1 instead of 1/2F1'),
    'm25_loopmat': ([('x = mul(x, mats[e][c][c])', 'x = mul(x, mats[e][c][(c + 1) % 5])')], 'loop edge contracted off the diagonal'),
    'm26_leaf_from3': ([('for x in range(min(cap, left - 2 * (k - 1)), 1, -1):', 'for x in range(min(cap, left - 2 * (k - 1)), 2, -1):')],
                       'leaves start at exponent 3 (the psi^2 leaf dropped)'),
}


def mutant_source(name):
    src = open(os.path.join(PROJ, "fq4_split.py"), encoding="utf-8").read()
    for old, new in MUTANTS[name][0]:
        assert src.count(old) == 1, ("replacement must match exactly once", name, old, src.count(old))
        src = src.replace(old, new)
    for old, new in ((POP_OLD, POP_NEW), (IMP_OLD, IMP_NEW)):
        assert src.count(old) == 1, ("patch must match exactly once", name, old)
        src = src.replace(old, new)
    return src


def write_mutants(names=None):
    """write the copies (main process only); returns {name: path}"""
    os.makedirs(MUT, exist_ok=True); os.makedirs(OUT, exist_ok=True); paths = {}
    for name in (names or sorted(MUTANTS)):
        path = os.path.join(MUT, "fq4_split_%s.py" % name); src = mutant_source(name)
        if not (os.path.exists(path) and open(path, encoding="utf-8").read() == src):
            with open(path, "w", encoding="utf-8", newline="\n") as f: f.write(src)
        paths[name] = path
    return paths


_MODS = {}


def mutant_module(name):
    mod = _MODS.get(name)
    if mod is None:
        path = os.path.join(MUT, "fq4_split_%s.py" % name)
        spec = importlib.util.spec_from_file_location("fq4_split_%s" % name, path)
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); _MODS[name] = mod
    return mod
