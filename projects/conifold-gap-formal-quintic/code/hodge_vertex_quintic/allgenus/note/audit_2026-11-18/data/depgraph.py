import re, sys, json, os
NOTE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
BS = chr(92)

def load(path, out):
    full = os.path.join(NOTE, path + ('' if path.endswith('.tex') else '.tex'))
    rel = os.path.relpath(full, NOTE).replace(BS, '/')
    with open(full, encoding='utf-8') as f:
        lines = f.read().split('\n')
    for i, l in enumerate(lines, 1):
        m = re.match(r'\s*\\input\{([^}]*)\}', l)
        if m and m.group(1) not in ('preamble', 'macros'):
            load(m.group(1), out)
        else:
            out.append((rel, i, l))

out = []
load('main', out)
ENVS = ['thm', 'lem', 'prop', 'cor', 'defn', 'rmk', 'conv', 'notation', 'example']
ENVRE = '|'.join(ENVS)
items = []
cur = None
state = None

def refs(s):
    r = []
    for m in re.finditer(r'\\(?:[cC]ref|ref|eqref|[cC]refrange)\{([^}]*)\}', s):
        r += [x.strip() for x in m.group(1).split(',')]
    return r

last_item = None
proof_of = None
eq_owner = {}
for (f, i, l) in out:
    ls = '' if l.strip().startswith('%') else l
    m = re.search(r'\\begin\{(' + ENVRE + r')\}', ls)
    if m:
        cur = dict(label=None, env=m.group(1), file=f, line=i, stmt_refs=[], proof_refs=[], labels=[])
        items.append(cur)
        state = 'stmt'
    labs = re.findall(r'\\label\{([^}]*)\}', ls)
    if state == 'stmt':
        for lab in labs:
            cur['labels'].append(lab)
            if cur['label'] is None:
                cur['label'] = lab
        cur['stmt_refs'] += refs(ls)
        if re.search(r'\\end\{(' + ENVRE + r')\}', ls):
            state = None
            last_item = cur
    m = re.search(r'\\begin\{proof\}(\[[^\]]*\])?', ls)
    if m:
        state = 'proof'
        tgt = last_item
        if m.group(1):
            rr = refs(m.group(1))
            if rr:
                for it in items:
                    if it['label'] == rr[0]:
                        tgt = it
        proof_of = tgt
        if proof_of is not None:
            proof_of['has_proof'] = True
            proof_of.setdefault('proof_at', '%s:%d' % (f, i))
    if state == 'proof':
        if proof_of is not None:
            proof_of['proof_refs'] += refs(ls)
            for lab in labs:
                eq_owner[lab] = proof_of['label']
        if re.search(r'\\end\{proof\}', ls):
            state = None

json.dump(items, open('items.json', 'w'), indent=1)
by = {it['label']: it for it in items if it['label']}
for it in items:
    for lab in it['labels']:
        eq_owner.setdefault(lab, it['label'])
print(len(items), 'items;', sum(1 for it in items if it.get('has_proof')), 'with proof env')
for it in items:
    if it['env'] in ('thm', 'lem', 'prop', 'cor') and not it.get('has_proof'):
        print('NO PROOF ENV:', it['label'], it['file'], it['line'])

def res(y):
    if y in by:
        return y
    return eq_owner.get(y)

def closure(root):
    seen = []
    stack = [root]
    while stack:
        x = stack.pop()
        if x is None or x in seen or x not in by:
            continue
        seen.append(x)
        for y in set(by[x]['proof_refs'] + by[x]['stmt_refs']):
            stack.append(res(y))
    return seen

for root in ['thm:mai-main', 'cor:mai-gap', 'thm:gam-D', 'thm:mod-MLK', 'thm:red-main', 'thm:flh-main']:
    c = closure(root)
    print(root, 'closure size', len(c))
allc = set(closure('thm:mai-main')) | set(closure('cor:mai-gap'))
notin = [it['label'] for it in items if it['label'] and it['label'] not in allc and it['env'] in ('thm', 'lem', 'prop', 'cor', 'defn', 'conv')]
print('NOT IN CLOSURE of A/B (%d):' % len(notin), notin)
TH = ('thm', 'lem', 'prop', 'cor')
G = {}
for it in items:
    if it['label'] and it['env'] in TH:
        G[it['label']] = sorted(set(r for r in (res(y) for y in it['proof_refs']) if r and r in by and by[r]['env'] in TH and r != it['label']))
order = {it['label']: k for k, it in enumerate(items) if it['label']}
fw = [(a, b) for a in G for b in G[a] if order[b] > order[a]]
print('FORWARD proof refs (%d):' % len(fw))
for a, b in fw:
    print('   ', a, '->', b)
idx = {}
low = {}
st = []
on = set()
sccs = []
c = [0]
sys.setrecursionlimit(10000)

def sc(v):
    idx[v] = low[v] = c[0]
    c[0] += 1
    st.append(v)
    on.add(v)
    for w in G.get(v, []):
        if w not in idx:
            sc(w)
            low[v] = min(low[v], low[w])
        elif w in on:
            low[v] = min(low[v], idx[w])
    if low[v] == idx[v]:
        comp = []
        while True:
            w = st.pop()
            on.discard(w)
            comp.append(w)
            if w == v:
                break
        if len(comp) > 1:
            sccs.append(comp)

for v in G:
    if v not in idx:
        sc(v)
print('CYCLES (SCCs):', sccs)
with open('deps.txt', 'w', encoding='utf-8') as fo:
    for it in items:
        if it['label']:
            fo.write('%s [%s] %s:%d proof@%s\n  stmt-refs: %s\n  proof-refs: %s\n' % (it['label'], it['env'], it['file'], it['line'], it.get('proof_at', '-'), sorted(set(it['stmt_refs'])), sorted(set(it['proof_refs']))))
