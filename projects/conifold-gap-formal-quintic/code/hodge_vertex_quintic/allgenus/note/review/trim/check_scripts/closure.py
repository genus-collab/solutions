import re, sys, os, json
BS = chr(92)
STMT = ['thm', 'lem', 'prop', 'cor', 'defn', 'rmk', 'notation', 'example', 'conv']


def strip_comments(s):
    return '\n'.join(re.sub(r'(?<!\\)%.*$', '', l) for l in s.split('\n'))


def expand(root, name):
    p = os.path.join(root, name if name.endswith('.tex') else name + '.tex')
    s = strip_comments(open(p, encoding='utf-8').read())
    return re.sub(r'\\input\{([^}]*)\}', lambda m: expand(root, m.group(1)), s)


REF = re.compile(r'\\(?:[cC]ref|ref|eqref|labelcref|[cC]refrange|autoref|nameref|cpageref)\*?\{([^}]*)\}')


def refs(t):
    out = []
    for m in REF.finditer(t):
        out += [x.strip() for x in m.group(1).split(',') if x.strip()]
    return out


def norm(t):
    return ' '.join(t.split())


def parse(root):
    s = expand(root, 'main')
    tok = re.compile(r'\\(begin|end)\{(' + '|'.join(STMT) + r'|proof)\}')
    stmts = {}; order = []; proofs = []; lab2stmt = {}
    stack = []
    for m in tok.finditer(s):
        if m.group(1) == 'begin':
            stack.append((m.group(2), m.end(), m.start()))
        else:
            env, b, st = stack.pop()
            assert env == m.group(2), (env, m.group(2), s[st:st + 200])
            body = s[b:m.start()]
            if env == 'proof':
                if not any(e == 'proof' for e, _, _ in stack):
                    proofs.append((st, m.end(), body))
            else:
                labs = re.findall(r'\\label\{([^}]*)\}', body)
                key = labs[0] if labs else 'UNLABELLED@%d' % st
                stmts[key] = dict(env=env, text=body, start=st, end=m.end(), labels=labs, proofs=[])
                for l in labs:
                    lab2stmt[l] = key
                order.append(key)
    order.sort(key=lambda k: stmts[k]['start'])
    for st, en, body in proofs:
        mo = re.match(r'\s*\[(.*?)\]\s', body, re.S)
        tgt = None
        if mo and refs(mo.group(1)):
            for r in refs(mo.group(1)):
                if r in lab2stmt:
                    tgt = lab2stmt[r]; break
        if tgt is None:
            prev = [k for k in order if stmts[k]['end'] <= st]
            tgt = prev[-1] if prev else None
        if tgt:
            stmts[tgt]['proofs'].append(body)
    other = {}
    for m in re.finditer(r'\\label\{([^}]*)\}', s):
        l = m.group(1)
        if l not in lab2stmt:
            other[l] = (m.start(), s)
    return s, stmts, lab2stmt, other


def env_text(s, pos):
    """text of the display environment (or paragraph) containing pos"""
    b = s.rfind(BS + 'begin{', 0, pos)
    e = s.find(BS + 'end{', pos)
    e2 = s.find('}', e)
    return s[b:e2 + 1]


def closure(stmts, lab2stmt, seeds, use_stmt_refs):
    seen = set(); todo = list(seeds); edges = {}; nonstmt = set()
    while todo:
        k = todo.pop()
        if k in seen:
            continue
        seen.add(k)
        d = stmts[k]; rs = []
        for p in d['proofs']:
            rs += refs(p)
        if use_stmt_refs:
            rs += refs(d['text'])
        deps = set()
        for r in rs:
            if r in lab2stmt:
                deps.add(lab2stmt[r])
            else:
                nonstmt.add(r)
        deps.discard(k); edges[k] = deps
        todo += [x for x in deps if x not in seen]
    return seen, edges, nonstmt


if __name__ == '__main__':
    NOTE = sys.argv[1]
    O = parse(os.path.join(NOTE, 'archive_full_162pp')); N = parse(NOTE)
    sO, stO, l2O, othO = O; sN, stN, l2N, othN = N
    seeds = ['thm:mai-main', 'cor:mai-gap']
    print('statements old', len(stO), 'new', len(stN))
    for usr in (False, True):
        cO, eO, mO = closure(stO, l2O, seeds, usr)
        cN, eN, mN = closure(stN, l2N, seeds, usr)
        print('==== mode', 'proof+statement refs' if usr else 'proof refs only',
              '| closure old', len(cO), 'new', len(cN))
        for k in sorted(cO):
            d = stO[k]
            if k not in stN:
                if k in l2N:
                    print('MOVED-LABEL', k, '-> inside', l2N[k]); continue
                print('MISSING', k, d['env']); continue
            n = stN[k]
            if norm(d['text']) != norm(n['text']):
                print('TEXTDIFF', k)
            if d['env'] != n['env']:
                print('ENVDIFF', k, d['env'], n['env'])
            if d['proofs'] and not n['proofs']:
                print('PROOFLOST', k)
            if not n['proofs'] and d['env'] in ('thm', 'lem', 'prop', 'cor'):
                print('NOPROOF-ENV (old %s)' % bool(d['proofs']), k)
            if set(d['labels']) != set(n['labels']):
                print('LABELSET', k, set(d['labels']) ^ set(n['labels']))
        print('new-only in closure:', sorted(cN - cO))
        print('old-only in closure:', sorted(cO - cN))
        for k in sorted(cO & cN):
            lost = eO[k] - eN[k]
            if lost:
                print('DEPS-DROPPED', k, sorted(lost))
        for r in sorted(mO):
            if r not in othN:
                print('NONSTMT-LABEL-MISSING', r)
            else:
                to = norm(env_text(othO[r][1], othO[r][0])); tn = norm(env_text(othN[r][1], othN[r][0]))
                if to != tn:
                    print('NONSTMT-TEXTDIFF', r)
        print('non-stmt refs old', len(mO), 'new', len(mN))
    if len(sys.argv) > 2:
        json.dump({'old': {k: [norm(stO[k]['text'])] for k in stO}, 'new': {k: [norm(stN[k]['text'])] for k in stN}},
                  open(sys.argv[2], 'w'), indent=0)
