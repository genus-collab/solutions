import sys, os, re, difflib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from closure import strip_comments, norm
NOTE = sys.argv[1]
KEYS = ['modulo the quoted', 'conditionally on', 'conditional on', 'every other step', 'enters through', 'enter only in',
        'quoted results of', 'inputs from the literature']


def paras(root, files):
    out = []
    for f in files:
        p = os.path.join(root, 'sections', f + '.tex')
        if not os.path.exists(p):
            continue
        s = strip_comments(open(p, encoding='utf-8').read())
        for para in re.split(r'\n\s*\n', s):
            # split itemize items as separate units too
            if any(k in para for k in KEYS) or re.search(r'\\item\[\((Q|QRR|M|G0|T)\)\]', para):
                out.append((f, norm(para)))
    return out


files = ['s01_intro', 's01_intro_2', 's01_intro_3', 's10_main', 's10_main_2', 's10_main_3', 's10_main_4']
O = paras(os.path.join(NOTE, 'archive_full_162pp'), files)
N = paras(NOTE, files)
Ntexts = [t for _, t in N]
for f, t in O:
    if t in Ntexts:
        print('IDENTICAL', f, t[:90]); continue
    best = max(N, key=lambda x: difflib.SequenceMatcher(None, t, x[1], autojunk=False).ratio())
    print('CHANGED', f, '->', best[0])
    wa, wb = t.split(), best[1].split()
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, wa, wb, autojunk=False).get_opcodes():
        if op != 'equal':
            print('   ', op, '| OLD:', ' '.join(wa[i1:i2])[:700], '| NEW:', ' '.join(wb[j1:j2])[:700])
print('new-only paragraphs:', [(f, t[:80]) for f, t in N if t not in [x for _, x in O]])
