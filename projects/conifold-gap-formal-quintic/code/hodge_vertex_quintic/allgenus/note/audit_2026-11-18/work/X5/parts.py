import re, glob, os, sys
d = sys.argv[1]
BS = chr(92)
files = sorted(glob.glob(os.path.join(d, '*.tex')))
texts = {os.path.basename(f): open(f, encoding='utf-8').read() for f in files}
ENVS = ['thm', 'lem', 'prop', 'cor', 'defn', 'rmk', 'conv', 'notation', 'example']
env_re = re.compile(re.escape(BS) + r'begin\{(' + '|'.join(ENVS) + r')\}')
end_re = {e: re.compile(re.escape(BS) + r'end\{' + e + r'\}') for e in ENVS}
lab_re = re.compile(re.escape(BS) + r'label\{([^}]*)\}')
# map label -> body text of environment
bodies = {}
for fn, t in texts.items():
    for m in env_re.finditer(t):
        e = m.group(1)
        me = end_re[e].search(t, m.end())
        if not me:
            continue
        body = t[m.end():me.start()]
        lm = lab_re.search(body[:400])
        if lm:
            bodies[lm.group(1)] = (fn, body)
roman = ['i', 'ii', 'iii', 'iv', 'v', 'vi', 'vii', 'viii', 'ix', 'x']
def parts_of(body):
    # explicit item labels
    explicit = set(re.findall(re.escape(BS) + r'item\[\(?([^\]\)]*)\)?\]', body))
    # count plain items inside enumerate (top-level approx)
    n_plain = len(re.findall(re.escape(BS) + r'item(?!\[)', body))
    s = set(explicit)
    for k in range(min(n_plain, 10)):
        s.add(roman[k])
    # inline labelled parts like (a) (b) in description-free text: also accept "(i)" occurrences at start
    for x in re.findall(r'\((i{1,3}|iv|v|vi{0,3}|[a-h])\)', body):
        s.add(x)
    return s
ref_re = re.compile(re.escape(BS) + r'[cC]ref\{([^}]*)\}((?:~\([a-z0-9]+\)(?:,\s*\([a-z0-9]+\))*(?:\s*(?:and|or|--)\s*\([a-z0-9]+\))*))')
problems = []
for fn, t in texts.items():
    for i, line in enumerate(t.split('\n'), 1):
        for m in ref_re.finditer(line):
            labs = [x.strip() for x in m.group(1).split(',')]
            parts = re.findall(r'\(([a-z0-9]+)\)', m.group(2))
            lab = labs[-1]
            if lab not in bodies:
                continue
            avail = parts_of(bodies[lab][1])
            for p in parts:
                if p not in avail:
                    problems.append((fn, i, lab, p, sorted(avail)))
for p in problems:
    print(p)
print('checked; problems:', len(problems))
