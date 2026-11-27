import re, glob, os, sys
d = sys.argv[1]
files = sorted(glob.glob(os.path.join(d, '*.tex'))) + [os.path.join(d, '..', 'main.tex')]
labels = {}
refs = {}
BS = chr(92)
lab_re = re.compile(re.escape(BS) + r'label\{([^}]*)\}')
ref_re = re.compile(re.escape(BS) + r'(?:cref|Cref|ref|eqref|autoref|pageref|labelcref|crefrange)\*?\{([^}]*)\}')
com_re = re.compile(r'(?<!' + re.escape(BS) + r')%.*')
for f in files:
    txt = open(f, encoding='utf-8').read()
    lines = txt.split('\n')
    for i, line in enumerate(lines, 1):
        s = com_re.sub('', line)
        for m in lab_re.finditer(s):
            labels.setdefault(m.group(1), []).append((os.path.basename(f), i))
        for m in ref_re.finditer(s):
            for lab in m.group(1).split(','):
                refs.setdefault(lab.strip(), []).append((os.path.basename(f), i))
print('labels', len(labels), 'distinct refs', len(refs))
print('DUPLICATE LABELS:')
for k, v in labels.items():
    if len(v) > 1:
        print(' ', k, v)
print('REFS WITHOUT LABEL:')
for k, v in sorted(refs.items()):
    if k not in labels:
        print(' ', k, v)
print('LABELS NEVER REFERENCED:')
for k in sorted(labels):
    if k not in refs:
        print(' ', k, labels[k])
# tags (\tag) and environment-less refs
