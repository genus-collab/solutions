# number the displayed equations of EO appendix A (appproofproperties.tex) up to a given line (reading aid)
# each \eeq or \eea closes one numbered display (EO uses \cr inside eqnarray, so one number per display)
import re, sys
lines = open(sys.argv[1], encoding='latin-1').read().split('\n')
upto = int(sys.argv[2])
n = 0
for k, line in enumerate(lines[:upto], 1):
    s = line.split('%')[0] if not line.strip().startswith('%') else ''
    for m in re.finditer(r'\\eeq(?![a-zA-Z])|\\eea(?![a-zA-Z])|\\nonumber|\\\\', s):
        t = m.group(0)
        if t in ('\\eeq', '\\eea'):
            n += 1
            print('A-%d closes at line %d: %s' % (n, k, lines[k - 1][:70].strip()))
        else:
            print('   note: %s at line %d' % (t, k))
