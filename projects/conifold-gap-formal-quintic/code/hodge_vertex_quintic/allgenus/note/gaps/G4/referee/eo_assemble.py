# assemble the EO source (texte9.tex with its \input files) and count equation numbers (reading aid only)
# usage: py -I -B eo_assemble.py <dir> <out.tex>
import re, sys, os
d = sys.argv[1]
main = open(os.path.join(d, 'texte9.tex'), encoding='latin-1').read()


def repl(m):
    fn = m.group(1)
    p = os.path.join(d, fn)
    if not os.path.exists(p):
        p = os.path.join(d, fn + '.tex')
    return open(p, encoding='latin-1').read()


lines = []
for line in main.split('\n'):
    if line.strip().startswith('%'):
        lines.append(line); continue
    lines.append(re.sub(r'\\input\{([^}]*)\}', repl, line))
src = '\n'.join(lines)
src = src.replace('\\beq', '\\begin{equation}').replace('\\eeq', '\\end{equation}')
src = re.sub(r'\\bea(?![a-zA-Z])', r'\\begin{eqnarray}', src)
src = re.sub(r'\\eea(?![a-zA-Z])', r'\\end{eqnarray}', src)
open(sys.argv[2], 'w', encoding='latin-1').write(src)
print('written', len(src))
