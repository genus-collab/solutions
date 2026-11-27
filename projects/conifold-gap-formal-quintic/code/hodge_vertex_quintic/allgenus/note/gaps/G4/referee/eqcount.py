# count numbered equations per section in a LaTeX file (approximate reading aid; not evidence)
# usage: py -I -B eqcount.py file.tex [mode]; mode 'sec' (default): numbering section-equation
import re, sys
src = open(sys.argv[1], encoding='latin-1').read()
# strip comments
src = re.sub(r'(?<!\\)%[^\n]*', '', src)
start = src.find('\\begin{document}')
body = src[start:]
sec = 0; eq = 0; app = False
pat = (r'\\section\*?\{|\\begin\{(equation|align|eqnarray|gather|multline)\}|'
       r'\\end\{(equation|align|eqnarray|gather|multline)\}|\\beq(?![a-zA-Z])|\\eeq(?![a-zA-Z])|'
       r'\\label\{([^}]*)\}|\\\\|\\notag|\\nonumber|\\appendix|\\begin\{(equation|align|eqnarray|gather|multline)\*\}|'
       r'\\end\{(equation|align|eqnarray|gather|multline)\*\}')
inenv = None; star = False; linenum = True; out = []; curlabels = []
for m in re.finditer(pat, body):
    t = m.group(0)
    if t.startswith('\\section'):
        if t.startswith('\\section*'):
            continue
        sec += 1; eq = 0; continue
    if t == '\\appendix':
        app = True; sec = 0; out.append(('APPENDIX', 'start')); continue
    if m.group(4):
        inenv = m.group(4); star = True; continue
    if m.group(5):
        inenv = None; star = False; continue
    if t.startswith('\\begin') or t == '\\beq':
        inenv = m.group(1) or 'equation'; star = False; linenum = True; curlabels = []; continue
    if inenv and star:
        continue
    if inenv and (t == '\\notag' or t == '\\nonumber'):
        linenum = False; continue
    if inenv and t.startswith('\\label'):
        curlabels.append(m.group(3)); continue
    if inenv and t == '\\\\' and inenv in ('align', 'eqnarray', 'gather'):
        if linenum:
            eq += 1; out.append(('%s%d-%d' % ('A' if app else '', sec, eq), ','.join(curlabels)))
        linenum = True; curlabels = []; continue
    if t.startswith('\\end') or t == '\\eeq':
        if inenv and linenum:
            eq += 1; out.append(('%s%d-%d' % ('A' if app else '', sec, eq), ','.join(curlabels)))
        inenv = None; linenum = True; curlabels = []; continue
for o in out:
    print(o[0], o[1])
