# Minimal stdlib PDF text extractor (reading aid only). Usage: py -I -B pdftext.py in.pdf out.txt
# Decompresses FlateDecode streams and prints the strings shown by Tj/TJ/'/" operators, page by page
# (pages in the order of their content streams). Font encodings are ignored, so math is garbled.
import re, sys, zlib

def unesc(s):
    out = []; i = 0
    while i < len(s):
        c = s[i]
        if c == '\\' and i + 1 < len(s):
            d = s[i + 1]
            if d in 'nrtbf':
                out.append({'n': '\n', 'r': '', 't': ' ', 'b': '', 'f': ''}[d]); i += 2; continue
            if d in '()\\':
                out.append(d); i += 2; continue
            m = re.match(r'[0-7]{1,3}', s[i + 1:])
            if m:
                code = int(m.group(0), 8)
                out.append(chr(code) if 32 <= code < 127 else '?'); i += 1 + len(m.group(0)); continue
            i += 2; continue
        out.append(c if (32 <= ord(c) < 127) else '?'); i += 1
    return ''.join(out)

def strings_in(op_text):
    res = []; i = 0; depth = 0; cur = []
    while i < len(op_text):
        c = op_text[i]
        if depth == 0:
            if c == '(':
                depth = 1; cur = []
            elif c == ']' or c == '\n':
                pass
            elif c in '-0123456789.':
                m = re.match(r'-?\d+(\.\d+)?', op_text[i:])
                if m:
                    if float(m.group(0)) < -250: res.append(' ')
                    i += len(m.group(0)); continue
        else:
            if c == '\\':
                cur.append(op_text[i:i + 2]); i += 2; continue
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    res.append(unesc(''.join(cur))); i += 1; continue
            cur.append(c)
        i += 1
    return ''.join(res)

def page_text(content):
    out = []
    for m in re.finditer(r'(\[(?:[^\]\\]|\\.)*\]\s*TJ|\((?:[^)\\]|\\.)*\)\s*Tj|T\*|Td|TD|Tm|ET)', content, re.S):
        tok = m.group(0)
        if tok in ('T*', 'ET'):
            out.append('\n'); continue
        if tok in ('Td', 'TD', 'Tm'):
            out.append(' '); continue
        out.append(strings_in(tok))
    return ''.join(out)

def main():
    data = open(sys.argv[1], 'rb').read()
    chunks = []
    for m in re.finditer(rb'stream\r?\n', data):
        start = m.end(); end = data.find(b'endstream', start)
        if end < 0: continue
        raw = data[start:end]
        try:
            dec = zlib.decompress(raw)
        except Exception:
            try:
                dec = zlib.decompressobj().decompress(raw)
            except Exception:
                continue
        txt = dec.decode('latin-1')
        if 'BT' in txt and ('Tj' in txt or 'TJ' in txt):
            chunks.append(page_text(txt))
    with open(sys.argv[2], 'w', encoding='ascii', errors='replace') as f:
        for k, c in enumerate(chunks):
            f.write('\n===== chunk %d =====\n' % k); f.write(c)
    print('chunks', len(chunks))

main()
