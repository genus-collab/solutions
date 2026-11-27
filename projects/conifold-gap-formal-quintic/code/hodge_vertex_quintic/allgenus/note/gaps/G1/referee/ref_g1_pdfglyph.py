# -*- coding: utf-8 -*-
"""ref_g1_pdfglyph.py -- referee record: which glyph is printed in the exponent of Givental's Cor. 11.7 (alg-geom/9603021v2, p.41).

Reads a TeX-generated PDF as DATA only (zlib-decompressed content streams, Type1 font names), no rendering, no execution.
usage: py -I -B ref_g1_pdfglyph.py PATH_TO_alg-geom-9603021v2.pdf      (the path is not printed; writes ref_g1_pdfglyph.log next to itself)
In the cmmi fonts the letter l is code 0x6C ('l') and lambda is code 0x15; alpha is 0x0B.  The script prints, for pages 40 and 41, the
font-tagged text of the multiplication factors of Thm 11.6 (2) and Cor. 11.7 (2), with control codes written as \\xNN.
"""
import os, re, sys, zlib

HERE = os.path.dirname(os.path.abspath(__file__))
LOGNAME = "ref_g1_pdfglyph.log"


def page_texts(path, want):
    data = open(path, 'rb').read()
    objs = {int(m.group(1)): m.group(2) for m in re.finditer(rb'(\d+)\s+0\s+obj(.*?)endobj', data, re.S)}
    def stream(o):
        body = objs.get(o, b'')
        m = re.search(rb'stream\r?\n(.*?)endstream', body, re.S)
        if not m:
            return b''
        raw = m.group(1)
        return zlib.decompressobj().decompress(raw) if b'FlateDecode' in body[:m.start()] else raw
    root = [k for k, b in objs.items() if re.search(rb'/Type\s*/Pages', b) and not re.search(rb'/Parent', b)][0]
    pages = []
    def walk(k):
        b = objs[k]
        if re.search(rb'/Type\s*/Pages', b):
            for kk in re.findall(rb'(\d+)\s+0\s+R', re.search(rb'/Kids\s*\[(.*?)\]', b, re.S).group(1)):
                walk(int(kk))
        else:
            pages.append(k)
    walk(root)
    out = {}
    for p in want:
        pb = objs[pages[p - 1]]
        res = pb
        m = re.search(rb'/Resources\s+(\d+)\s+0\s+R', pb)
        if m:
            res = objs[int(m.group(1))]
        fd = re.search(rb'/Font\s*<<(.*?)>>', res, re.S)
        fd = fd.group(1) if fd else b''
        fi = re.search(rb'/Font\s+(\d+)\s+0\s+R', res)
        if fi:
            fd = objs[int(fi.group(1))]
        fonts = {}
        for name, ref in re.findall(rb'/(\w+)\s+(\d+)\s+0\s+R', fd):
            bf = re.search(rb'/BaseFont\s*/([\w+\-]+)', objs.get(int(ref), b''))
            fonts[name.decode()] = bf.group(1).decode().split('+')[-1] if bf else '?'
        cont = re.search(rb'/Contents\s+(\d+)\s+0\s+R', pb)
        cs = stream(int(cont.group(1))).decode('latin-1') if cont else ''
        s = ''
        for tok in re.finditer(r'/(\w+)\s+[\d.]+\s+Tf|\[(.*?)\]\s*TJ', cs, re.S):
            if tok.group(1):
                s += ' <%s> ' % fonts.get(tok.group(1), tok.group(1))
            else:
                s += ''.join(re.findall(r'\(((?:[^()\\]|\\.)*)\)', tok.group(2)))
        out[p] = ''.join(ch if 32 <= ord(ch) < 127 else '\\x%02x' % ord(ch) for ch in s)
    return out


def main():
    txt = page_texts(sys.argv[1], (40, 41))
    lines = ["ref_g1_pdfglyph: font-tagged text (cmmi: l = 'l', lambda = \\x15, alpha = \\x0b)"]
    i = txt[40].find('2\\)multiply')
    lines.append("p.40 Thm 11.6 (2): " + txt[40][i:i + 720])
    j = txt[41].find('2\\)multiplicationby')
    lines.append("p.41 Cor. 11.7 (1)-(2): " + txt[41][txt[41].find('1\\)change'):j + 330])
    for l in lines:
        print(l)
    with open(os.path.join(HERE, LOGNAME), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("log written to " + LOGNAME)


main()
