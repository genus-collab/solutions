# pdftext.py -- minimal stdlib PDF text extractor (FlateDecode content streams, Tj/TJ/'/" operators).
# usage: py -I -B pdftext.py IN.pdf OUT.txt
import sys, re, zlib

OT1 = {0x0B: 'ff', 0x0C: 'fi', 0x0D: 'fl', 0x0E: 'ffi', 0x0F: 'ffl', 0x10: 'i', 0x7B: '-', 0x7C: '--', 0x19: 'ss',
       0x1A: 'ae', 0x1B: 'oe', 0x1C: 'o', 0x5C: '"', 0x22: '"'}


def streams(data):
    out = []
    for m in re.finditer(rb'stream\r?\n', data):
        s = m.end(); e = data.find(b'endstream', s)
        if e < 0: continue
        raw = data[s:e]
        head = data[max(0, m.start() - 600):m.start()]
        k = head.rfind(b' obj'); head = head[k:] if k >= 0 else head
        try:
            d = zlib.decompressobj().decompress(raw)
        except Exception:
            d = raw if b'FlateDecode' not in head else b''
        out.append((head, d))
    return out


def unesc(s):
    out = bytearray(); i = 0
    while i < len(s):
        c = s[i]
        if c == 0x5C and i + 1 < len(s):
            n = s[i + 1]
            if n in b'01234567':
                j = i + 1; v = 0; k = 0
                while j < len(s) and k < 3 and s[j] in b'01234567':
                    v = v * 8 + s[j] - 48; j += 1; k += 1
                out.append(v & 255); i = j; continue
            mp = {ord('n'): 10, ord('r'): 13, ord('t'): 9, ord('b'): 8, ord('f'): 12}
            if n in mp: out.append(mp[n])
            elif n in (10, 13): pass
            else: out.append(n)
            i += 2; continue
        out.append(c); i += 1
    return bytes(out)


def dec(b, twobyte):
    if twobyte:
        return ''.join(chr((b[i] << 8) | b[i + 1]) if i + 1 < len(b) else '' for i in range(0, len(b), 2))
    return ''.join(OT1.get(x, chr(x) if 32 <= x < 127 else ('?' if x >= 127 else '')) for x in b)


TOK = re.compile(rb'\((?:\\.|[^\\)])*\)|<[0-9A-Fa-f\s]*>|\[|\]|/[^\s/\[\]()<>]+|-?\d*\.?\d+|[A-Za-z\'"*]+', re.S)


def parse_paren(data, i):
    depth = 0; j = i
    while j < len(data):
        c = data[j]
        if c == 0x5C: j += 2; continue
        if c == 0x28: depth += 1
        elif c == 0x29:
            depth -= 1
            if depth == 0: return data[i + 1:j], j + 1
        j += 1
    return data[i + 1:], len(data)


THR = 7.5


def text_of(cs):
    out = []; stack = []; i = 0; n = len(cs); arr = None; yc = [0.0, 0.0, 1.0]
    while i < n:
        c = cs[i]
        if c in b' \t\r\n': i += 1; continue
        if c == 0x28:
            s, i = parse_paren(cs, i); tok = ('s', unesc(s))
        elif c == 0x3C and i + 1 < n and cs[i + 1] != 0x3C:
            j = cs.find(b'>', i); hx = re.sub(rb'\s', b'', cs[i + 1:j]); i = j + 1
            if len(hx) % 2: hx += b'0'
            tok = ('h', bytes.fromhex(hx.decode()))
        elif c == 0x5B: arr = []; i += 1; continue
        elif c == 0x5D: stack.append(('a', arr)); arr = None; i += 1; continue
        else:
            m = re.match(rb'-?\d*\.?\d+|/[^\s/\[\]()<>]+|<<|>>|[A-Za-z\'"*]+|.', cs[i:], re.S)
            t = m.group(0); i += len(t)
            if re.match(rb'-?\d*\.?\d+$', t): tok = ('n', float(t))
            elif t.startswith(b'/') or t in (b'<<', b'>>'): tok = ('o', t)
            else:
                op = t
                if op == b'Tj' and stack and stack[-1][0] in 'sh':
                    k, v = stack[-1]; out.append(dec(v, k == 'h'))
                elif op in (b"'", b'"') and stack and stack[-1][0] in 'sh':
                    k, v = stack[-1]; out.append('\n' + dec(v, k == 'h'))
                elif op == b'TJ' and stack and stack[-1][0] == 'a':
                    for k, v in stack[-1][1] or []:
                        if k in 'sh': out.append(dec(v, k == 'h'))
                        elif k == 'n' and v < -180: out.append(' ')
                elif op in (b'Td', b'TD'):
                    nums = [v for k, v in stack[-2:] if k == 'n']
                    if len(nums) == 2:
                        yc[0] += nums[1] * yc[2]
                        if abs(yc[0] - yc[1]) > THR: out.append('\n'); yc[1] = yc[0]
                        else: out.append(' ')
                elif op == b'T*': out.append('\n')
                elif op == b'Tm':
                    nums = [v for k, v in stack[-6:] if k == 'n']
                    if len(nums) == 6:
                        yc[0] = nums[5]; yc[2] = abs(nums[3]) if abs(nums[3]) > 1e-9 else 1.0
                        if abs(yc[0] - yc[1]) > THR: out.append('\n'); yc[1] = yc[0]
                        else: out.append(' ')
                elif op == b'ET': out.append(' ')
                stack = []; continue
        if arr is not None: arr.append(tok)
        else: stack.append(tok)
    return ''.join(out)


def main():
    data = open(sys.argv[1], 'rb').read()
    pages = []
    for head, d in streams(data):
        if b'BT' in d and (b'Tj' in d or b'TJ' in d) and b'/Subtype' not in head:
            try: pages.append(text_of(d))
            except Exception as ex: pages.append('[[parse error %s]]' % ex)
    txt = '\n\n=====PAGE=====\n'.join(pages)
    txt = re.sub(r'[ ]{2,}', ' ', txt)
    open(sys.argv[2], 'w', encoding='utf-8').write(txt)
    print('streams with text:', len(pages), 'chars:', len(txt))


if __name__ == '__main__':
    main()
