import sys, re, html
from html.parser import HTMLParser

class P(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out = []
        self.skip = 0  # depth inside math
        self.skiptags = 0
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'math':
            alt = a.get('alttext', '')
            self.out.append(' $' + alt + '$ ')
            self.skip += 1
            return
        if tag in ('script', 'style'):
            self.skiptags += 1
        if tag in ('p', 'div', 'br', 'tr', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'li', 'table'):
            self.out.append('\n')
    def handle_endtag(self, tag):
        if tag == 'math':
            self.skip -= 1
            return
        if tag in ('script', 'style'):
            self.skiptags -= 1
        if tag in ('p', 'div', 'tr', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'li'):
            self.out.append('\n')
        if tag in ('td',):
            self.out.append(' ')
    def handle_data(self, data):
        if self.skip or self.skiptags:
            return
        self.out.append(data)

src, dst = sys.argv[1], sys.argv[2]
with open(src, encoding='utf-8', errors='replace') as f:
    s = f.read()
p = P()
p.feed(s)
t = ''.join(p.out)
t = re.sub(r'[ \t]+', ' ', t)
t = re.sub(r'\n\s*\n+', '\n\n', t)
with open(dst, 'w', encoding='utf-8') as f:
    f.write(t)
