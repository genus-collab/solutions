import io, os, shutil, sys
BK = os.path.dirname(os.path.abspath(__file__)) + '/backup'
os.makedirs(BK, exist_ok=True)

class Doc:
    def __init__(self, path):
        self.path = path
        b = BK + '/' + os.path.basename(path)
        if not os.path.exists(b):
            shutil.copyfile(path, b)
        with io.open(path, 'r', encoding='utf-8', newline='') as f:
            self.s = f.read()
        self.n = 0
        self.fail = []

    def rep(self, old, new, count=1):
        c = self.s.count(old)
        if c != count:
            self.fail.append((c, old[:110]))
            return
        self.s = self.s.replace(old, new)
        self.n += 1

    def save(self):
        if self.fail:
            for c, o in self.fail:
                print('FAIL', c, ascii(o))
            print('NOT SAVED', self.path, len(self.fail), 'failures of', self.n + len(self.fail))
            return False
        with io.open(self.path, 'w', encoding='utf-8', newline='') as f:
            f.write(self.s)
        print('saved', self.path, self.n, 'replacements')
        return True
