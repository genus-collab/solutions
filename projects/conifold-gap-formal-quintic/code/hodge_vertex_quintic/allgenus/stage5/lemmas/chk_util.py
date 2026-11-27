# -*- coding: utf-8 -*-
"""chk_util.py -- logging for the Stage-5 lemma checks (logs print relative paths only)."""
import os, time


class Log:
    def __init__(self, here, name):
        self.rel = 'logs/%s.log' % name
        self.f = open(os.path.join(here, self.rel), 'w', encoding='utf-8')
        self.fails, self.t0 = [], time.time()

    def out(self, s):
        print(s); self.f.write(s + '\n'); self.f.flush()

    def check(self, name, ok, info=''):
        self.out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
        if not ok:
            self.fails.append(name)

    def control(self, name, differs, info=''):
        """a deliberate control: it must differ from the true identity (PASS = the control fails as predicted)"""
        self.check('control (must fail): ' + name, differs, info)

    def finish(self, tag):
        self.out('%s: %s   (%.0f s)' % (tag, 'ALL CHECKS PASSED' if not self.fails else 'FAILED: %s' % self.fails,
                                       time.time() - self.t0))
        self.f.close()


def first_diff(A, B, top):
    """first z-order <= top at which the z-series dicts A, B differ (None if they agree through top)"""
    for k in range(min(list(A) + list(B) + [0]), top + 1):
        if A.get(k, 0) != B.get(k, 0):
            return k
    return None
