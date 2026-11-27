# -*- coding: utf-8 -*-
"""cy_log.py -- referee cyc_adversary: check bookkeeping; the log goes to logs/<name>.log next to this file
(relative path printed; no absolute path is written into any log)."""
import os, time

HERE = os.path.dirname(os.path.abspath(__file__))


class Log(object):
    def __init__(self, name):
        self.name = name; self.lines = []; self.fails = []; self.t0 = time.time()

    def __call__(self, s=''):
        print(s); self.lines.append(s)

    def check(self, label, ok, info='', control=False):
        """control=True: the check is a falsifier and MUST fail"""
        if control:
            tag = 'CONTROL FAILS (as required)' if not ok else 'CONTROL PASSED (ERROR: it must fail)'
            if ok: self.fails.append(label)
        else:
            tag = 'PASS' if ok else 'FAIL'
            if not ok: self.fails.append(label)
        self('[%s] %s%s' % (tag, label, ('   ' + info) if info else ''))
        return ok

    def close(self):
        self('time %.1f s' % (time.time() - self.t0))
        self('ALL CHECKS PASSED' if not self.fails else 'FAILURES: %s' % '; '.join(self.fails))
        os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
        with open(os.path.join(HERE, 'logs', self.name + '.log'), 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(self.lines) + '\n')
        print('log written: logs/%s.log' % self.name)
