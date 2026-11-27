"""lf_log: logging and check bookkeeping for the lmi_formal scripts (relative paths only)."""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


class Log:
    def __init__(self, name):
        self.rel = "logs/%s.log" % name
        os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
        self.fh = open(os.path.join(HERE, "logs", name + ".log"), "w", encoding="utf-8")
        self.n_pass = 0
        self.n_fail = 0
        self.t0 = time.time()
        self.say("log: %s   (python %s)" % (self.rel, sys.version.split()[0]))

    def say(self, msg):
        print(msg, flush=True)
        self.fh.write(msg + "\n")
        self.fh.flush()

    def check(self, label, ok, detail=""):
        """A check that must hold."""
        if ok:
            self.n_pass += 1
            self.say("PASS  %s%s" % (label, ("  [" + detail + "]") if detail else ""))
        else:
            self.n_fail += 1
            self.say("FAIL  %s%s" % (label, ("  [" + detail + "]") if detail else ""))

    def control(self, label, identity_holds, detail=""):
        """A failing control: the identity must FAIL; passing means it failed."""
        if not identity_holds:
            self.n_pass += 1
            self.say("PASS  control: %s -- identity FAILS as expected%s"
                     % (label, ("  [" + detail + "]") if detail else ""))
        else:
            self.n_fail += 1
            self.say("FAIL  control: %s -- identity HOLDS although it should fail%s"
                     % (label, ("  [" + detail + "]") if detail else ""))

    def close(self):
        dt = time.time() - self.t0
        self.say("checks: %d passed, %d failed; %.1f s" % (self.n_pass, self.n_fail, dt))
        self.say("ALL CHECKS PASSED" if self.n_fail == 0 else "SOME CHECKS FAILED")
        self.fh.close()


def first_diff(a, b):
    """Index of the first difference of two coefficient lists (None if equal on the common range)."""
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return i
    return None
