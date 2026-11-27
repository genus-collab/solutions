import io, sys
p = sys.argv[1]
s = io.open(p, encoding="utf-8", newline="").read()
old = "**Lemma 5.7** [P mod [BS] prop:ResConifold at genus 1, [GP99], [M83]]. e1 e2 GW_{1,d} = P_1/d"
new = "**Lemma 5.7** [P mod [BS] prop:ResConifold at genus 1 for d >= 2, [GP99], [M83]; at d = 1 without [BS], see the direct localisation below]. e1 e2 GW_{1,d} = P_1/d"
assert s.count(old) == 1
s = s.replace(old, new)
io.open(p, "w", encoding="utf-8", newline="").write(s)
print("ok")
