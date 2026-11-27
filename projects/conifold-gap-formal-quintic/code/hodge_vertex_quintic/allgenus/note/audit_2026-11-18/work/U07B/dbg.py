exec(open('check_lemE_T.py').read().split('# Lemma T(b):')[0])
print("qT[0] low:", sorted([(k,c) for k,c in qT[0].items() if k[0]<=2 and k[1]<=2]))
print("P0vi[0][0]:", P0vi[0][0], " d:", d)
print("Pi_0(-z1):", Pi(0, mz1))
print("Pi_0(qT) low:", sorted(Pi(0,qT).items())[:3])
for a in n5:
    print(a, "qT[a] terms z^1:", sorted([(k,c) for k,c in qT[a].items() if k[0]==1 and k[1]<=2]))
