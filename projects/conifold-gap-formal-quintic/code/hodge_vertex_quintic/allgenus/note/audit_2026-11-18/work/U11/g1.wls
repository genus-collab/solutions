v = 1; c = 1; ord = 10;
x[p_] := p^4/4 + v/p;
(* involution: s(t) with x(c+s)=x(c+t), s=-t+... *)
sAns = -t + Sum[a[k] t^k, {k, 2, ord}];
eqs = CoefficientList[Normal[Series[x[c + sAns] - x[c + t], {t, 0, ord + 1}]], t];
sol = Solve[Thread[Drop[eqs, 3] == 0], Table[a[k], {k, 2, ord}]][[1]];
iq = c + (sAns /. sol); q = c + t;
diq = D[iq, t];
(* genus 1: Res_t (1/2)(-Log[q/iq]) iq'/((q-iq)^3 x'(q)) *)
xp = D[x[p], p] /. p -> q;
expr1 = (1/2) (-Log[q/iq]) diq/((q - iq)^3 xp);
r1 = SeriesCoefficient[Series[expr1, {t, 0, 2}], -1];
Print["Loc_1[-dq/q] = ", Simplify[r1]];
(* check Loc formula pieces: chi expansion *)
Print["x(1+t) series: ", Series[x[1 + t], {t, 0, 6}]];
