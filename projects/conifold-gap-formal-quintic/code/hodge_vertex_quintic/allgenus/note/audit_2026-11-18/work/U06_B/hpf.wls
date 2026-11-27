(* sanity check of lem:gau-poly on eq:gau-HPF: polynomiality and resonance, with symbolic free coefficients *)
Oe[e_][F_] := F + zh (1 - u) ((4/5) zh D[F, zh] - (e/5) F - u D[F, u]);
lhs[Phi_] := Module[{P1, P2, P3, P4, P5},
  P1 = Oe[-1][Phi]; P2 = Oe[-2][P1]; P3 = Oe[-3][P2]; P4 = Oe[-4][P3]; P5 = Oe[-5][P4];
  P5 - Phi - (1 - u) (3 zh P4 + 17/5 u zh^2 P3 + 9/5 u^2 zh^3 P2 + 274/625 u^3 zh^4 P1 + 24/625 u^4 zh^5 Phi)];
NN = 16;
ph[0] = 1;
Do[
  Phi = Sum[ph[m] zh^m, {m, 0, n - 1}];
  cf = Coefficient[Expand[lhs[Phi]], zh, n + 1];
  rn = Together[cf/(1 - u)];
  If[! PolynomialQ[rn, u], Print["r_", n, " not polynomial: ", rn]];
  rn = Expand[rn];
  deg = Exponent[rn, u];
  new = 0;
  Do[
    If[5 k == 4 n,
      res = Coefficient[rn, u, k];
      Print["n=", n, " resonance coefficient [u^", k, "] r_n = ", Simplify[res]];
      new += fc[n/5] u^k,
      new += -Coefficient[rn, u, k]/(4 n - 5 k) u^k],
    {k, 0, Max[deg, 4 n/5]}];
  ph[n] = Expand[new];
  , {n, 1, NN}];
Print["ph1 = ", ph[1]];
Print["ph2 = ", ph[2]];
Print["ph3 = ", ph[3]];
Print["degrees: ", Table[Exponent[ph[n], u], {n, 1, NN}]];
(* check that the solution indeed satisfies HPF through zh^(NN) *)
Phi = Sum[ph[m] zh^m, {m, 0, NN}];
chk = Table[Simplify[Coefficient[Expand[lhs[Phi]], zh, j]], {j, 0, NN}];
Print["HPF residuals through zh^NN: ", chk];
