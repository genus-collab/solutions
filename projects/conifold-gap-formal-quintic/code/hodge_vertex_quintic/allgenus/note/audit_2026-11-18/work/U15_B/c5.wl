N0 = 5;
S0pp[s_] := 5/(s (1 - s^5));
P1[s_] := 5 s^4/(1 - s^5);
bk[k_] := BernoulliB[2 k]/(2 k (2 k - 1));
S1p[s_] := 1/(2 s) + P1[s]/2;
cn[Lv_] := Module[{CC, ser0, F3, ser1, logA, expo, ser, wick, x, S2k, psum},
  psum[r_, s_] := (D[P1[x], {x, r - 1}]/(r - 1)!) /. x -> s;
  S2k[k_, s_] := bk[k] ((5 s)^(1 - 2 k) + psum[2 k - 1, s]);
  CC = -1/S0pp[Lv];
  ser0 = Series[S0pp[Lv + t], {t, 0, 2 N0}];
  F3 = Sum[SeriesCoefficient[ser0, j - 2] t^j/(j (j - 1)), {j, 3, 2 N0 + 2}];
  ser1 = Series[S1p[Lv + t], {t, 0, 2 N0}];
  logA = Sum[SeriesCoefficient[ser1, k - 1] t^k/k, {k, 1, 2 N0}] +
     Sum[z^(2 k - 1) Normal[Series[S2k[k, Lv + t], {t, 0, 2 N0}]], {k, 1, N0}];
  expo = Expand[(F3/z + logA) /. {t -> e y, z -> e^2}];
  ser = Normal[Series[Exp[expo], {e, 0, 2 N0}]];
  wick[expr_] := Module[{ex = Expand[expr]}, (ex /. y^p_. :> If[EvenQ[p], (p - 1)!! CC^(p/2), 0]) /. y -> 0];
  Table[wick[Coefficient[ser, e, 2 n]], {n, 1, N0}]];
pts = {2, 3, 1/2, 1/3, -2};
vals = cn /@ pts;
c5 = vals[[All, 5]];
sol = Solve[Table[Sum[a[i] pts[[k]]^(5 i), {i, 0, 4}] == c5[[k]], {k, 1, 5}], Table[a[i], {i, 0, 4}]];
Print["c5 coefficients: ", sol];
Print["c5(0) target B6(1/5)/6 = ", BernoulliB[6, 1/5]/6];
Print["sum a_i (c5(1)) = ", Sum[a[i], {i, 0, 4}] /. sol];
h = Exp[-Sum[bk[k] z^(2 k - 1) (Sum[(1 - Exp[2 Pi I j/5])^(1 - 2 k), {j, 1, 4}] + (-5)^(1 - 2 k)), {k, 1, 3}]];
Print["h5 = ", FullSimplify[SeriesCoefficient[h, {z, 0, 5}]]];
Print["c4 at pts / L^16 : ", Table[vals[[k, 4]]/pts[[k]]^16, {k, 1, 5}]];
