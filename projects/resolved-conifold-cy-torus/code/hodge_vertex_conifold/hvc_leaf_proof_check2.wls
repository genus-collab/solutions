(* hvc_leaf_proof_check2.wls -- Mathematica checks requested by leaf_proof/leaf_structure_proof.md (s.12.3: M8, M10, M13-M15) and by
   its completeness critic (Q1, Q2), against the stored factors. Conventions of hvc_lib.m: eps1 = x, eps2 = 1, c = 1 + x, t = Log Q,
   tau = 1 - e^-t, sst = (1 - U1)/tau. Conjugate encoding: rc_lib.m stores the data of idempotent 2 as cconj of those of idempotent 1, so
   the true value of every stored factor is its value at the singular root U1; its value at U2 is the conjugate quantity.
   (Q2)  E_22 in the point frame and in frame A (genus-3 factors, z^6): t-orders at U1 (the true E_22) and at U2 (what line 24 of
         hvc_asymptotics.wls and hvc_asymptotics_g4.wls printed), m, n <= 2; E_11, E_12, E_21, Dil_2 for comparison.
   (M14) ratCN[a,2] = cconj ratCN[a,1], DilN[2] = cconj DilN[1], EdN[2,2] = cconj EdN[1,1], EdN[2,1] = cconj EdN[1,2], exactly;
         E_22 at t = 0 against (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w) (Proposition 8.6(c)), m + n <= 8; E_22(0,0) = -x^2 (1+x+x^2)/(6 c^2);
         the 2x2 minor of Remark 8.8(1); EdN[2,2] read at U2 equals the E_11 value -(x^3/2)/c^(m+n+2) (the erratum to F3 of the note).
   (M15) E_21(m, n) = E_12(n, m) exactly (m + n <= 8); smallest t-order at U1 of Dil_i and E_ij (expected >= 0).
   (M8)  v_2(z; 0) = h_w(z)^-2 (1, x/c) to z^9, [t^1] r_{0,2}[m] = 0, Delta_2 = x^3 delta sst^2 to t^10, Delta_2(0) = x^3/c, Delta_2'(0) = 0.
   (Q1)  determinant of the E_11 block over the psi-orders {0, 1}: expected 0, 0, (1-x)^4 x^6/(18 c^10) at t^0, t^1, t^2.
   (M10b) Theorem 11.4: rho^2 E_11(psi, psi') = Sum_{i,l} e_il(t) lam_i(psi) lam_l(psi'), e_il = (i l/(c^2 t^2)) Int_0^t Theta Ahat_i' Ahat_l',
         against EdN[1,1] (m + n <= 8, t^0..t^KT).
   (M13) (9.6): rho E_12(psi, w) = -Sum_j j <S_j(t), v_2(-w; t)> lam_j(psi) lam_j(-w), S_j from Corollary 11.2(f), against EdN[1,2]
         (m + n <= 8, t^0..t^KT2).
   Usage: wolframscript -file hvc_leaf_proof_check2.wls [KT] [KT2]     Output: hvc_leaf_proof_check2_out.txt (rewritten after every line) *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine]; KT = If[Length[args] >= 1, ToExpression[args[[1]]], 6]; KT2 = If[Length[args] >= 2, ToExpression[args[[2]]], 5];
Get["hvc_lib.m"]; t0 = AbsoluteTime[];
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]; Export["hvc_leaf_proof_check2_out.txt", StringRiffle[out, "\n"], "Text"]);
tm[] := " [" <> ToString[Round[AbsoluteTime[] - t0]] <> " s]";
c = 1 + x;
tidy[s_SeriesData] := MapAt[Together, s, 3]; tidy[s_] := Together[s];
evC[cc_C] := cc[[1]] + cc[[2]] U;
zeroC[C[a_, b_]] := Together[a] === 0 && Together[b] === 0;
diffC[a_C, b_C] := zeroC[cadd[a, cscale[-1, b]]];
(* leading t-order of a series ("zero" if every coefficient up to the truncation vanishes) *)
ordS[s_SeriesData] := Module[{cs = Together /@ s[[3]], p = 0}, Do[If[cs[[i]] =!= 0, p = i; Break[]], {i, Length[cs]}]; If[p == 0, "zero", (s[[4]] + p - 1)/s[[6]]]];
(* t^k coefficient of a series that may have collapsed to a constant *)
sc[s_SeriesData, k_] := SeriesCoefficient[s, k]; sc[s_, k_] := If[k == 0, s, 0];
ordS[s_] := If[Together[s] === 0, "zero", 0];
minOrd[l_List] := Module[{n = DeleteCases[Flatten[l], "zero"]}, If[n === {}, "all zero", Min[n]]];
(* exact roots (the branch with sq -> (1 + x)/x at Q = 1, so U1 -> 1, U2 -> -1/x) *)
roots[] := Module[{S0 = Together[SsumQ], sq}, sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; {(S0 + sq)/2, (S0 - sq)/2}];
serAt[cc_C, root_, n_] := tidy[Series[evC[cc] /. U -> root /. Q -> Exp[t], {t, 0, n}]];

(* ---------- (Q2) point frame and frame A, genus-3 factors ---------- *)
Do[Get["hvc_factors_g3_" <> tag <> ".mx"]; Ssum = SsumQ; Pprod = PprodQ;
  Module[{r = roots[], tab},
   tab[ij_, rt_] := Table[ordS[serAt[EdN[ij[[1]], ij[[2]]][m, n], rt, 3]], {m, 0, 2}, {n, 0, 2}];
   emit["(Q2) frame " <> tag <> " (genus-3 factors, ZMAX = " <> ToString[ZMAX] <> "), t-orders for m, n = 0..2: E_22 at U1 (its true value) " <> ToString[tab[{2, 2}, r[[1]]]] <>
     "; EdN[2,2] read at U2 (as printed by hvc_asymptotics*.wls) " <> ToString[tab[{2, 2}, r[[2]]]] <> "; E_11 at U1 " <> ToString[tab[{1, 1}, r[[1]]]] <>
     "; E_12 at U1 " <> ToString[tab[{1, 2}, r[[1]]]] <> "; E_21 at U1 " <> ToString[tab[{2, 1}, r[[1]]]] <>
     "; Dil_2(n), n = 1..6, at U1 " <> ToString[Table[ordS[serAt[DilN[2][n], r[[1]], 3]], {n, 1, 6}]] <> tm[]]],
  {tag, {"point", "Am"}}];

(* ---------- frame B, genus-4 factors (z^9) ---------- *)
Get["hvc_factors_g4_Bm.mx"]; ZM = ZMAX; Ssum = SsumQ; Pprod = PprodQ;
{U1, U2} = roots[];
ser1[cc_C, n_] := serAt[cc, U1, n];

(* (M14) conjugate encoding *)
conjRows = And @@ Flatten[Table[diffC[ratCN[a, 2][m], cconj[ratCN[a, 1][m]]], {a, 0, 1}, {m, 0, ZM}]];
conjDil = And @@ Table[diffC[DilN[2][n], cconj[DilN[1][n]]], {n, 1, ZM + 1}];
conjE22 = And @@ Flatten[Table[If[m + n <= ZM - 1, diffC[EdN[2, 2][m, n], cconj[EdN[1, 1][m, n]]], True], {m, 0, ZM - 1}, {n, 0, ZM - 1}]];
conjE21 = And @@ Flatten[Table[If[m + n <= ZM - 1, diffC[EdN[2, 1][m, n], cconj[EdN[1, 2][m, n]]], True], {m, 0, ZM - 1}, {n, 0, ZM - 1}]];
emit["(M14a) conjugate encoding, exact in the quadratic algebra: ratCN[a,2] = cconj ratCN[a,1] (a = 0, 1; z^0..z^" <> ToString[ZM] <> "): " <> ToString[conjRows] <>
  "; DilN[2] = cconj DilN[1]: " <> ToString[conjDil] <> "; EdN[2,2] = cconj EdN[1,1]: " <> ToString[conjE22] <> "; EdN[2,1] = cconj EdN[1,2]: " <> ToString[conjE21] <> tm[]];
(* E_22 at t = 0 *)
h2c = hCoeffs[{c, -x, -1}, 2, 12]; hm2c = hCoeffs[{c, -x, -1}, -2, 12];
H2[v_] := Sum[h2c[[j + 1]] v^j, {j, 0, 12}];
num9 = Expand[Normal[Series[1 - H2[lm z] H2[lm w], {lm, 0, 9}]] /. lm -> 1];
quot = PolynomialQuotient[num9, z + w, z]; rem = Together[PolynomialRemainder[num9, z + w, z]];
E22pred[m_, n_] := Together[(x^3/c) Coefficient[Coefficient[Expand[quot], z, m], w, n]];
E22ser = Table[If[m + n <= ZM - 1, ser1[EdN[2, 2][m, n], 2], 0], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
E22t0 = Table[If[m + n <= ZM - 1, Together[sc[E22ser[[m + 1, n + 1]], 0]], 0], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
okE22 = And @@ Flatten[Table[If[m + n <= ZM - 1, Together[E22t0[[m + 1, n + 1]] - E22pred[m, n]] === 0, True], {m, 0, ZM - 1}, {n, 0, ZM - 1}]];
St[k_] := 1 + x^(1 - 2 k) - c^(1 - 2 k);
minorPred = (x^6/c^2) St[1] (72 St[2] + 5 St[1]^3)/77760;
minorVal = Together[E22t0[[1, 1]] E22t0[[2, 2]] - E22t0[[1, 2]] E22t0[[2, 1]]];
E22U2 = Table[If[m + n <= ZM - 1, Together[sc[serAt[EdN[2, 2][m, n], U2, 0], 0]], 0], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
okU2 = And @@ Flatten[Table[If[m + n <= ZM - 1, Together[E22U2[[m + 1, n + 1]] + (x^3/2)/c^(m + n + 2)] === 0, True], {m, 0, ZM - 1}, {n, 0, ZM - 1}]];
emit["(M14b) E_22 at t = 0 (read at U1) = (x^3/c)[z^m w^n](1 - h_w(z)^2 h_w(w)^2)/(z + w) for all m + n <= " <> ToString[ZM - 1] <> " (division remainder " <> ToString[rem] <> "): " <> ToString[okE22] <>
  "; E_22(0,0) = -x^2 (1+x+x^2)/(6 c^2): " <> ToString[Together[E22t0[[1, 1]] + x^2 (1 + x + x^2)/(6 c^2)] === 0] <>
  "; 2x2 minor = (x^6/c^2) St_1 (72 St_2 + 5 St_1^3)/77760: " <> ToString[Together[minorVal - minorPred] === 0] <> ", value at x = 1: " <> ToString[minorVal /. x -> 1] <>
  "; EdN[2,2] read at U2 = -(x^3/2)/c^(m+n+2) (the E_11 value; what F3 of conifold_regularity_note.md reports for E_22): " <> ToString[okU2] <> tm[]];
emit["      E_22(m, n) at t = 0, m + n <= 2: " <> ToString[InputForm[Table[If[m + n <= 2, Factor[E22t0[[m + 1, n + 1]]], "-"], {m, 0, 2}, {n, 0, 2}]]]];

(* (M15) transposition and regularity *)
transp = And @@ Flatten[Table[If[m + n <= ZM - 1, diffC[EdN[2, 1][m, n], EdN[1, 2][n, m]], True], {m, 0, ZM - 1}, {n, 0, ZM - 1}]];
E11ser = Table[If[m + n <= ZM - 1, ser1[EdN[1, 1][m, n], KT], 0], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
E12ser = Table[If[m + n <= ZM - 1, ser1[EdN[1, 2][m, n], KT2], 0], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
E21ser = Table[If[m + n <= ZM - 1, ser1[EdN[2, 1][m, n], 2], 0], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
ordsE[tab_] := minOrd[Table[If[m + n <= ZM - 1, ordS[tab[[m + 1, n + 1]]], "zero"], {m, 0, ZM - 1}, {n, 0, ZM - 1}]];
dilOrd[i_] := minOrd[Table[ordS[ser1[DilN[i][n], 2]], {n, 1, ZM + 1}]];
emit["(M15) E_21(m, n) = E_12(n, m) exactly in the quadratic algebra, m + n <= " <> ToString[ZM - 1] <> ": " <> ToString[transp] <>
  "; smallest t-order at U1 of Dil_1, Dil_2 (d = 1.." <> ToString[ZM + 1] <> "): " <> ToString[{dilOrd[1], dilOrd[2]}] <>
  "; of E_11, E_12, E_21, E_22 (m + n <= " <> ToString[ZM - 1] <> "): " <> ToString[ordsE /@ {E11ser, E12ser, E21ser, E22ser}] <> tm[]];

(* (M8) the regular rows at t = 0 and Delta_2 *)
r02 = Table[ser1[ratCN[0, 2][m], 2], {m, 0, ZM}]; r12 = Table[ser1[ratCN[1, 2][m], 2], {m, 0, ZM}];
m8a = Table[{ordS[r02[[m + 1]]], ordS[r12[[m + 1]]], Together[sc[r02[[m + 1]], 0] - hm2c[[m + 1]]] === 0,
    Together[sc[r12[[m + 1]], 0] - (x/c) hm2c[[m + 1]]] === 0, Together[sc[r02[[m + 1]], 1]] === 0}, {m, 0, ZM}];
tauE = 1 - Exp[-t]; sstE = (1 - (U1 /. Q -> Exp[t]))/tauE; delE = c Sqrt[1 - 4 x tauE/c^2];
d2ser = ser1[DeltaCN[2], 10]; d2pred = tidy[Series[x^3 delE sstE^2, {t, 0, 10}]];
emit["(M8) regular rows, m = 0.." <> ToString[ZM] <> ": {t-order of r_{0,2}[m], of r_{1,2}[m], r_{0,2}[m](0) = [z^m] h_w^-2, r_{1,2}[m](0) = (x/c)[z^m] h_w^-2, [t^1] r_{0,2}[m] = 0}: " <>
  ToString[Union[m8a]] <> "; Delta_2 = x^3 delta sst^2 through t^10: " <> ToString[Normal[tidy[d2ser - d2pred]] === 0] <>
  "; Delta_2(0) = x^3/c, Delta_2'(0) = 0: " <> ToString[{Together[SeriesCoefficient[d2ser, 0] - x^3/c] === 0, Together[SeriesCoefficient[d2ser, 1]] === 0}] <> tm[]];

(* (Q1) the E_11 block over psi-orders {0, 1} *)
detS = tidy[Series[Normal[E11ser[[1, 1]]] Normal[E11ser[[2, 2]]] - Normal[E11ser[[1, 2]]] Normal[E11ser[[2, 1]]], {t, 0, 3}]];
emit["(Q1) det of the E_11 block over psi-orders {0, 1}, t^0..t^3: " <> ToString[InputForm[Table[Factor[SeriesCoefficient[detS, k]], {k, 0, 3}]]] <>
  "; t^2 coefficient - (1-x)^4 x^6/(18 c^10) = " <> ToString[Together[SeriesCoefficient[detS, 2] - (1 - x)^4 x^6/(18 c^10)]] <> tm[]];

(* closed forms: A_j, Theta, rho, l *)
NK = Max[KT, KT2] + 3;
tauS = Normal[Series[1 - Exp[-t], {t, 0, NK + 2}]];
sst = Normal[tidy[Series[sstE, {t, 0, NK + 2}]]];
dphi = Integrate[Normal[tidy[Series[1 - x sst, {t, 0, NK + 1}]]], {t, 0, t}];          (* dphi' = 1 - x sst, dphi(0) = 0 *)
hyp[j_] := Sum[Pochhammer[j x/c, m]^2/(Pochhammer[j, m] m!) tauS^m, {m, 0, NK + 1}];
Cj[j_] := c^2 Pochhammer[1 - j/c, j - 1]^2/(j! (j - 1)!);
zetaS = Normal[tidy[Series[tauS Exp[dphi/c], {t, 0, NK + 2}]]];
Ahat[j_] := Ahat[j] = Normal[tidy[Series[(-1)^j Cj[j] zetaS^j hyp[j], {t, 0, NK + 1}]]];   (* Ahat_j = t A_j *)
ThetaS = Normal[tidy[Series[-x^3 tauS sst^2, {t, 0, NK + 1}]]];
rhoS = Normal[tidy[Series[(tauS/t) 2/((1 - uu)^(1/4) (1 + (1 - uu)^(1/2))) /. uu -> 4 x tauS/c^2, {t, 0, NK}]]];
rhoD = tidy[Series[(Normal[ser1[DeltaCN[1], KT + 2]] t^2/(-x c^3))^(-1/2), {t, 0, KT}]];
emit["      rho: the branch-free closed form of Corollary 11.2(d) equals (Delta_1 t^2/(-x c^3))^(-1/2) through t^" <> ToString[KT] <> ": " <>
  ToString[Normal[tidy[rhoD - Series[rhoS, {t, 0, KT}]]] === 0]];

(* (M10b) Theorem 11.4 against EdN[1,1] *)
eil[i_, l_] := eil[i, l] = Module[{ig = Normal[tidy[Series[ThetaS D[Ahat[i], t] D[Ahat[l], t], {t, 0, KT + 1}]]]},
    Expand[(i l/c^2) Integrate[ig, {t, 0, t}]/t^2]];
rho2 = Normal[tidy[Series[rhoS^2, {t, 0, KT}]]];
m10b = Table[If[m + n <= ZM - 1, Module[{lhs = tidy[Series[rho2 Normal[E11ser[[m + 1, n + 1]]], {t, 0, KT}]]},
      And @@ Table[Together[SeriesCoefficient[lhs, k] - Sum[If[i + l - 2 <= k, Coefficient[eil[i, l], t, k] i^m l^n/c^(m + n + 2), 0], {i, 1, k + 1}, {l, 1, k + 1}]] === 0, {k, 0, KT}]], True],
   {m, 0, ZM - 1}, {n, 0, ZM - 1}];
emit["(M10b) Theorem 11.4, rho^2 E_11 = Sum e_il(t) lam_i lam_l' with e_il = (i l/(c^2 t^2)) Int_0^t Theta Ahat_i' Ahat_l', against EdN[1,1] for m + n <= " <> ToString[ZM - 1] <>
  ", t^0..t^" <> ToString[KT] <> ": " <> ToString[And @@ Flatten[m10b]] <> tm[]];

(* (M13) (9.6) against EdN[1,2] *)
lS = Normal[tidy[Series[c (1/tauS - 1/t) - x sst, {t, 0, NK}]]];
Aj[j_] := Expand[Ahat[j]/t];
Sj[j_] := Sj[j] = Expand[{-(j/c) Aj[j], D[Aj[j], t] - ((j - 1)/t) Aj[j] - (j/c) lS Aj[j]}];
v2 = Table[{Normal[ser1[ratCN[0, 2][m], KT2]], Normal[ser1[ratCN[1, 2][m], KT2]]}, {m, 0, ZM}];
pr[a_, b_] := x^2 (a[[1]] b[[2]] + a[[2]] b[[1]]) - c x a[[2]] b[[2]];
pairS[j_, mm_] := pairS[j, mm] = Expand[pr[Sj[j], v2[[mm + 1]]]];
pred12[a_, b_] := Expand[-Sum[j (j^a/c^(a + 1)) Sum[(-1)^mm pairS[j, mm] (-j)^(b - mm)/c^(b - mm + 1), {mm, 0, b}], {j, 1, KT2 + 2}]];
m13 = Table[If[m + n <= ZM - 1, Module[{lhs = tidy[Series[Normal[rhoS] Normal[E12ser[[m + 1, n + 1]]], {t, 0, KT2}]], p = pred12[m, n]},
      And @@ Table[Together[SeriesCoefficient[lhs, k] - Coefficient[p, t, k]] === 0, {k, 0, KT2}]], True], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
emit["(M13) (9.6), rho E_12(psi, w) = -Sum_j j <S_j(t), v_2(-w; t)> lam_j(psi) lam_j(-w), against EdN[1,2] for m + n <= " <> ToString[ZM - 1] <> ", t^0..t^" <> ToString[KT2] <> ": " <>
  ToString[And @@ Flatten[m13]] <> tm[]];
emit["total time " <> ToString[Round[AbsoluteTime[] - t0]] <> " s"];
