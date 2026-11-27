(* hvc_leaf_proof_check.wls -- checks of the all-order proof of the leaf structure (Lemma C) against the stored frame-B factors
   (hvc_factors_g4_Bm.mx, z^9). Conventions of hvc_lib.m: eps1 = x, eps2 = 1, c = 1 + x, t = Log Q, tau = 1 - e^-t, singular root U1.
   Claims checked:
   (1) Delta_1 = x^2 W2(U1), W2 = (U d/dU)^2 W, W = x Log[1 - Q U] - Log[U] + Log[1 - U]  (identity in the quadratic algebra; also Delta_2 at U2);
   (2) exact singular row: G(z; t) := rho(t) r_{0,1}(z; t) = c/(c+z) (tau/t) 2F1(1/z+1, 1/z+1; c/z+2; tau) Exp[-dphi(tau)/z],
       dphi = phi_tau(sst) - phi_0(1/c), phi_tau(s) = Log[s] + x Log[1-s] - Log[1 - tau s], sst = (1 - U1)/tau;
       and rho = tau/(t alg), alg = (1 - tau sst) Sqrt[phi2_tau(sst)/phi2_0(1/c)], phi2 = second s-derivative of phi;
   (3) L2 = x sst, L1 = c/tau - x sst (eigenvalues of H-multiplication);
   (4) second row of the quantum differential equation for both idempotents: (D_i r_{1,i})[m-1] = L_j (r_{1,i}[m] - L_i r_{0,i}[m]),
       D_i = d/dt - (d/dt Delta_i)/(2 Delta_i);
   (5) symplecticity: Sum_{p+q=n} (-1)^q <v_i[p], v_j[q]> = delta_{n0} delta_ij Delta_i (exact in Q);
   (6) E11 flow: d/dt (t^2 rho^2 E11(z, w)) = Theta(t) d/dt(t g(z)) d/dt(t g(w)), g(z) = G(-z), Theta = (L2 - L1)(c x L2 - x^2)^2/Delta_2;
   (7) leaf structure of the closed form far beyond the stored data: [t^k] (rho z r_{0,1}) is a combination of 1/(c + j z), j = 1..k+1, with no
       constant, no pole at z = 0 (k <= KX, exact in x).
   Usage: wolframscript -file hvc_leaf_proof_check.wls [KT] [KX]      Output: hvc_leaf_proof_check_out.txt *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine]; KT = If[Length[args] >= 1, ToExpression[args[[1]]], 6]; KX = If[Length[args] >= 2, ToExpression[args[[2]]], 10];
Get["hvc_lib.m"]; Get["hvc_factors_g4_Bm.mx"]; t0 = AbsoluteTime[];
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]);
ZM = ZMAX; c = 1 + x;
S0 = Together[SsumQ]; P0 = Together[PprodQ];
sq = ((1 + x)/x) Sqrt[1 + 4 x (1/Q - 1)/(1 + x)^2]; U1 = (S0 + sq)/2; U2 = (S0 - sq)/2;
NT = KT + 3;
tidy[s_SeriesData] := MapAt[Together, s, 3]; tidy[s_] := Together[s];
U1t = Normal[tidy[Series[U1 /. Q -> Exp[t], {t, 0, NT}]]];
ser[e_, n_] := tidy[Series[e /. Q -> Exp[t] /. U -> U1t, {t, 0, n}]];
evC[cc_C] := cc[[1]] + cc[[2]] U;
(* (1) Delta identity *)
spot = x Log[1 - Q U] - Log[U] + Log[1 - U]; W2 = Together[U D[U D[spot, U], U]];
genSym = U; Ssum = SsumQ; Pprod = PprodQ;
w2C = toC[x^2 W2];
emit["(1) Delta_1 = x^2 W2(U1) as an identity in the quadratic algebra (so Delta_2 = x^2 W2(U2) too): " <> ToString[Together[w2C[[1]] - DeltaCN[1][[1]]] === 0 && Together[w2C[[2]] - DeltaCN[1][[2]]] === 0]];
(* t-series of the stored data *)
tau = 1 - Exp[-t];
del1 = ser[evC[DeltaCN[1]], KT + 1]; del2 = ser[evC[DeltaCN[2]], KT + 1];
rho = tidy[Series[(Normal[ser[evC[DeltaCN[1]] t^2, KT + 2]]/(-x c^3))^(-1/2), {t, 0, KT + 1}]];
L1 = ser[-(Q x U + Q)/(1 - Q), KT + 1]; L2 = ser[-(Q x (S0 - U) + Q)/(1 - Q), KT + 1];
sstar = tidy[Series[(1 - U1t)/tau, {t, 0, KT + 2}]];
emit["(3) L2 = x s*, L1 = c/tau - x s*: " <> ToString[{Normal[tidy[L2 - x sstar]] === 0, Normal[tidy[L1 - (c/tau - x sstar)]] === 0}]];
(* (3b) critical equation x tau s^2 - c s + 1 = 0, symmetric functions of the eigenvalues, metric in terms of L *)
emit["(3b) x tau s*^2 - c s* + 1 = 0: " <> ToString[Normal[tidy[Series[x tau Normal[sstar]^2 - c Normal[sstar] + 1, {t, 0, KT + 1}]]] === 0] <>
   "; L1 + L2 = c/tau, L1 L2 = x/tau: " <> ToString[{Normal[tidy[L1 + L2 - Series[c/tau, {t, 0, KT + 1}]]] === 0, Normal[tidy[L1 L2 - Series[x/tau, {t, 0, KT + 1}]]] === 0}] <>
   "; Delta_1 = x^3 (L1 - L2)/(x - c L2), Delta_2 = x^3 (L2 - L1)/(x - c L1): " <>
   ToString[{Normal[tidy[del1 - x^3 (L1 - L2)/(x - c L2)]] === 0, Normal[tidy[del2 - x^3 (L2 - L1)/(x - c L1)]] === 0}]];
(* (2) closed form of the singular row *)
phiT[s_] := Log[s] + x Log[1 - s] - Log[1 - tau s];
dphi = Series[phiT[sstar] - (Log[1/c] + x Log[x/c]), {t, 0, KT + 1}];
dphi = MapAt[Simplify[PowerExpand[#], x > 0] &, dphi, 3];
emit["dphi(t) = " <> ToString[InputForm[Normal[dphi]]]];
F21 = Sum[(tau/z)^k Product[(1 + j z)^2, {j, 1, k}]/(k! Product[c + j z, {j, 2, k + 1}]), {k, 0, KT + 1}];
Gz = tidy[Series[(c/(c + z)) (tau/t) F21 Exp[-Normal[dphi]/z], {t, 0, KT}]];
Gk = Table[Together[SeriesCoefficient[Gz, k]], {k, 0, KT}];
emit["(2a) pole at z = 0 absent in [t^k] G, k = 0.." <> ToString[KT] <> ": " <> ToString[Table[(Denominator[Gk[[k + 1]]] /. z -> 0) =!= 0, {k, 0, KT}]]];
r01 = Table[ser[evC[ratCN[0, 1][m]], KT], {m, 0, ZM}];
lhs = Table[Together[SeriesCoefficient[rho r01[[m + 1]], k]], {k, 0, KT}, {m, 0, ZM}];
rhs = Table[Together[SeriesCoefficient[Series[Gk[[k + 1]], {z, 0, ZM}], m]], {k, 0, KT}, {m, 0, ZM}];
emit["(2b) rho r_{0,1}(z; t) = closed form, all t^k (k <= " <> ToString[KT] <> "), z^m (m <= " <> ToString[ZM] <> "): " <> ToString[Union[Flatten[Together[lhs - rhs]]] === {0}]];
alg = tidy[Series[(1 - tau sstar) Sqrt[(-1/sstar^2 - x/(1 - sstar)^2 + tau^2/(1 - tau sstar)^2)/(-c^3/x)], {t, 0, KT + 1}]];
emit["(2c) rho = tau/(t alg): " <> ToString[Normal[tidy[rho - Series[tau/(t Normal[alg]), {t, 0, KT + 1}]]] === 0]];
Print["time ", Round[AbsoluteTime[] - t0], " s"];
(* (4) second QDE row, both idempotents; rows of root 2 are the stored ratCN[a, 2] evaluated at U1 (conjugate encoding) *)
rows[i_] := rows[i] = Table[{ser[evC[ratCN[0, i][m]], KT], ser[evC[ratCN[1, i][m]], KT]}, {m, 0, ZM}];
Dop[f_, dl_] := tidy[D[f, t] - f D[dl, t]/(2 dl)];
qde[i_] := Module[{dl = If[i == 1, del1, del2], Li = If[i == 1, L1, L2], Lj = If[i == 1, L2, L1], rr = rows[i]},
   Table[Normal[tidy[Dop[rr[[m, 2]], dl] - Lj (rr[[m + 1, 2]] - Li rr[[m + 1, 1]])]] === 0, {m, 1, ZM}]];
emit["(4) second QDE row, idempotent 1 (m = 1..9): " <> ToString[qde[1]] <> "; idempotent 2: " <> ToString[qde[2]]];
(* (5) symplecticity, exact in Q in the algebra *)
pairC[i_, j_, p_, q_] := cadd @@ Flatten[Table[cscale[etainvCur[[a + 1, b + 1]], cmul[ratCN[a, i][p], ratCN[b, j][q]]], {a, 0, 1}, {b, 0, 1}]];
symp = Table[Module[{s = cadd @@ Table[cscale[(-1)^q, pairC[i, j, n - q, q]], {q, 0, n}], target},
     target = If[n == 0 && i == j, DeltaCN[i], C[0, 0]];
     Together[s[[1]] - target[[1]]] === 0 && Together[s[[2]] - target[[2]]] === 0], {i, 1, 2}, {j, 1, 2}, {n, 0, ZM}];
emit["(5) symplecticity <v_i(w), v_j(-w)> = delta_ij Delta_i, exact in Q, w^n for n <= 9: " <> ToString[Union[Flatten[symp]]]];
Print["time ", Round[AbsoluteTime[] - t0], " s"];
(* (6) E11 flow *)
gser = Sum[(-z)^m SeriesCoefficient[Series[Normal[Gz], {z, 0, ZM}], m], {m, 0, ZM}];   (* g(z) = G(-z), t-polynomial to KT, z-polynomial to z^9 *)
Theta = tidy[(L2 - L1) (c x L2 - x^2)^2/del2];
emit["Theta(t) = " <> ToString[InputForm[Factor /@ Normal[Theta]]]];
dtg = D[t gser, t];
flowR = Expand[Normal[Theta] dtg (dtg /. z -> w)];
E11s = Table[If[m + n <= ZM - 1, ser[evC[EdN[1, 1][m, n]], KT + 1], 0], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
calE = Normal[tidy[Series[t^2 Normal[rho]^2, {t, 0, KT + 2}]]];
flowOK = Table[If[m + n <= ZM - 1,
     Module[{lhsS = tidy[Series[D[Normal[tidy[Series[calE Normal[E11s[[m + 1, n + 1]]], {t, 0, KT + 1}]]], t], {t, 0, KT}]],
        rhsS = tidy[Series[Coefficient[Coefficient[flowR, z, m], w, n], {t, 0, KT}]]},
       Normal[tidy[lhsS - rhsS]] === 0], True], {m, 0, ZM - 1}, {n, 0, ZM - 1}];
emit["(6) E11 flow d/dt(t^2 rho^2 E11) = Theta d/dt(t g(z)) d/dt(t g(w)) at every (m, n) with m + n <= 8, t^0..t^" <> ToString[KT] <> ": " <> ToString[And @@ Flatten[flowOK]]];
Print["time ", Round[AbsoluteTime[] - t0], " s"];
(* (7) leaf structure of the closed form for k <= KX (exact in x, no stored data) *)
tauS = Normal[Series[1 - Exp[-t], {t, 0, KX + 1}]];
U1x = Normal[tidy[Series[U1 /. Q -> Exp[t], {t, 0, KX + 3}]]];
sstarX = Normal[tidy[Series[(1 - U1x)/tauS, {t, 0, KX + 1}]]];
dphiX = MapAt[Simplify[PowerExpand[#], x > 0] &, Series[(Log[s] + x Log[1 - s] - Log[1 - tauS s] /. s -> sstarX) - (Log[1/c] + x Log[x/c]), {t, 0, KX + 1}], 3];
F21X = Sum[(tauS/z)^k Product[(1 + j z)^2, {j, 1, k}]/(k! Product[c + j z, {j, 2, k + 1}]), {k, 0, KX + 1}];
GzX = tidy[Series[(c z/(c + z)) (tauS/t) F21X Exp[-Normal[dphiX]/z], {t, 0, KX}]];
leafX = Table[Module[{f = Together[SeriesCoefficient[GzX, k]], den, ok0, okInf, pf},
    den = Denominator[f];
    ok0 = (den /. z -> 0) =!= 0;
    okInf = Together[Limit[f, z -> Infinity] - If[k == 0, c, 0]] === 0;
    pf = Cases[FactorList[den], {p_, e_} /; ! FreeQ[p, z] :> {Together[p/Coefficient[p, z, 0]], e}];
    {k, ok0, okInf, pf}], {k, 0, KX}];
emit["(7) closed form, [t^k] rho z r_{0,1}: {k, regular at z = 0, limit at infinity = c delta_k0, z-factors of the denominator (normalised, multiplicity)}:"];
Do[emit["   " <> ToString[InputForm[l]]], {l, leafX}];
emit["total time " <> ToString[Round[AbsoluteTime[] - t0]] <> " s"];
Export["hvc_leaf_proof_check_out.txt", StringRiffle[out, "\n"], "Text"];
