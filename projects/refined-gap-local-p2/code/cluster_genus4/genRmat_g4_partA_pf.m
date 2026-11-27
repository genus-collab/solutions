(* Part A1 (numeric-sampling variant): v15's PF ansatz on the delta-curve family
     ee1 = 3/(1-u), ee2 = 3/(1-u), ee3 = 1/(1-u)  (alpha = 1),
   solved at rational sample values u_s (all lower levels numeric at the same u_s), then each pp[k][l][m](u)
   is reconstructed as P(u)/(u^a (1-u)^b) by polynomial interpolation with verification on extra samples.
   Usage: wolframscript -file genRmat_g4_partA_pf.m KMAX   (this file: the Part A1 cell of localp2final_CY5_v16.nb, unchanged except for the argument handling and the portable PFu<k>.m checkpoints)
   Output: PFu<k>.mx (PFu rules, exact rational functions of u) after each level k. *)

KMAX = ToExpression[$ScriptCommandLine[[2]]];
SetDirectory[DirectoryName[$InputFileName]];
log[m_] := Module[{s = OpenAppend["partA_pf_num.log"]}, WriteString[s, DateString[{"Hour", ":", "Minute", ":", "Second"}], "  ", m, "\n"]; Close[s]; Print[m]];

ee1 = 3/(1 - u); ee2 = 3/(1 - u); ee3 = 1/(1 - u);
polPFL[yy_, X_] := ee1 X^2 + ee2 X + ee3 + (1 + 27 yy) X^3;
ySol[X_] := -(ee1 X^2 + ee2 X + ee3 + X^3)/(27 X^3);
Module[{Pf, r1, r2, r3, Xf},
  Pf[a_, b_] := polPFL[a, b];
  r1[x_, X_] = Together[-D[Pf[x, X], x]/D[Pf[x, X], X]];
  r2[x_, X_] = Together[D[r1[x, Xf[x]], x] /. Xf'[x] -> r1[x, Xf[x]] /. Xf[x] -> X];
  r3[x_, X_] = Together[D[r2[x, Xf[x]], x] /. Xf'[x] -> r1[x, Xf[x]] /. Xf[x] -> X];
  Xd1[x_, X_] = r1[x, X]; Xd2[x_, X_] = r2[x, X]; Xd3[x_, X_] = r3[x, X];];

ZMAX = KMAX;
IAs[y_] = Exp[mm[y]/z] Sum[S0[k][y] z^k, {k, 0, ZMAX + 1}];
pf = Expand[Exp[-mm[y]/z] Block[{coeff = {ee3, ee2 z, ee1 z^2, z^3},
      coeff2 = CoefficientList[Series[(3 z DY) (3 z DY + z) (3 z DY + 2 z), {DY, 0, 3}], DY]},
     Sum[coeff[[n + 1]] D[IAs[Exp[Y]], {Y, n}], {n, 0, 3}] + y Sum[coeff2[[n + 1]] D[IAs[Exp[Y]], {Y, n}], {n, 0, 3}]] /. Y -> Log[y]
   /. {mm'[y] -> XX[y]/y, mm''[y] -> D[XX[y]/y, y], mm'''[y] -> D[XX[y]/y, {y, 2}]}];
pfk[k_] := pfk[k] = Coefficient[pf, z, k];
log["pf built"];

S0ans[k_][yy_, X_] := Sum[pp[k][l][m] X^l yy^m, {l, 0, 2}, {m, 0, 2 k}]/(3 ee3 + ee2 X)^(k + 2)/((1 + 27 yy) (-ee2^3 + 27 (1 + 27 yy) ee3^2))^k;
Sder[k_, n_] := Sder[k, n] = Module[{f = S0ans[k][yy, Xf[yy]]}, D[f, {yy, n}] /. {Derivative[3][Xf][yy] -> X3, Derivative[2][Xf][yy] -> X2, Derivative[1][Xf][yy] -> X1, Xf[yy] -> X0}];
subsAnsatz[k_] := subsAnsatz[k] = Module[{e = pfk[k]},
   e = e /. {Derivative[n_][S0[j_]][y] :> Sder[j, n], S0[j_][y] :> Sder[j, 0]};
   e /. {Derivative[n_][XX][y] :> {X1, X2, X3}[[n]], XX[y] -> X0, y -> yy}];

samplesX = Table[-2/3 - 1/(7 + 3 s), {s, 0, 80}];
usamples = Table[1/(4 + 3 s) + 1/(100 + 7 s), {s, 0, 140}];    (* rational u in (0,1) *)
NUS[k_] := Max[24, 2 (7 k - 1) + 8];  (* initial number of u samples per level; increased if reconstruction fails *)
known = <||>;
(* reuse the highest available PFu<k>.mx (each contains all levels <= k) *)
Module[{avail = Select[Range[0, KMAX], FileExistsQ["PFu" <> ToString[#] <> ".m"] &]},   (* portable checkpoints PFu<k>.m (PFu6.m is provided) *)
  If[avail =!= {}, PFu = Get["PFu" <> ToString[Max[avail]] <> ".m"];
   Do[known[k] = Select[PFu, #[[1, 0, 0, 1]] == k &], {k, 0, Max[avail]}]]];
Print["levels loaded from disk: ", Keys[known]];

(* numeric values: val[k][u_s] = association pp[k][l][m] -> rational *)
val[k_][us_] := val[k][us] = If[KeyExistsQ[known, k], Thread[known[k][[All, 1]] -> (known[k][[All, 2]] /. u -> us)], Module[{unknowns, nU, e, eqs = {}, sol, lower},
   unknowns = Flatten[Table[pp[k][l][m], {l, 0, 2}, {m, 0, 2 k}]]; nU = Length[unknowns];
   lower = Flatten[Table[val[j][us], {j, 0, k - 1}]];
   e = subsAnsatz[k + 1] /. u -> us /. lower;
   Do[Module[{Xs = samplesX[[s]], ys, vals}, ys = ySol[Xs] /. u -> us;
     vals = {yy -> ys, X0 -> Xs, X1 -> Xd1[ys, Xs] /. u -> us, X2 -> Xd2[ys, Xs] /. u -> us, X3 -> Xd3[ys, Xs] /. u -> us};
     AppendTo[eqs, e /. vals]], {s, 1, nU + 3}];
   sol = Solve[Join[Thread[eqs == 0], {pp[k][0][0] == KroneckerDelta[k, 0]}], unknowns];
   If[Length[sol] != 1 || Length[Flatten[sol]] != nU, log["level " <> ToString[k] <> " u=" <> ToString[us] <> ": bad solve " <> ToString[Length[sol]]]; Abort[]];
   Flatten[sol]]];

(* rational reconstruction from exact samples {u_s, v_s}: N(u)/D(u) with deg N + deg D + 1 = n - 5 (balanced),
   one NullSpace computation; verified on the 5 spare points *)
reconstruct[pts_] := Module[{n = Length[pts], us = pts[[All, 1]], vs = pts[[All, 2]], m, dd, dn, A, ns, v, N0, D0, ok},
   If[Union[vs] === {0}, Return[0]];
   m = n - 5; dd = Floor[(m - 1)/2]; dn = m - 1 - dd;
   A = Table[Join[Table[us[[s]]^j, {j, 0, dn}], Table[-vs[[s]] us[[s]]^j, {j, 0, dd}]], {s, 1, m}];
   ns = NullSpace[A]; If[ns === {}, Return[$Failed]];
   v = First[ns]; N0 = Sum[v[[j + 1]] u^j, {j, 0, dn}]; D0 = Sum[v[[dn + 2 + j]] u^j, {j, 0, dd}];
   If[D0 === 0, Return[$Failed]];
   ok = And @@ Table[((N0 /. u -> us[[s]]) - vs[[s]] (D0 /. u -> us[[s]])) == 0, {s, m + 1, n}];
   If[ok, Together[N0/D0], $Failed]];

PFu = {};
Do[
  Module[{t0 = AbsoluteTime[], keys, nus = NUS[k], recon, done = False},
   keys = Flatten[Table[pp[k][l][m], {l, 0, 2}, {m, 0, 2 k}]];
   If[KeyExistsQ[known, k], PFu = Join[PFu, known[k]]; log["level " <> ToString[k] <> " loaded from disk"]; Continue[]];
   While[! done,
    Do[val[k][usamples[[s]]], {s, 1, nus}];
    recon = Table[key -> reconstruct[Table[{usamples[[s]], key /. val[k][usamples[[s]]]}, {s, 1, nus}]], {key, keys}];
    If[FreeQ[recon, $Failed], done = True,
     log["level " <> ToString[k] <> ": reconstruction failed with " <> ToString[nus] <> " u-samples; increasing"]; nus += 8;
     If[nus > Length[usamples], log["too many samples needed -- ABORT"]; Abort[]]]];
   PFu = Join[PFu, recon];
   Put[PFu, "PFu" <> ToString[k] <> ".m"];
   log["level " <> ToString[k] <> " done: " <> ToString[Length[keys]] <> " coefficients, " <> ToString[nus] <> " u-samples, " <>
     ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s; max LeafCount " <> ToString[Max[LeafCount /@ recon[[All, 2]]]]];
   ], {k, 0, KMAX}];
log["saved PFu" <> ToString[KMAX] <> ".m"];

If[FileExistsQ["v15data.mx"],
  Module[{ok, save = PFu},
   Get["v15data.mx"]; PFu = save;
   ok = Table[Module[{v15v = r[[2]] /. {ee[2] -> 3/(1 - u), ee[3] -> 1/(1 - u)}}, Simplify[(r[[1]] /. PFu) - v15v] === 0],
     {r, Select[PFAsRepl, MemberQ[PFu[[All, 1]], #[[1]]] &]}];
   log["agreement with v15 PFAsRepl on the curve: " <> ToString[Union[ok]] <> " (" <> ToString[Length[ok]] <> " coefficients)"]]];