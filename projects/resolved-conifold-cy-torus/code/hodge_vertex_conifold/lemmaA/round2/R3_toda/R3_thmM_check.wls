(* R3_thmM_check.wls -- independent symbolic (in tau) check of Theorem M (leading order at the dilaton point of the first
   modified-KP equation of the 2D Toda hierarchy for the Marino-Vafa tau function), for one-part targets R = (j).
   Data: lemmaA/data/q_g<g>_r<r>.txt (Mumford-formula Hodge tables, c-normalisation, x symbolic); genus 0 and 1 closed forms.
   No topological recursion anywhere.  Run from the project folder hodge_vertex_conifold:
     wolframscript -file lemmaA/round2/R3_toda/R3_thmM_check.wls
   Expected output: every line "... 0" and a final "ALL ZERO: True". *)
dataDir = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "data"}];
readQ[g_, r_] := Module[{t = Import[FileNameJoin[{dataDir, "q_g" <> ToString[g] <> "_r" <> ToString[r] <> ".txt"}], "Text"], p},
   p = StringSplit[t, " | ", 3]; ToExpression[p[[3]]]];
Q[2, 1] = readQ[2, 1]; Q[2, 2] = readQ[2, 2]; Q[3, 1] = readQ[3, 1]; Q[3, 2] = readQ[3, 2]; Q[4, 1] = readQ[4, 1]; Q[4, 2] = readQ[4, 2];
sig = tau (1 + tau);
Pi2[k_] := Product[k tau + a, {a, 1, k - 1}]/(k - 1)!;
(* X_g(T), c = 1 normalisation: Xi = c^(3g-3+2r) Q(s=-1), c = 1+x = 1/(1+tau), x = -tau/(1+tau) *)
X[0, T_List] := Module[{r = Length[T], J = Total[T]}, Which[r == 1, 1/(J (J + 1)), r == 2, 1/J, True, Product[J - i, {i, 1, r - 3}]]];
X[1, {j_}] := (j - 2)/24 - 1/(24 sig);
X[1, {j1_, j2_}] := 1/12 - (j1 + j2)/8 + (j1^2 + j1 j2 + j2^2)/24 - (j1 + j2 - 1)/(24 sig);
X[g_ /; g >= 2, T_List] := Module[{r = Length[T], e},
   e = Q[g, r] /. Join[{s -> -1, x -> -tau/(1 + tau)}, Table[ToExpression["j" <> ToString[i]] -> T[[i]], {i, r}]];
   Together[(1/(1 + tau))^(3 g - 3 + 2 r) e]];
ff[a_, k_] := Product[a - i, {i, 0, k - 1}];
dl[g_, T_List] := Total[T] - Length[T] - 2 g + 2;
Phi[g_, k_, T_List] := If[T === {}, Which[g > 0, 0, k == 1, tau, k == 2, tau^2/2, True, 0],
   If[ff[dl[g, T], k] == 0, 0, sig^g tau^k Total[T] ff[dl[g, T], k] X[g, T]/k!]];
(* Theorem M for R = (j): x-power n = 2G + eps *)
M[G_, eps_, j_] := Module[{n = 2 G + eps, t = 0, R2 = {2, j}, rest},
   Do[With[{k = n + 1 - 2 g}, If[k >= 1 && ff[dl[g, R2], k] != 0, t += -2 Pi2[2] sig^g tau^k ff[dl[g, R2], k] X[g, R2]/k!]], {g, 0, G + 1}];
   If[eps == 0, t += 2 sig^(G + 1) j (j + 1) X[G, {j}]];
   Do[With[{k = n - 2 g}, If[k >= 1 && ff[dl[g, {j}], k] != 0, t += sig^(g + 1) tau^k j (j + 1) ff[dl[g, {j}], k] X[g, {j}]/k!]], {g, 0, G}];
   Do[Do[With[{rem = n + 2 - 2 g1 - k1}, Do[With[{k2 = rem - 2 g2}, If[k2 >= 1, t += Phi[g1, k1, sp[[1]]] Phi[g2, k2, sp[[2]]]]], {g2, 0, Floor[(rem - 1)/2]}]],
       {g1, 0, Floor[(n + 1)/2]}, {k1, 1, n + 1 - 2 g1}], {sp, {{{j}, {}}, {{}, {j}}}}];
   Together[t]];
res = Table[With[{v = M[G, eps, j]}, Print["G=", G, " eps=", eps, " R=(", j, "): ", v]; v], {G, 2, 4}, {eps, 0, 1}, {j, 2, 12}];
Print["ALL ZERO: ", AllTrue[Flatten[res], # === 0 &]];
(* control: perturb X[3,{2,5}] -> must give non-zero for (G,R) = (3,(5)) *)
Xsave = X[3, {2, 5}]; X[3, {2, 5}] = Xsave + 1/1000;
Print["control (perturbed Xi_3(2,5)): ", M[3, 0, 5]];
