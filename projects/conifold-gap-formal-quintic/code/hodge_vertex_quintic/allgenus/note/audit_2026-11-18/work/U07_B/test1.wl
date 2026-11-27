(* U07_B sanity tests: tower values, frames, Lemma E / Cor E' formula, Lemma T(b),(d) for a random graded symplectic P *)
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
M = 9;        (* work modulo v^M *)
NN = 10;      (* highest z^{-n} needed mod v^M *)
ZK = 18;      (* truncation order of P in z *)
trv[e_] := Normal[Series[e, {v, 0, M - 1}]];
(* vanishing solution YC_b, b=0..4, as polynomial in zi = 1/z *)
YC[b_] := Sum[If[Mod[n - (nu[[b + 1]] + 1), 5] == 0,
    eps[[b + 1]] (-1)^n v^(n - Floor[n/5])/(4^Floor[n/5] Floor[n/5]! (n - Floor[n/5])!) z^(-n), 0], {n, 0, NN}];
YCvec = Table[YC[b], {b, 0, 4}];
(* check against the residue definition IC_m *)
IC[m_] := Sum[If[m + 4 i + 1 >= 0, (-1/(4 z))^i/i! (-v/z)^(m + 4 i + 1)/(m + 4 i + 1)!, 0], {i, 0, 4}];
chk0 = Table[trv[Expand[YC[b] - eps[[b + 1]] IC[nu[[b + 1]]]] /. z -> 1/zi] , {b, 0, 4}];
Print["YC closed form vs residue def (mod v^M, should be 0 up to z^-NN): ",
  Table[PolynomialMod[Sum[Coefficient[chk0[[b]], zi, n] zi^n, {n, 0, NN}], v^M], {b, 1, 5}]];
eta = Table[If[Mod[a + b - 3, 5] == 0, 5, 0], {a, 0, 4}, {b, 0, 4}];
etainv = Inverse[eta];
Print["eta^-1 == eta/25: ", etainv == eta/25];
(* Pi_i on q = list of 5 polynomials in z *)
coefz[e_, k_] := Coefficient[Expand[e /. z -> 1/zi] /. zi -> 1/z0, z0, k] ;
lc[e_, k_] := Module[{ee = Expand[e]}, If[k == 0, ee /. z -> 0 /. (zz_^p_ /; False) -> 0, Coefficient[ee, z, k]]];
(* robust Laurent coefficient extraction *)
LCoef[e_, k_] := Module[{ee = Expand[e z^40]}, Coefficient[ee, z, k + 40]];
PiT[i_, q_] := trv[-LCoef[Sum[(YCvec[[b]] /. z -> -z) q[[b]], {b, 1, 5}], i]];
H[a_] := UnitVector[5, a + 1];
Print["Pi0(H)=", PiT[0, H[1]], "  Pi0(z 1)=", PiT[0, z H[0]], "  Pi0(z^6 1)=", PiT[0, z^6 H[0]]];
Print["Pi0(z^5 H)=", PiT[0, z^5 H[1]], " expected ", v^4/(4 24)];
(* frames: solve for f_j *)
frame[k_, a_] := Module[{f, q, eqs, sol},
   q = z^k H[a] + Sum[f[j] z^j, {j, 0, k - 1}] H[1];
   eqs = Table[PiT[i, q] == 0, {i, 0, k}];
   sol = Solve[eqs, Table[f[j], {j, 0, k - 1}]];
   trv[q /. First[sol]]];
Print["frame(1,0)=", frame[1, 0]]; Print["frame(2,4)=", frame[2, 4]]; Print["frame(3,3)=", frame[3, 3]];
Print["frame(4,2)=", frame[4, 2]]; Print["frame(6,0)=", frame[6, 0]];
trivial = Select[Flatten[Table[{k, a}, {a, {0, 2, 3, 4}}, {k, 0, 6}], 1], (frame[#[[1]], #[[2]]] === z^#[[1]] H[#[[2]]]) &];
Print["trivial frames: ", trivial];
(* random graded symplectic P *)
SeedRandom[20261117];
E5[a_, b_] := Table[If[i == Mod[a, 5] + 1 && j == Mod[b, 5] + 1, 1, 0], {i, 5}, {j, 5}];
amat[k_] := Module[{r, rs},
   r = Table[RandomInteger[{-3, 3}]/RandomInteger[{1, 3}], {b, 0, 4}];
   rs = Table[(r[[b + 1]] - (-1)^k r[[Mod[3 - b - k, 5] + 1]])/2, {b, 0, 4}];
   Sum[rs[[b + 1]] E5[b + k, b], {b, 0, 4}]];
ks = {1, 2, 3, 4, 5, 7};
A = Sum[z^k amat[k], {k, ks}];
tz[m_] := Map[Normal[Series[#, {z, 0, ZK}]] &, m, {2}];
expA = IdentityMatrix[5]; term = IdentityMatrix[5];
Do[term = tz[term . A/j]; expA = expA + term, {j, 1, ZK}];
P0 = DiagonalMatrix[{2, 3, 1/3, 1/2, -1}];
P = tz[P0 . expA];
symp = tz[Transpose[P /. z -> -z] . eta . P] - eta;
Print["(Sp) defect mod z^(ZK+1): ", Union[Flatten[Expand[symp]]]];
(* graded check: P_k has entries only at (b+k,b) *)
gradedQ = And @@ Flatten[Table[If[Mod[i - j - k, 5] != 0, Coefficient[P[[i + 1, j + 1]], z, k] === 0, True], {i, 0, 4}, {j, 0, 4}, {k, 0, ZK}]];
Print["graded: ", gradedQ];
Pvec = etainv . P . eta; P0vec = etainv . P0 . eta; d = P0vec[[1, 1]];
Print["d = ", d, "  (P0)_33 = ", P0[[4, 4]]];
(* inverse of Pvec mod z^(ZK+1) *)
Rvec = tz[Inverse[P0vec] . Pvec];
X = IdentityMatrix[5] - Rvec; Rinv = IdentityMatrix[5]; term = IdentityMatrix[5];
Do[term = tz[term . X]; Rinv = Rinv + term, {j, 1, ZK}];
Pvecinv = tz[Rinv . Inverse[P0vec]];
Print["inverse check: ", Union[Flatten[tz[Pvec . Pvecinv] - IdentityMatrix[5]]]];
(* edge E_P via numerator, using R^vec *)
Nnum = Expand[etainv - (Rinv /. z -> x) . etainv . Transpose[Rinv /. z -> y]];
(* coefficient E_{k,l} = sum_{j=0}^l (-1)^j N_{k+1+j, l-j} *)
Ncoef[k_, l_] := Map[Coefficient[Coefficient[#, x, k], y, l] &, Nnum, {2}];
Ecoef[k_, l_] := Sum[(-1)^j Ncoef[k + 1 + j, l - j], {j, 0, l}];
Lmax = 3; Kmax = NN + 3;
Etab = Table[Ecoef[k, l], {k, 0, Kmax}, {l, 0, Lmax}];
(* same edge with Pvec instead of Rvec *)
Nnum2 = Expand[etainv - (Pvecinv /. z -> x) . etainv . Transpose[Pvecinv /. z -> y]];
Print["edge with P^vec == edge with R^vec (numerators, low degree): ",
  Union[Flatten[Table[Map[Coefficient[Coefficient[#, x, k], y, l] &, Nnum - Nnum2, {2}], {k, 0, 8}, {l, 0, 8}]]]];
(* symmetry E(psi,psi')^T = E(psi',psi) *)
Print["edge symmetry: ", Union[Flatten[Table[Transpose[Etab[[k + 1, l + 1]]] - Etab[[l + 1, k + 1]], {k, 0, Lmax}, {l, 0, Lmax}]]]];
(* (Pi_i (x) 1) E : vector valued polynomial in y *)
PiE[i_] := Table[Sum[y^l PiT[i, Table[Sum[z^k Etab[[k + 1, l + 1]][[a, b]], {k, 0, Kmax}], {a, 1, 5}]], {l, 0, Lmax}], {b, 1, 5}];
(* prediction from Lemma E / Cor E' : (Pi_{k-1} (x) 1)E_P = (-1)^k Rvec^{-1}(y) [ z^{-k} P0vec^{-1} etainv (Y')_{>=k} ](y) *)
Yp = Map[trv, Expand[P . YCvec]];
YpGE[k_] := Table[Sum[LCoef[Yp[[a]], m] z^m, {m, k, ZK}], {a, 1, 5}];
pred[k_] := Module[{w},
   w = Inverse[P0vec] . etainv . Expand[YpGE[k]/z^k];
   Map[trv, Normal[Series[(Rinv /. z -> y) . (w /. z -> y), {y, 0, Lmax}]]] (-1)^k];
Do[Print["CorE' formula, i=", k - 1, ": diff = ", Union[Flatten[CoefficientList[#, {v, y}] & /@ Expand[PiE[k - 1] - pred[k]]]],
   "   (nonzero test: ", (PiE[k - 1] /. {y -> 1, v -> 1}) =!= {0, 0, 0, 0, 0}, ")"], {k, 1, 3}];
(* Lemma T: q_T = Pvec^{-1}(-z 1) *)
qT = Expand[Pvecinv . (-z H[0])];
Y0 = Yp[[1]];
Do[Print["Pi_", i, "(q_T) = ", PiT[i, qT], "   (-1)^(i+1)[z^(i-1)]Y'_0 = ", trv[(-1)^(i + 1) LCoef[Y0, i - 1]]], {i, 0, 5}];
TP = trv[-LCoef[Y0, -1]];
Print["T_P - p00 v = ", Expand[TP - P0[[1, 1]] v], "   d*p00 = ", d P0[[1, 1]]];
(* leaf: q_T = -(z/d) 1 + T>=2 *)
Print["q_T + (z/d) 1 low orders: ", Table[Coefficient[qT + z/d H[0], z, j], {j, 0, 2}]];
