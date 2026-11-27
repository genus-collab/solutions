(* E8 exact checks: eta^{-1} bivector vs matrix, hatted unit, Psi_Q unit identity, dilaton binomial *)
z5 = Exp[2 Pi I/5];
J = Table[If[Mod[a + b, 5] == 3, 1, 0], {a, 0, 4}, {b, 0, 4}];
eta = 5 J; etaInv = Inverse[eta];
Print["etaInv == eta/25: ", etaInv == eta/25];
nu = {0, -1, 3, 2, 1}; eps = {1, -1, -1, 1, 1};
Lam[s_] := DiagonalMatrix[s^(1 - nu)];
Print["Lambda_s etaInv Lambda_s == etaInv: ", Simplify[Lam[s].etaInv.Lam[s] == etaInv]];
Psi0 = Table[eps[[b + 1]] z5^(i (b + 1)), {b, 0, 4}, {i, 0, 4}];
Print["Psi0^T etaInv Psi0 == 1: ", FullSimplify[Transpose[Psi0].etaInv.Psi0 == IdentityMatrix[5]]];
(* Psi_Q unit identity in N: m_a, Delta_j^{-1/2} = zeta^j L/C0 *)
mm = {L/C0, L^2/(C0 C1), C0 C1/L^2, C0/L, 1};
PsiQ = Table[z5^(j (a + 1)) mm[[a + 1]], {a, 0, 4}, {j, 0, 4}];
unitQ = Table[z5^j L/C0, {j, 0, 4}];
Print["PsiQ . 1_TFT (lowered) = ", FullSimplify[PsiQ.unitQ]];
Print["PsiQ^T etaInv PsiQ == 1: ", FullSimplify[Transpose[PsiQ].etaInv.PsiQ == IdentityMatrix[5]]];
(* symplectic R' = Psi X with X(z) = exp(zB), B symmetric rational; truncate at order K *)
K = 5;
B = {{1, 2, 0, -1, 3}, {2, -1, 1, 0, 1}, {0, 1, 2, 1, -2}, {-1, 0, 1, 0, 1}, {3, 1, -2, 1, 1}}/7;
Xs[t_] := Sum[MatrixPower[B, k] t^k/k!, {k, 0, K}];
Xinv[t_] := Sum[MatrixPower[-B, k] t^k/k!, {k, 0, K}];
(* isometry: Psi^m-like: Lambda_varsigma^{-1} Psi0 with varsigma = c (a number) *)
c = 3/2;
Psi = Inverse[Lam[c]].Psi0;
Print["Psi^T etaInv Psi == 1: ", FullSimplify[Transpose[Psi].etaInv.Psi == IdentityMatrix[5]]];
Rp[t_] := Psi.Xs[t];
RpInv[t_] := Xinv[t].Inverse[Psi];
trunc[e_] := Normal[Series[e /. {p -> x p, q -> x q}, {x, 0, K - 1}]] /. x -> 1;
Nbiv = trunc[IdentityMatrix[5] - RpInv[p].eta.Transpose[RpInv[q]]];
NX = trunc[IdentityMatrix[5] - Xinv[p].Transpose[Xinv[q]]];
Print["N (bivector reading, lowered matrix eta) == N_X up to degree K-1: ", FullSimplify[Nbiv - NX] == ConstantArray[0, {5, 5}]];
Nmat = trunc[IdentityMatrix[5] - RpInv[p].etaInv.Transpose[RpInv[q]]];
Print["N(0,0) misread: ", FullSimplify[Nmat /. {p -> 0, q -> 0}][[1, 1]]];
anti[NN_] := NN - (1/2) (NN /. {p -> -q}) - (1/2) (NN /. {q -> -p});
Ebiv = Simplify[Cancel[Together[anti[Nbiv]/(p + q)]]];
Emat = Simplify[Cancel[Together[anti[Nmat]/(p + q)]]];
lowdeg[e_] := Normal[Series[e /. {p -> x p, q -> x q}, {x, 0, K - 2}]] /. x -> 1;
Print["E_misread == E_true/25 (deg <= K-2): ", FullSimplify[lowdeg[Emat - Ebiv/25]] == ConstantArray[0, {5, 5}]];
Print["matrix form rmk:flh-fail(iii) agrees: ", FullSimplify[trunc[IdentityMatrix[5] - Transpose[Rp[-p]].etaInv.Rp[-q]] - Nbiv] == ConstantArray[0, {5, 5}]];
(* hatted unit: reading A (same map, hatted lowered coords) vs reading B (Lambda_s as endomorphism, hat1 = s 1 as vector) *)
one = {0, 0, 0, 5, 0}; (* lowered components of H^0 = 1 *)
qR[t_] := Xinv[t].Inverse[Psi].(-t one);
hatRinvA[t_] := RpInv[s^4 t].Inverse[Lam[s]]; (* inverse of hat R' = Lam_s R'(s^4 zh) *)
hat1A = Lam[s].(s one); (* hatted lowered comps of the vector s*1 *)
Print["hatted lowered comps of s*1: ", Simplify[hat1A]];
Print["reading A: q_hat - s^-3 q(s^4 zh) = 0: ", Simplify[trunc[hatRinvA[w].(-w hat1A) - s^-3 qR[s^4 w]] /. p -> 0] == {0, 0, 0, 0, 0}];
Print["reading B (insert flat s*1): q_hat / q(s^4 zh) ratio: ", Simplify[(hatRinvA[w].(-w s one))[[1]]/(qR[s^4 w])[[1]]]];
(* dilaton binomial: sum_b (N)_b/b! e^b = (1-e)^{-N} *)
Print["binomial: ", And @@ Table[Simplify[Normal[Series[(1 - e)^(-NN), {e, 0, 8}]] - Sum[Pochhammer[NN, b]/b! e^b, {b, 0, 8}]] == 0, {NN, 1, 6}]];
