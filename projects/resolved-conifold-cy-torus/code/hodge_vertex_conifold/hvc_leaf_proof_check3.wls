(* hvc_leaf_proof_check3.wls -- symbolic checks requested by leaf_proof/leaf_structure_proof.md (s.12.3: M2-M7, M9-M11) and by its
   completeness critic (Q3). No stored factors are used; hvc_lib.m supplies wickDirect and hCoeffs. Conventions as in
   hvc_leaf_proof_check.wls: eps1 = x, eps2 = 1, c = 1 + x, t = Log Q, tau = 1 - e^-t, sst the critical point of
   phi_tau(s) = Log[s] + x Log[1 - s] - Log[1 - tau s] with sst(0) = 1/c, delta = c (1 - 4 x tau/c^2)^(1/2) = c - 2 x tau sst.
   (M5)  (3.16) identically in (s, tau, z); (3.19) for generic S, Phi; the eikonal equation at dphi_tau = (1 - x sst)/(1 - tau);
         P_tr = delta; R_tr = -x (c + 3 delta)/(delta (c + delta)) (modulo x tau s^2 - c s + 1).
   (M6)  Lemma 3.2, k = 2..20: tau^k W_k(1 - tau s) at Q = 1/(1 - tau) is regular at (s, tau) = (1/c, 0) with value (-1)^k phi_0^(k)(1/c).
   (Q3)  Sum_j j A_j(t) = -c^2 rho(t) (that is, Dil_1(1) = -1 at every t), through t^12.
   (M7)  Remark 4.9: tau^2 <PhiF(z), PhiF(-z)> = -x c (c^2 - z^2), PhiF = (F, (z (1 - tau) d/dtau + omega) F), omega = (c + z)/tau - (1 + z).
   (M4)  Proposition 3.19: the transport recursion, m <= 5, through tau^10, against [z^m](Exp[-dphi/z] F), and Phi_0 = 1/alg.
   (M3)  Proposition 3.15 through z^15: the code's Wick coefficients (wickDirect) at W_k -> phi_0^(k)(1/c) are [z^m] h_w(z) c/(c + z);
         Proposition 8.6(a): at W_k -> (U d/dU)^k W_0 at U = -1/x, W_0 = c Log[1 - U] - Log[U], they are [z^m] h_w(z)^-1.
   (M9)  (8.5): Gauss_Dop[L2] = 0 at z^-2..z^5 as Laurent series in tau, L2 = Lnorm[phi_tau, (1 - tau s)^-1; s~] computed in the
         coordinate sigma = tau s (where all Taylor data are analytic at tau = 0).
   (M2)  (3.17): U_1 Lnorm[phi_tau, (1 - tau s)^-1; sst] = h_w(z) c/(c + z) F Exp[-dphi/z] alg, the left side by the Gaussian (Wick)
         formula with amplitude in the coordinate s; tau^0..tau^8, z^0..z^9 at x = 2, 1/3, -5/2, 3/7, and symbolically in x to (tau^5, z^6).
   (M10) E_11 from the closed form alone: rho v_1 = (G, (L_1 + z/t) G + z dG/dt), the code's difference quotient (evalQ) for m + n <= 14,
         t^-2..t^8, against (7.5); (M11) the leaf decomposition of (7.5) (second line, hence (7.7)), (7.6), and Theorem 11.4; at x = 2, 3/7, -5/2.
   Lnorm is computed by lnormS: Exp[Rr_f/z] = Sum_n En[n] by the exponential recursion (a term v^p z^q has weight p/2 + q), times the
   amplitude, then the Gaussian moments v^(2n) -> (2n-1)!! (-z/f2)^n (Definition 3.4, (D1)).
   Usage: wolframscript -file hvc_leaf_proof_check3.wls     Output: hvc_leaf_proof_check3_out.txt (rewritten after every line) *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_lib.m"]; t0 = AbsoluteTime[];
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]; Export["hvc_leaf_proof_check3_out.txt", StringRiffle[out, "\n"], "Text"]);
tm[] := " [" <> ToString[Round[AbsoluteTime[] - t0]] <> " s]";
c = 1 + x;
tidyC[s_SeriesData] := MapAt[Collect[#, {vv, zz, z}, Together] &, s, 3];
tidyC[s_] := Collect[s, {vv, zz, z}, Together];
coefV[s_SeriesData, p_] := MapAt[Coefficient[#, vv, p] &, s, 3];
coefV[s_, p_] := Coefficient[s, vv, p];
sc[s_SeriesData, k_] := SeriesCoefficient[s, k]; sc[s_, k_] := If[k == 0, s, 0];   (* t^k coefficient of a series that may be a constant *)
(* Lnorm[f, am; u0] from ffk[k] = f^(k)(u0) (k >= 3), f2 = f''(u0), aa[a] = am^(a)(u0)/a! (numbers or tau-series); {[z^0], ..., [z^M]} *)
lnormS[ffk_, f2_, aa_, M_] := Module[{nE = 2 M, g, En, iF2pow},
   g[n_] := g[n] = tidyC[ffk[n + 2] vv^(n + 2)/((n + 2)! zz)];
   En[0] = 1;
   Do[En[n] = tidyC[(1/n) Sum[k g[k] En[n - k], {k, 1, n}]], {n, 1, nE}];
   iF2pow[0] = 1; Do[iF2pow[n] = tidyC[iF2pow[n - 1]/f2], {n, 1, 3 M + 1}];
   Table[Module[{Ke = tidyC[Sum[aa[a] vv^a En[2 m - a], {a, 0, 2 m}]], gs},
      gs = tidyC[Sum[If[OddQ[p], 0, coefV[Ke, p] (p - 1)!! (-zz)^(p/2) iF2pow[p/2]], {p, 0, 6 m}]];
      gs = tidyC[gs/zz^m];
      If[! FreeQ[Normal[gs], vv | zz], Print["lnormS: unexpected v or z dependence at z^", m]];
      gs], {m, 0, M}]];
withX[x0_, body_] := If[x0 === x, body[], Block[{x = x0}, body[]]];

(* ---------- (M5) ---------- *)
Module[{f, am0, Dam, kE, al, be, ga, lhs, rhs, res1, res2, gaussOp, eik, Ptr, Rtr, Ttr, crit, red, sT, sTT, delS, PtrS, RtrS},
  al = 1/z + 1; be = 1/z + 1; ga = c/z + 2;
  f = Log[s] + x Log[1 - s] - Log[1 - tau s]; am0 = 1/(1 - tau s);
  Dam[q_] := D[q, tau] + q D[f, tau]/z; kE = s (1 - s)/(1 - tau s)^2;
  lhs = tau (1 - tau) Dam[Dam[am0]] + (ga - (al + be + 1) tau) Dam[am0] - al be am0;
  rhs = -(1 + 1/z) (D[kE, s] + D[f, s] kE/z);
  res1 = Together[lhs - rhs] === 0;
  gaussOp[Y_] := tau (1 - tau) D[Y, {tau, 2}] + (ga - (al + be + 1) tau) D[Y, tau] - al be Y;
  eik[p_] := tau (1 - tau) p^2 + (c - 2 tau) p - 1;
  Ptr = 2 tau (1 - tau) SS'[tau] + c - 2 tau; Rtr = tau (1 - tau) SS''[tau] + (2 - 3 tau) SS'[tau] - 2;
  Ttr[Y_] := tau (1 - tau) D[Y, {tau, 2}] + (2 - 3 tau) D[Y, tau] - Y;
  res2 = Simplify[Exp[-SS[tau]/z] gaussOp[Exp[SS[tau]/z] PH[tau]] - (z^-2 eik[SS'[tau]] PH[tau] + z^-1 (Ptr PH'[tau] + Rtr PH[tau]) + Ttr[PH[tau]])] === 0;
  crit = x tau s^2 - c s + 1;
  red[e_] := Together[PolynomialRemainder[Numerator[Together[e]], crit, s]] === 0;
  sT = (1 - x s)/(1 - tau);                       (* dphi_tau along s = sst *)
  delS = c - 2 x tau s;                           (* delta at s = sst *)
  sTT = D[sT, tau] + D[sT, s] x s^2/delS;         (* d/dtau along s = sst, dsst/dtau = x sst^2/delta *)
  PtrS = 2 tau (1 - tau) sT + c - 2 tau; RtrS = tau (1 - tau) sTT + (2 - 3 tau) sT - 2;
  emit["(M5) amplitude identity (3.16), identically in (s, tau, z): " <> ToString[res1] <> "; (3.19) for generic S, Phi: " <> ToString[res2] <>
    "; eikonal (1 - tau) Eik(dphi_tau) = 0 mod x tau s^2 - c s + 1: " <> ToString[red[(1 - tau) eik[sT]]] <>
    "; P_tr = c - 2 x tau sst identically and (c - 2 x tau s)^2 = c^2 - 4 x tau mod the critical equation (so P_tr = delta): " <>
    ToString[{Together[PtrS - delS] === 0, red[delS^2 - (c^2 - 4 x tau)]}] <>
    "; R_tr = -x (c + 3 delta)/(delta (c + delta)) mod the critical equation: " <> ToString[red[RtrS + x (c + 3 delS)/(delS (c + delS))]] <> tm[]]];

(* ---------- (M6) ---------- *)
Module[{WW, ok = {}, e, den, val, target},
  WW[1] = Together[U D[x Log[1 - Q U] - Log[U] + Log[1 - U], U]];
  Do[WW[k] = Together[U D[WW[k - 1], U]], {k, 2, 20}];
  Do[e = Together[tau^k WW[k] /. U -> 1 - tau s /. Q -> 1/(1 - tau)];
    den = Denominator[e];
    val = Together[(e /. tau -> 0) /. s -> 1/c];
    target = (-1)^k (k - 1)! c^k ((-1)^(k - 1) - x^(1 - k));
    AppendTo[ok, {Together[(den /. tau -> 0) /. s -> 1/c] =!= 0, Together[val - target] === 0}], {k, 2, 20}];
  emit["(M6) Lemma 3.2, k = 2..20, {denominator of tau^k W_k(1 - tau s) non-zero at (s, tau) = (1/c, 0), value there = (-1)^k phi_0^(k)(1/c)}: " <>
    ToString[Union[ok]] <> tm[]]];

(* ---------- (Q3) ---------- *)
Module[{NQ = 12, tauS, U1e, sst, dphi, zetaS, hyp, Cj, ahat, rhoS, lhs, ok},
  tauS = Normal[Series[1 - Exp[-t], {t, 0, NQ + 3}]];
  U1e = ((x - 1)/x + ((1 + x)/x) Sqrt[1 + 4 x (Exp[-t] - 1)/(1 + x)^2])/2;
  sst = Normal[tidyC[Series[(1 - U1e)/(1 - Exp[-t]), {t, 0, NQ + 2}]]];
  dphi = Integrate[Normal[tidyC[Series[1 - x sst, {t, 0, NQ + 2}]]] /. t -> uu, {uu, 0, t}];
  zetaS = Normal[tidyC[Series[tauS Exp[dphi/c], {t, 0, NQ + 2}]]];
  hyp[j_] := Sum[Pochhammer[j x/c, m]^2/(Pochhammer[j, m] m!) tauS^m, {m, 0, NQ + 1}];
  Cj[j_] := c^2 Pochhammer[1 - j/c, j - 1]^2/(j! (j - 1)!);
  ahat[j_] := Normal[tidyC[Series[(-1)^j Cj[j] zetaS^j hyp[j], {t, 0, NQ + 1}]]];
  lhs = Expand[Sum[j ahat[j], {j, 1, NQ + 2}]/t];
  rhoS = tidyC[Series[(tauS/t) 2/((1 - uu)^(1/4) (1 + (1 - uu)^(1/2))) /. uu -> 4 x tauS/c^2, {t, 0, NQ}]];
  ok = Table[Together[Coefficient[lhs, t, k] + c^2 SeriesCoefficient[rhoS, k]] === 0, {k, 0, NQ}];
  emit["(Q3) Sum_j j A_j(t) = -c^2 rho(t), t^0..t^" <> ToString[NQ] <> ": " <> ToString[ok] <> tm[]]];

(* ---------- (M7) ---------- *)
Module[{NW = 12, FsS, P1, P2, w0},
  FsS = Series[Sum[tau^n Product[(1 + j z)^2, {j, 1, n}]/(n! z^n Product[c + j z, {j, 2, n + 1}]), {n, 0, NW + 2}], {tau, 0, NW + 2}];
  P1[zz0_] := FsS /. z -> zz0;
  P2[zz0_] := Module[{F = FsS /. z -> zz0}, zz0 (1 - tau) D[F, tau] + ((c + zz0)/tau - (1 + zz0)) F];
  w0 = MapAt[Together, tau^2 (x^2 (P1[z] P2[-z] + P2[z] P1[-z]) - c x P2[z] P2[-z]), 3];
  emit["(M7) Remark 4.9, tau^2 <PhiF(z), PhiF(-z)> = -x c (c^2 - z^2) (symbolic x, z), through tau^" <> ToString[w0[[5]] - 1] <> ": " <>
    ToString[Together[Normal[w0] + x c (c^2 - z^2)] === 0] <> tm[]]];

(* ---------- (M4) ---------- *)
Module[{NTr = 10, MM = 5, delT, sstT, algT, RtrT, TtrOp, Phi, ok0, dphiT, Fs, prodS, Hn, okReg, okm},
  delT = Series[c Sqrt[1 - 4 x tau/c^2], {tau, 0, NTr + 2}];
  sstT = tidyC[2/(c + delT)];
  algT = tidyC[Sqrt[delT/c] (c + delT)/(2 c)];
  RtrT = tidyC[-x (c + 3 delT)/(delT (c + delT))];
  Phi[0] = tidyC[Series[Exp[-Integrate[Normal[tidyC[RtrT/delT]] /. tau -> uu, {uu, 0, tau}]], {tau, 0, NTr}]];
  ok0 = Normal[tidyC[Phi[0] - 1/algT]] === 0;
  TtrOp[Y_] := tau (1 - tau) D[Y, {tau, 2}] + (2 - 3 tau) D[Y, tau] - Y;
  Do[Phi[m] = Module[{ig = Normal[tidyC[Series[TtrOp[Normal[Phi[m - 1]]], {tau, 0, NTr - 1}]/(delT Phi[0])]]},
      tidyC[Series[-Normal[Phi[0]] Integrate[ig /. tau -> uu, {uu, 0, tau}], {tau, 0, NTr}]]], {m, 1, MM}];
  dphiT = Integrate[Normal[tidyC[sstT/(1 - tau sstT)]] /. tau -> uu, {uu, 0, tau}];
  Fs = Sum[tau^n Product[(1 + j z)^2, {j, 1, n}]/(n! z^n Product[c + j z, {j, 2, n + 1}]), {n, 0, NTr}];
  prodS = Series[Exp[-dphiT/z] Fs, {tau, 0, NTr}];
  Hn = Table[Together[SeriesCoefficient[prodS, n]], {n, 0, NTr}];
  okReg = And @@ Table[(Denominator[Hn[[n + 1]]] /. z -> 0) =!= 0, {n, 0, NTr}];
  okm = Table[Together[Sum[Coefficient[Normal[Series[Hn[[n + 1]], {z, 0, MM}]], z, m] tau^n, {n, 0, NTr}] - Normal[Phi[m]]] === 0, {m, 0, MM}];
  emit["(M4) Proposition 3.19: Phi_0 = Exp[-Int R_tr/delta] = 1/alg: " <> ToString[ok0] <> "; [tau^n](Exp[-dphi/z] F) regular at z = 0, n <= " <> ToString[NTr] <> ": " <>
    ToString[okReg] <> "; transport solutions Phi_m = [z^m](Exp[-dphi/z] F) through tau^" <> ToString[NTr] <> ", m = 0.." <> ToString[MM] <> ": " <> ToString[okm] <> tm[]]];

(* ---------- (M3) ---------- *)
Module[{M = 15, wz, phi0k, lhs1, rhs1, th, subs2, lhs2, hw, hwm},
  wz = wickDirect[M];
  phi0k[k_] := (k - 1)! c^k ((-1)^(k - 1) - x^(1 - k));
  lhs1 = Table[Together[wz[[m + 1]] /. Table[Wk[k] -> phi0k[k], {k, 2, 2 M + 2}]], {m, 0, M}];
  hw = hCoeffs[{c, -x, -1}, 1, M];
  rhs1 = Table[Together[Sum[hw[[j + 1]] (-1/c)^(m - j), {j, 0, m}]], {m, 0, M}];
  th[0] = c Log[1 - U] - Log[U]; Do[th[k] = Together[U D[th[k - 1], U]], {k, 1, 2 M + 2}];
  subs2 = Table[Wk[k] -> Together[th[k] /. U -> -1/x], {k, 2, 2 M + 2}];
  lhs2 = Table[Together[wz[[m + 1]] /. subs2], {m, 0, M}];
  hwm = hCoeffs[{c, -x, -1}, -1, M];
  emit["(M3) Proposition 3.15: wickDirect at W_k -> phi_0^(k)(1/c) equals [z^m] h_w(z) c/(c + z), m = 0.." <> ToString[M] <> ": " <> ToString[Union[Together[lhs1 - rhs1]] === {0}] <>
    "; Proposition 8.6(a): at the regular point U = -1/x (Q = 1) it equals [z^m] h_w(z)^-1: " <> ToString[Union[Together[lhs2 - hwm]] === {0}] <>
    "; z^1 coefficient " <> ToString[InputForm[Factor[lhs1[[2]]]]] <> tm[]]];

(* ---------- (M9) ---------- *)
m9body[NT9_, M9_] := Module[{delT, sstT, sig, U2T, fk, f2t, aaT, L2c, L2, f0T, f2T2, DopS, al, ga, gaussD},
   delT = Series[c Sqrt[1 - 4 x tau/c^2], {tau, 0, NT9}];
   sstT = tidyC[2/(c + delT)];
   sig = tidyC[1/(x sstT)];                         (* sigma~ = tau s~ = (c + delta)/(2 x) *)
   U2T = tidyC[1 - sig];
   fk[k_] := fk[k] = tidyC[(k - 1)! ((-1)^(k - 1) sig^-k - x (tau - sig)^-k + (1 - sig)^-k)];   (* d^k/dsigma^k of phi_tau(sigma/tau) *)
   f2t = fk[2];
   aaT[a_] := aaT[a] = tidyC[U2T^(-a - 1)];
   L2c = lnormS[fk, f2t, aaT, M9];
   L2 = Sum[L2c[[m + 1]] z^m, {m, 0, M9}];
   f0T = tidyC[(sig/tau)/U2T];                      (* f0_tau = s~/(1 - tau s~) *)
   f2T2 = tidyC[1/tau + D[f2t, tau]/(2 f2t)];       (* f2_tau/(2 f2), f2 = phi_tau''(s~) = tau^2 f2t *)
   DopS[Y_] := tidyC[D[Y, tau] + (f0T/z - f2T2) Y];
   al = 1/z + 1; ga = c/z + 2;
   gaussD = tidyC[tau (1 - tau) DopS[DopS[L2]] + (ga - (2 al + 1) tau) DopS[L2] - al^2 L2];
   {gaussD[[4]], gaussD[[5]] - 1, Table[Normal[MapAt[Together[Coefficient[Expand[#], z, j]] &, gaussD, 3]] === 0, {j, -2, M9 - 2}]}];
Module[{r = m9body[8, 7]},
  emit["(M9) (8.5), Gauss_Dop[L2] = 0 (symbolic x), tau-orders " <> ToString[r[[1]]] <> ".." <> ToString[r[[2]]] <> ", z^-2..z^5: " <> ToString[r[[3]]] <> tm[]]];

(* ---------- (M2) ---------- *)
m2body[NT2_, MZ2_] := Module[{delT, sstT, U1T, fk, aaT, Lc, Nc, dphiT, algT, Fs, prodS, Psi, hwz, abz, regOK, ok},
   delT = Series[c Sqrt[1 - 4 x tau/c^2], {tau, 0, NT2 + 1}];
   sstT = tidyC[2/(c + delT)];
   U1T = tidyC[1 - tau sstT];
   fk[k_] := fk[k] = tidyC[(k - 1)! ((-1)^(k - 1) sstT^-k - x (1 - sstT)^-k + tau^k U1T^-k)];   (* (3.6) at sst *)
   aaT[a_] := aaT[a] = tidyC[tau^a U1T^(-a - 1)];                                               (* Taylor data of (1 - tau s)^-1 *)
   Lc = lnormS[fk, fk[2], aaT, MZ2];
   Nc = Table[tidyC[U1T Lc[[m + 1]]], {m, 0, MZ2}];
   dphiT = Integrate[Normal[tidyC[sstT/(1 - tau sstT)]] /. tau -> uu, {uu, 0, tau}];
   algT = tidyC[Sqrt[delT/c] (c + delT)/(2 c)];
   Fs = Sum[tau^n Product[(1 + j z)^2, {j, 1, n}]/(n! z^n Product[c + j z, {j, 2, n + 1}]), {n, 0, NT2}];
   prodS = Series[Fs Exp[-dphiT/z] Normal[algT], {tau, 0, NT2}];
   Psi = Table[Together[SeriesCoefficient[prodS, n]], {n, 0, NT2}];
   regOK = And @@ Table[(Denominator[Psi[[n + 1]]] /. z -> 0) =!= 0, {n, 0, NT2}];
   hwz = hCoeffs[{c, -x, -1}, 1, MZ2];
   abz = Sum[hwz[[j + 1]] z^j, {j, 0, MZ2}] Sum[(-z/c)^j, {j, 0, MZ2}];
   ok = And @@ Flatten[Table[Together[sc[Nc[[m + 1]], n] - Coefficient[Normal[Series[abz Psi[[n + 1]], {z, 0, MZ2}]], z, m]] === 0, {n, 0, NT2}, {m, 0, MZ2}]];
   {regOK, ok, "sample [tau^1 z^1] N = " <> ToString[InputForm[Together[sc[Nc[[2]], 1]]]] <> " (expected (1-x)^2/c^3 = " <> ToString[InputForm[Together[(1 - x)^2/c^3]]] <> ")"}];
Do[Module[{r = withX[x0, m2body[8, 9] &]},
   emit["(M2) (3.17) at x = " <> ToString[InputForm[x0]] <> ", tau^0..tau^8, z^0..z^9: {Psi_n regular at z = 0, U_1 Lnorm in the coordinate s = h_w c/(c + z) F Exp[-dphi/z] alg}: " <> ToString[r] <> tm[]]],
  {x0, {2, 1/3, -5/2, 3/7}}];

(* ---------- (M10), (M11) ---------- *)
m10body[NE_, MZ_] := Module[{tauS, U1e, sst, dphi, F21, Gz, Gk, Gzc, u0, u1, L1s, pr, Ed, ThetaS, Th, okNeg, rhs75, ok75, zetaS, hyp, Cj, ahat, aj, gj, ejk,
    okLeaf, Pk, okAnti, eilT, okInt},
   tauS = Normal[Series[1 - Exp[-t], {t, 0, NE + 4}]];
   U1e = ((x - 1)/x + ((1 + x)/x) Sqrt[1 + 4 x (Exp[-t] - 1)/(1 + x)^2])/2;
   sst = Normal[tidyC[Series[(1 - U1e)/(1 - Exp[-t]), {t, 0, NE + 3}]]];
   dphi = Integrate[Normal[tidyC[Series[1 - x sst, {t, 0, NE + 2}]]] /. t -> uu, {uu, 0, t}];
   F21 = Sum[(tauS/z)^k Product[(1 + j z)^2, {j, 1, k}]/(k! Product[c + j z, {j, 2, k + 1}]), {k, 0, NE + 3}];
   Gz = MapAt[Together, Series[(c/(c + z)) (tauS/t) F21 Exp[-dphi/z], {t, 0, NE + 2}], 3];
   Gk = Table[Together[SeriesCoefficient[Gz, k]], {k, 0, NE + 2}];
   Gzc = Table[PadRight[CoefficientList[Normal[Series[Gk[[k + 1]], {z, 0, MZ}]], z], MZ + 1], {k, 0, NE + 2}];   (* Gzc[[k+1, p+1]] = [t^k z^p] G *)
   u0 = Table[Sum[Gzc[[k + 1, p + 1]] t^k, {k, 0, NE + 2}], {p, 0, MZ}];
   L1s = Normal[tidyC[Series[c/(1 - Exp[-t]) - x sst, {t, 0, NE + 1}]]];
   u1 = Table[Expand[L1s u0[[p + 1]] + If[p >= 1, u0[[p]]/t + D[u0[[p]], t], 0]], {p, 0, MZ}];
   pr[p_, q_] := pr[p, q] = Expand[x^2 (u0[[p + 1]] u1[[q + 1]] + u1[[p + 1]] u0[[q + 1]]) - c x u1[[p + 1]] u1[[q + 1]]];
   Ed[m_, n_] := Ed[m, n] = Expand[-(-1)^(m + n + 1) Sum[(-1)^k pr[m + 1 + k, n - k], {k, 0, n}]];   (* rho^2 EdN[1,1][m, n]; exact for t^-2..t^NE *)
   ThetaS = Normal[tidyC[Series[-x^3 tauS sst^2, {t, 0, NE + 1}]]]; Th[k_] := Coefficient[ThetaS, t, k];
   okNeg = And @@ Flatten[Table[If[m + n <= MZ - 1, {Together[Coefficient[Ed[m, n], t, -2]] === 0, Together[Coefficient[Ed[m, n], t, -1]] === 0}, True], {m, 0, MZ - 1}, {n, 0, MZ - 1}]];
   rhs75[m_, n_, k_] := (1/(k + 2)) Sum[Th[k0] (k1 + 1) (k + 2 - k0 - k1) (-1)^(m + n) Gzc[[k1 + 1, m + 1]] Gzc[[k + 2 - k0 - k1, n + 1]], {k0, 1, k + 1}, {k1, 0, k + 1 - k0}];
   ok75 = And @@ Flatten[Table[If[m + n <= MZ - 1, Table[Together[Coefficient[Ed[m, n], t, k] - rhs75[m, n, k]] === 0, {k, 0, NE}], True], {m, 0, MZ - 1}, {n, 0, MZ - 1}]];
   zetaS = Normal[tidyC[Series[tauS Exp[dphi/c], {t, 0, NE + 3}]]];
   hyp[j_] := Sum[Pochhammer[j x/c, mm]^2/(Pochhammer[j, mm] mm!) tauS^mm, {mm, 0, NE + 2}];
   Cj[j_] := c^2 Pochhammer[1 - j/c, j - 1]^2/(j! (j - 1)!);
   ahat[j_] := ahat[j] = Normal[tidyC[Series[(-1)^j Cj[j] zetaS^j hyp[j], {t, 0, NE + 2}]]];
   aj[j_, k_] := Coefficient[ahat[j], t, k + 1];
   gj[j_, k_] := -(j/c) aj[j, k];
   ejk[i_, l_, k_] := ejk[i, l, k] = Together[(1/(k + 2)) Sum[Th[k0] (k1 + 1) (k + 2 - k0 - k1) gj[i, k1] gj[l, k + 1 - k0 - k1], {k0, 1, k + 1}, {k1, 0, k + 1 - k0}]];
   okLeaf = And @@ Flatten[Table[If[m + n <= MZ - 1, Table[Together[Coefficient[Ed[m, n], t, k] - Sum[ejk[i, l, k] i^m l^n/c^(m + n + 2), {i, 1, k + 1}, {l, 1, k + 1}]] === 0, {k, 0, NE}], True],
       {m, 0, MZ - 1}, {n, 0, MZ - 1}]];
   Pk[mm_] := Product[(mm + 1 - i c)^2, {i, 1, mm}];
   okAnti = And @@ Flatten[Table[Together[ejk[k1 + 1, k - k1 + 1, k] - (-1)^(k + 1) x^3 (k1 + 1) (k - k1 + 1) Pk[k1] Pk[k - k1]/((k + 2) (k1!)^2 ((k - k1)!)^2 c^(2 k))] === 0,
       {k, 0, NE}, {k1, 0, k}]];
   eilT[i_, l_] := Expand[(i l/c^2) Integrate[Normal[tidyC[Series[ThetaS D[ahat[i], t] D[ahat[l], t], {t, 0, NE + 1}]]] /. t -> uu, {uu, 0, t}]/t^2];
   okInt = And @@ Flatten[Table[If[i + l - 2 <= NE, Module[{e = eilT[i, l]}, Table[Together[Coefficient[e, t, k] - ejk[i, l, k]] === 0, {k, 0, NE}]], True], {i, 1, NE + 1}, {l, 1, NE + 1}]];
   {okNeg, ok75, okLeaf, okAnti, okInt, "sample [t^0..t^2] rho^2 E_11(0,0) = " <> ToString[InputForm[Table[Together[Coefficient[Ed[0, 0], t, k]], {k, 0, 2}]]] <> " ([t^0] expected -x^3/(2 c^2) = " <> ToString[InputForm[Together[-x^3/(2 c^2)]]] <> ")"}];
Do[Module[{r = withX[x0, m10body[8, 15] &]},
   emit["(M10/M11) x = " <> ToString[InputForm[x0]] <> ": rho^2 E_11 built from the closed form alone by the code's difference quotient, m + n <= 14: {no t^-2, t^-1 terms, (7.5) first line t^0..t^8, " <>
     "(7.5) leaf decomposition with a_j^(k) from Theorem 11.1, anti-diagonal (7.6), Theorem 11.4 integral = e^(k)_il}: " <> ToString[r] <> tm[]]],
  {x0, {2, 3/7, -5/2}}];

(* ---------- (M2) symbolic ---------- *)
Module[{r = m2body[5, 6]},
  emit["(M2s) (3.17) symbolically in x, tau^0..tau^5, z^0..z^6: " <> ToString[r] <> tm[]]];
emit["total time " <> ToString[Round[AbsoluteTime[] - t0]] <> " s"];
