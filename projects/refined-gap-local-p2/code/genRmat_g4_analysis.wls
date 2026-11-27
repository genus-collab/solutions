(* Conifold-frame expansion and gap conditions for a twisted free energy F_{g,0} produced by the cluster pipeline.
   Usage:  wolframscript -file genRmat_g4_analysis.wls genus
   Input:  cluster_genus4/genRmat_F<genus>0_cluster_result.m  (from genRmat_g4_collect.m; for genus <= 3 the notebook
           result genRmat_F<genus>0_result.mx is used if the cluster file is absent),
           genRmat_conifold.mx (t_c, invcf, C1cf, q1sol, q2sol from localp2CY5_genRmat_genus3.nb).
   The gap relations of the lower genera (Qrrc[0,1] from genus 2, Qrrc[0,2] from genus 3) are substituted; the
   remaining coefficients entering at this genus (Qrrc[0,k] with 2k+1 <= 3g-3) are solved for from the vanishing
   of t_c^-(2g-3) .. t_c^-1, the consistency of these conditions is checked, and the substitution observation
   (eps1+eps2)^2 -> (eps1+eps2)^2 + 12 Qrrc[0,0] eps1 eps2 (eps1+eps2) is tested on every coefficient.
   Output: genRmat_F<genus>0_gap_analysis.txt and printed log. *)
SetDirectory[DirectoryName[$InputFileName]];
gg = ToExpression[$ScriptCommandLine[[2]]];
out = {}; log[m_] := (AppendTo[out, m]; Print[m]);
cfile = "cluster_genus4/genRmat_F" <> ToString[gg] <> "0_cluster_result.m";
If[FileExistsQ[cfile], FgQ = Get[cfile][[4]]; log["loaded " <> cfile],
  Get["genRmat_F" <> ToString[gg] <> "0_result.mx"]; FgQ = Fgtw /. Table[Symbol["q" <> ToString[k]] -> Qrrc[0, k], {k, 0, 2}]; log["loaded notebook result genRmat_F" <> ToString[gg] <> "0_result.mx"]];
Get["genRmat_conifold.mx"];   (* tc (in Del), invcf (in tcv), C1cf, cf[2], cf[3], q1sol, q2sol *)
NC = 2 gg + 4;
If[Exponent[tc, Del] < NC, log["WARNING: t_c known to Del^" <> ToString[Exponent[tc, Del]] <> " only; NC = " <> ToString[NC] <> " needed -- recompute t_c with NC terms"];
  Module[{L2b, wc, ceq, csol, wcDel}, L2b[f_] := Module[{th = y D[#, y] &}, (th[th[f]] + 3 y (9 th[th[f]] + 9 th[f] + 2 f))];
   wc = Sum[a[n] Del^n, {n, 0, NC}] /. Del -> 1 + 27 y;
   ceq = CoefficientList[Normal[Series[L2b[wc] /. y -> (Del - 1)/27, {Del, 0, NC - 1}]], Del];
   csol = Solve[Thread[ceq == 0] /. a[0] -> 1, Table[a[n], {n, 1, NC}]];
   wcDel = Sum[a[n] Del^n, {n, 0, NC}] /. First[csol] /. a[0] -> 1;
   tc = Module[{ser = Normal[Series[wcDel/(Del - 1), {Del, 0, NC}]]}, Sum[Coefficient[ser, Del, n] Del^(n + 1)/(n + 1), {n, 0, NC}]];
   tc = Normal[Series[tc/SeriesCoefficient[tc, {Del, 0, 1}], {Del, 0, NC}]];
   invcf = InverseSeries[Series[tc, {Del, 0, NC}], tcv]; C1cf = Normal[y 27 D[tc, Del]] /. Del -> 1 + 27 y]];
cfexp[f_] := Series[f /. XX -> y D[C1cf, y]/C1cf /. y -> (Del - 1)/27 /. Del -> invcf, {tcv, 0, 1}];
cfg = cfexp[FgQ];
known = If[gg <= 3, First[q1sol], Join[First[q1sol], First[q2sol[-1]]]];   (* gap relations of the lower genera: genus 2 fixes Qrrc[0,1], genus 3 fixes Qrrc[0,2] *)
knownQ = known[[All, 1]];
newQ = Select[Table[Qrrc[0, k], {k, 0, Floor[(3 gg - 4)/2]}], ! MemberQ[knownQ, #] && # =!= Qrrc[0, 0] &];
log["genus " <> ToString[gg] <> ": twist coefficients entering: " <> ToString[Table[Qrrc[0, k], {k, 0, Floor[(3 gg - 4)/2]}], InputForm] <> "; fixed by lower genera: " <> ToString[knownQ, InputForm] <> "; to be determined: " <> ToString[newQ, InputForm]];
Do[log["coefficient of t_c^" <> ToString[n] <> ": " <> ToString[Factor[SeriesCoefficient[cfg, n]], InputForm]], {n, -(2 gg - 2), 1}];
Do[cs[n] = Factor[SeriesCoefficient[cfg, n] /. known], {n, -(2 gg - 2), 1}];
gapeqs = Table[cs[n] == 0, {n, -(2 gg - 3), -1}];
If[newQ =!= {},
  sol = Solve[gapeqs, newQ];
  log["gap conditions t_c^-" <> ToString[2 gg - 3] <> " .. t_c^-1 solved simultaneously for " <> ToString[newQ, InputForm] <> ": " <> ToString[FullSimplify[sol], InputForm]];
  If[sol === {}, log["INCONSISTENT gap conditions"];
    Do[log["  single condition t_c^" <> ToString[n] <> " = 0: " <> ToString[FullSimplify[Solve[cs[n] == 0, newQ]], InputForm]], {n, -(2 gg - 3), -1}],
    full = Join[known, First[sol]];
    log["all gap coefficients after substitution: " <> ToString[Simplify[Table[cs[n] /. First[sol], {n, -(2 gg - 3), -1}]], InputForm]];
    log["leading coefficient t_c^-" <> ToString[2 gg - 2] <> " with the gap relations: " <> ToString[Factor[cs[-(2 gg - 2)] /. First[sol]], InputForm]];
    log["regular part t_c^0, t_c^1 with the gap relations: " <> ToString[Factor /@ {cs[0] /. First[sol], cs[1] /. First[sol]}, InputForm]];
    toSP[f_] := Module[{g = Together[f]}, SymmetricReduction[Numerator[g], {\[Epsilon][1], \[Epsilon][2]}, {sig, pp}][[1]]/SymmetricReduction[Denominator[g], {\[Epsilon][1], \[Epsilon][2]}, {sig, pp}][[1]]];
    subcheck[c_] := Module[{tw = Together[c /. full], un}, un = tw /. Qrrc[0, 0] -> 0; Simplify[(toSP[un] /. sig -> Sqrt[sig^2 + 12 Qrrc[0, 0] pp sig]) - toSP[tw]]];
    log["substitution observation (eps1+eps2)^2 -> (eps1+eps2)^2 + 12 Qrrc[0,0] eps1 eps2 (eps1+eps2), t_c^-" <> ToString[2 gg - 2] <> " .. t_c^1 (0 = holds): " <> ToString[Table[subcheck[SeriesCoefficient[cfg, n]], {n, -(2 gg - 2), 1}], InputForm]];
    qrr[k_] := -BernoulliB[2 k + 2]/((2 k + 2) (2 k + 1)) (-12 Qrrc[0, 0])^(2 k + 1);
    log["pure-Qrrc[0,0] parts of the new relations minus the QRR coefficients of the weight -1/(12 Qrrc[0,0]): " <> ToString[Simplify[Table[Coefficient[Qrrc[0, k] /. full, Qrrc[0, 0], 2 k + 1] Qrrc[0, 0]^(2 k + 1) - qrr[k], {k, Select[Range[1, 4], MemberQ[newQ, Qrrc[0, #]] &]}]], InputForm]]]];
Export["genRmat_F" <> ToString[gg] <> "0_gap_analysis.txt", StringRiffle[out, "\n\n"], "Text"];
