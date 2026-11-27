(* Collector: assembles the sample files of all jobs, checks the eps-structure and reconstructs F_{g,0} exactly.
   Usage:  wolframscript -file genRmat_g4_collect.m genus
   Input:  samples_g<genus>_job*.m in the script directory (any number of jobs; a t-sample is used only if all its
           eps-ratio points are present; the log lists what is missing).
   Output: genRmat_F<genus>0_cluster_result.m  = {Fg, Fgred, mons, FgQrrc, FgredQrrc}
             Fg in y, XX = y c1'/c1, eps1, eps2, q0..q4;  Fgred = Fg/(eps1 eps2) in LL = (1+27y)^(-1/3);
             FgQrrc, FgredQrrc the same with q_k -> Qrrc[0,k];
           genRmat_F<genus>0_cluster_result.txt (FgredQrrc in InputForm), log in genRmat_g4_collect.log. *)
SetDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
If[Length[args] < 1, Print["usage: wolframscript -file genRmat_g4_collect.m genus"]; Exit[1]];
gg = ToExpression[args[[1]]];
Get["genRmat_g4_lib.m"];
log[m_] := Module[{s = OpenAppend["genRmat_g4_collect.log"]}, WriteString[s, DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}], "  [g=", gg, "] ", m, "\n"]; Close[s]; Print[m]];
wmax = 3 gg - 3; mons = qmonos[wmax];
xsamples = xsamplesFor[gg]; tall = tsamplesFor[gg];
files = FileNames["samples_g" <> ToString[gg] <> "_job*.m"];
If[files === {}, log["no sample files found"]; Exit[1]];
sample = <||>; diagmax = 0; badpts = {};
Do[Do[sample[r[[1]]] = r[[2]]; diagmax = Max[diagmax, "imag" /. r[[3]]]; If[("bad" /. r[[3]]) =!= {}, AppendTo[badpts, r[[1]]]], {r, ReadList[f]}], {f, files}];
log[ToString[Length[files]] <> " files, " <> ToString[Length[sample]] <> " of " <> ToString[Length[pointList[gg]]] <> " sample points; max |Im| " <> ToString[diagmax] <> "; points with unexpected monomials: " <> ToString[badpts]];
If[! And @@ (Dimensions[#] == {wmax + 1, Length[mons]} & /@ Values[sample]), log["ERROR: sample matrices of the wrong shape (different genus or KMAX?)"]; Exit[1]];
tsamples = Select[tall, And @@ Table[KeyExistsQ[sample, {#, x0}], {x0, xsamples}] &];
log[ToString[Length[tsamples]] <> " complete t-samples: " <> ToString[tsamples] <> "; incomplete: " <> ToString[Complement[tall, tsamples]]];
If[Length[tsamples] < 8, log["too few complete t-samples for a reconstruction"]; Exit[1]];

(* eps-structure: x * coeff(x) = sum_{j=0}^{2g+w} a_j x^j (eps2 = 1), palindromic; verified on the 2 spare x-samples *)
aCoef = <||>; warn = 0;
Do[Module[{w = mons[[m, 2]], DD, A, rhs, sol, chk},
   DD = 2 gg + w;
   A = Table[x0^j, {x0, xsamples[[1 ;; DD + 1]]}, {j, 0, DD}];
   rhs = Table[x0 sample[{t0, x0}][[k + 1, m]], {x0, xsamples[[1 ;; DD + 1]]}];
   sol = LinearSolve[A, rhs];
   chk = Table[Sum[sol[[j + 1]] x0^j, {j, 0, DD}] - x0 sample[{t0, x0}][[k + 1, m]], {x0, xsamples[[DD + 2 ;;]]}];
   If[DeleteCases[chk, 0] =!= {}, warn++; log["WARNING: eps-polynomial of degree " <> ToString[DD] <> " fails at t=" <> ToString[t0] <> " k=" <> ToString[k] <> " q-monomial " <> ToString[qmono[mons[[m, 1]]], InputForm] <> ": " <> ToString[N[Max[Abs[chk]], 3]]]];
   If[sol =!= Reverse[sol], warn++; log["WARNING: eps-polynomial not palindromic at t=" <> ToString[t0] <> " k=" <> ToString[k] <> " q-monomial " <> ToString[qmono[mons[[m, 1]]], InputForm]]];
   Do[aCoef[{t0, k, m, j}] = sol[[j + 1]], {j, 0, DD}]], {t0, tsamples}, {k, 0, wmax}, {m, Length[mons]}];
log["eps structure done, " <> ToString[warn] <> " warnings"];

(* y-dependence: rational reconstruction verified on >= 5 spare t-samples *)
Fg = 0; failed = {};
Do[Module[{w = mons[[m, 2]], pts, r},
   pts = Table[{t0^3, aCoef[{t0, k, m, j}]}, {t0, tsamples}];
   r = reconstruct[pts] /. u -> y;
   If[r === $Failed, AppendTo[failed, {k, m, j}], Fg += r XX^k \[Epsilon][1]^j \[Epsilon][2]^(2 gg + w - j) qmono[mons[[m, 1]]]]], {k, 0, wmax}, {m, Length[mons]}, {j, 0, 2 gg + mons[[m, 2]]}];
If[failed =!= {}, log["RECONSTRUCTION FAILED for (k, m, j) = " <> ToString[failed] <> " (more t-samples needed)"]];
Fg = Fg/(\[Epsilon][1] \[Epsilon][2]);
Fgred = Collect[Expand[Together[Fg/(\[Epsilon][1] \[Epsilon][2]) /. y -> (1 - LL^3)/(27 LL^3)]], Join[{XX}, qsyms], Factor];
FgQrrc = Fg /. qrule; FgredQrrc = Fgred /. qrule;
Put[{Fg, Fgred, mons, FgQrrc, FgredQrrc}, "genRmat_F" <> ToString[gg] <> "0_cluster_result.m"];
Export["genRmat_F" <> ToString[gg] <> "0_cluster_result.txt", ToString[InputForm[FgredQrrc]], "Text"];
log["F" <> ToString[gg] <> "0 reconstructed from " <> ToString[Length[tsamples]] <> " t-samples; denominators: " <> ToString[Union[Cases[Fg, Power[b_, e_?Negative] :> b, Infinity]]] <> "; LeafCount " <> ToString[LeafCount[Fg]] <> "; written to genRmat_F" <> ToString[gg] <> "0_cluster_result.m"];
