(* hvc_multiprobe_g4b.wls -- genus-4 irreducible points of Conjecture M with r = 3, 4, 5 leaves (all parts >= 2, D <= 5),
   exact in x, plus control points at D = 6 (predicted non-zero). Xi evaluated at the lattice point directly:
   sum over sorted exponent multisets a, symmetrised over the distinct permutations. Output: hvc_multiprobe_g4b_out.txt *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_mp_header.m"];
out = {}; emit[s_] := (AppendTo[out, s]; Print[s]; Export["hvc_multiprobe_g4b_out.txt", StringRiffle[out, "\n"], "Text"]);
XiAt[h_, jl_List] := Module[{r = Length[jl], amax, tot = 0},
   amax = 3 h - 3 + r;
   Do[tot += (Times @@ (1/c^(a + 1))) (Total[Table[Times @@ (jl^p), {p, Permutations[a]}]]) Lval[h, a], {a, msets[r, 0, amax, amax]}];
   Together[tot]];
pts = {{2, 2, 2}, {2, 2, 3}, {2, 2, 4}, {2, 3, 3}, {2, 2, 5}, {2, 3, 4}, {3, 3, 3},
   {2, 2, 2, 2}, {2, 2, 2, 3}, {2, 2, 2, 4}, {2, 2, 3, 3}, {2, 2, 2, 2, 2}, {2, 2, 2, 2, 3}};
t0 = AbsoluteTime[];
Do[Module[{t1 = AbsoluteTime[], v, dd = Total[p - 1]},
   v = XiAt[4, p];
   emit["genus 4, R = " <> ToString[p] <> ", D = " <> ToString[dd] <> " (" <> If[dd <= 5, "in range, predicted 0", "outside, predicted non-zero"] <>
     "): " <> If[v === 0, "0", "non-zero, " <> ToString[InputForm[Factor[v]]]] <> "   [" <> ToString[Round[AbsoluteTime[] - t1, 0.1]] <> " s]"]], {p, pts}];
emit["total " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s"];
