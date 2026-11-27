ord = 3;
F2 = -y1^2/2 - y2^2 + y1 y2/3 + y1^3/5 + y1^2 y2/2 - y1 y2^2/7 + y2^3/3 + y1^4/11 + y1^2 y2^2/13 - y1^3 y2/4;
ph1 = y1 + y1 y2 + y2^2 - y1^3/2; ph2 = 2 y2 + y1^2 - y1 y2/3 + y2^3 + y1^2 y2;
trunc[ex_] := Normal[Series[ex /. {y1 -> r y1, y2 -> r y2}, {r, 0, 2 ord + 4}]] /. r -> 1;
Fphi = trunc[F2 /. {y1 -> ph1, y2 -> ph2}];
Print[Fphi /. {y1 -> 0, y2 -> 0}];
Print[Short[Fphi, 5]];
