package RefHodge;
# Referee's independent Hodge-integral engine (exact, Math::BigRat).  No code shared with ../Hodge.pm, ../HodgePt.pm, ../WK.pm.
# No topological recursion of any kind.
#
#   I(g, mu, e) = Int_{Mbar_{g,n}} ch_1^{mu_1} ch_3^{mu_3} ch_5^{mu_5} ch_7^{mu_7} psi_1^{e_1} ... psi_n^{e_n}
#
# Differences from the author's engine (on purpose):
#  * NO kappa classes anywhere: the kappa_l term of Mumford's formula is written as an extra marked point,
#      Int_{Mbar_{g,n}} kappa_l X = Int_{Mbar_{g,n+1}} psi_{n+1}^{l+1} X   (X = ch-monomial x psi's; psi_{n+1} D_{i,n+1} = 0),
#  * string and dilaton equations (valid since ch(E) is pulled back) remove points with exponent 0 and 1,
#  * dilaton leaves are kept as explicit marked points (the definition of Xi in the brief), so the kappa form K, Prop. W
#    and the point recursion of HodgePt.pm are NOT used,
#  * lambda_k from Newton's identities (power sums p_l = l! ch_l, even ones zero), Bernoulli numbers computed here.
# Mumford (Faber-Pandharipande, Invent. Math. 139 (2000) s.2):  for l odd,
#   ch_l(E) = B_{l+1}/(l+1)! [ kappa_l - Sum_i psi_i^l + 1/2 Sum_{boundary} iota_*( Sum_{p=0}^{l-1} (-1)^p psi_*^p psi_o^{l-1-p} ) ],
# the boundary sum running over the irreducible map Mbar_{g-1,n+2} and ALL ordered splittings (g1, S | g2, S^c).
use strict; use warnings; use Math::BigRat; use Math::BigInt; use RefWK;
sub Q { Math::BigRat->new(shift) }
my $Z = Q(0);
our @CHL = (1, 3, 5, 7, 9)[0 .. ($ENV{REF_NCH} // 4) - 1];   # odd Chern characters available (REF_NCH=5 adds ch_9)
our $NCH = scalar @CHL;
# ---------- Bernoulli numbers (B_1 = -1/2 convention irrelevant: only even ones used) ----------
my @B = (Q(1));
sub bern { my $n = shift; while ($#B < $n) { my $m = scalar @B; my $s = Q(0); for my $k (0..$m-1) { $s += Q(binomZ($m + 1, $k)) * $B[$k]; } push @B, -$s / ($m + 1); } return $B[$n]; }
sub binomZ { my ($n, $k) = @_; return Math::BigInt->new(0) if $k < 0 || $k > $n; my $r = Math::BigInt->new(1); for my $i (1..$k) { $r->bmul($n - $k + $i); $r->bdiv($i); } return $r; }
sub fact { my $n = shift; my $r = Math::BigInt->new(1); $r->bmul($_) for 2..$n; return Q($r); }
my @MC = map { bern($_ + 1) / fact($_ + 1) } @CHL;   # B_{l+1}/(l+1)!
sub binom { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; my $r = 1; for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i } return $r; }
sub msplits { my @x = @_; my (@v, @c); for my $y (@x) { if (@v && $v[-1] == $y) { $c[-1]++ } else { push @v, $y; push @c, 1 } }
  my @out = ([[], [], 1]);
  for my $i (0..$#v) { my @new; for my $o (@out) { for my $s (0..$c[$i]) {
        push @new, [[@{$o->[0]}, ($v[$i]) x $s], [@{$o->[1]}, ($v[$i]) x ($c[$i] - $s)], $o->[2] * binom($c[$i], $s)]; } } @out = @new; }
  return @out; }
our %IM;
our $G1SHORT = 1;   # genus-1 shortcut: every ch-monomial of degree > 1 vanishes on Mbar_{1,n} (lambda_1^2 = 0, ch_l = lambda_1^l/l!)
sub I {
  my ($g, $mu, $e) = @_;
  my $n = scalar @$e;
  return $Z if $g < 0 || 2*$g - 2 + $n <= 0;
  my $se = 0; for (@$e) { return $Z if $_ < 0; $se += $_ }
  my $dm = 0; for my $i (0..$NCH-1) { $dm += $CHL[$i] * ($mu->[$i] // 0) }
  return $Z if $dm + $se != 3*$g - 3 + $n;
  my @ee = sort { $a <=> $b } @$e;
  return RefWK::wk($g, @ee) if $dm == 0;
  return $Z if $g == 0;
  return $Z if $G1SHORT && $g == 1 && $dm > 1;
  return $Z if $g >= 2 && $dm > 3*$g - 3;   # ch(E) is pulled back from Mbar_g
  my $key = "$g|" . join(',', map { $mu->[$_] // 0 } 0..$NCH-1) . "|" . join(',', @ee);
  return $IM{$key} if exists $IM{$key};
  my $res;
  if ($n > 0 && $ee[0] == 0 && 2*$g - 2 + $n - 1 > 0) {          # string
    my @x = @ee[1..$#ee]; $res = Q(0); my %done;
    for my $j (0..$#x) { next if $x[$j] == 0 || $done{$x[$j]}++; my $c = grep { $_ == $x[$j] } @x; my @y = @x; $y[$j]--; $res += I($g, $mu, \@y) * $c; } }
  elsif ($n > 0 && $ee[0] == 1 && 2*$g - 2 + $n - 1 > 0) {       # dilaton
    my @x = @ee[1..$#ee]; $res = I($g, $mu, \@x) * (2*$g - 2 + $n - 1); }
  else {                                               # Mumford on the highest ch_l
    my $li; for (my $i = $NCH - 1; $i >= 0; $i--) { if (($mu->[$i] // 0) > 0) { $li = $i; last } }
    my $l = $CHL[$li]; my @mp = map { $mu->[$_] // 0 } 0..$NCH-1; $mp[$li]--;
    my $t = I($g, \@mp, [@ee, $l + 1]);                # kappa_l as an extra point
    my %done;
    for my $j (0..$#ee) { next if $done{$ee[$j]}++; my $c = grep { $_ == $ee[$j] } @ee; my @y = @ee; $y[$j] += $l; $t -= I($g, \@mp, \@y) * $c; }
    my $half = Q('1/2');
    if ($g >= 1) { for my $p (0..$l-1) { my $v = I($g - 1, \@mp, [@ee, $p, $l - 1 - $p]); next if $v == 0; $t += ($p % 2 ? -$half : $half) * $v; } }
    # reducible boundary, ordered splittings
    my @musp = ([[], [], 1]);
    for my $i (0..$NCH-1) { my @new; for my $o (@musp) { for my $k (0..$mp[$i]) { push @new, [[@{$o->[0]}, $k], [@{$o->[1]}, $mp[$i] - $k], $o->[2] * binom($mp[$i], $k)]; } } @musp = @new; }
    my $dmp = $dm - $l;
    for my $sp (msplits(@ee)) { my ($S, $Sc, $sm) = @$sp; my $n1 = scalar(@$S) + 1; my $n2 = scalar(@$Sc) + 1; my $sS = 0; $sS += $_ for @$S;
      for my $g1 (0..$g) { my $g2 = $g - $g1; next if 2*$g1 - 2 + $n1 <= 0 || 2*$g2 - 2 + $n2 <= 0;
        for my $ms (@musp) { my ($m1, $m2, $mm) = @$ms; my $d1 = 0; $d1 += $CHL[$_] * $m1->[$_] for 0..$NCH-1;
          next if ($g1 == 0 && $d1 > 0) || ($g2 == 0 && $dmp - $d1 > 0);
          my $p = 3*$g1 - 3 + $n1 - $d1 - $sS; next if $p < 0 || $p > $l - 1;
          my $v1 = I($g1, $m1, [@$S, $p]); next if $v1 == 0;
          my $v2 = I($g2, $m2, [@$Sc, $l - 1 - $p]); next if $v2 == 0;
          $t += ($p % 2 ? -$half : $half) * $v1 * $v2 * ($sm * $mm); } } }
    $res = $t * $MC[$li]; }
  $IM{$key} = $res; return $res; }

# ---------- ch-polynomials: hash { "m1,m3,m5,m7" => BigRat } ----------
sub cdeg { my $k = shift; my @m = split /,/, $k; my $d = 0; $d += $CHL[$_] * $m[$_] for 0..$NCH-1; return $d; }
sub cmul { my ($p, $q, $maxd) = @_; my %r;
  for my $a (keys %$p) { my @ma = split /,/, $a; for my $b (keys %$q) { my @mb = split /,/, $b; my @mc = map { $ma[$_] + $mb[$_] } 0..$NCH-1; my $k = join(',', @mc);
      next if cdeg($k) > $maxd; my $c = $p->{$a} * $q->{$b}; if (exists $r{$k}) { $r{$k} += $c } else { $r{$k} = $c } } }
  for (keys %r) { delete $r{$_} if $r{$_} == 0 } return \%r; }
sub cadd { my ($p, $q, $cq) = @_; $cq //= Q(1); my %r = map { $_ => $p->{$_}->copy } keys %$p; for (keys %$q) { my $c = $q->{$_} * $cq; if (exists $r{$_}) { $r{$_} += $c } else { $r{$_} = $c } } for (keys %r) { delete $r{$_} if $r{$_} == 0 } return \%r; }
sub cone { return { join(',', (0) x $NCH) => Q(1) }; }
# lambda_0..lambda_h by Newton's identities: k lambda_k = Sum_{i=1}^k (-1)^{i-1} lambda_{k-i} p_i, p_i = i! ch_i (i odd), 0 (i even)
sub lambdas { my ($h, $maxd) = @_; my @lam = (cone());
  for my $k (1..$h) { my $acc = {};
    for my $i (1..$k) { next unless $i % 2; my ($ix) = grep { $CHL[$_] == $i } 0..$NCH-1; die "need ch_$i (NCH too small)" unless defined $ix;
      my @e = (0) x $NCH; $e[$ix] = 1; my $pi = { join(',', @e) => fact($i) };
      $acc = cadd($acc, cmul($lam[$k - $i], $pi, $maxd), Q(1)); }   # (-1)^{i-1} = +1 for odd i
    my %l = map { $_ => $acc->{$_} / $k } keys %$acc; push @lam, \%l; }
  return @lam; }
# Lambda^vee(w)/w^h = Sum_k (-1)^k lambda_k w^{-k};  Omega = prod over the weights (numeric BigRat list)
sub omega { my ($h, $w, $maxd) = @_; my @lam = lambdas($h, $maxd); my $om = cone();
  for my $wa (@$w) { my $f = {}; for my $k (0..$h) { $f = cadd($f, $lam[$k], (($k % 2) ? Q(-1) : Q(1)) / ($wa ** $k)); } $om = cmul($om, $f, $maxd); }
  return $om; }
# ---------- Xi with dilaton leaves as explicit points ----------
# xi(h, C, js, smax): returns arrayref q[m] = (1/m!) Int_{Mbar_{h,r+m}} C prod_i lam_{j_i}(psi_i) prod_{k=1}^m Tt(psi_{r+k}),
#   lam_j = 1/(1 - j psi), Tt = psi^2/(1 - psi) = Sum_{d>=2} psi^d (c = 1).  Xi_h(R) = Sum_m (-1)^m q[m]; s-polynomial Q(s) = Sum_m s^m q[m].
sub parts_ge1 { my ($n, $max) = @_; $max //= $n; return ([]) if $n == 0; my @o; for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for parts_ge1($n - $k, $k); } return @o; }
my %PC;
sub allparts { my $B = shift; return $PC{$B} //= [ map { parts_ge1($_) } 0..$B ]; }   # all partitions of 0..B (excesses d-1)
sub tuples { my ($r, $t) = @_; return ([]) if $r == 0; return ([$t]) if $r == 1; my @o; for my $a (0..$t) { push @o, [$a, @$_] for tuples($r - 1, $t - $a); } return @o; }
sub xi { my ($h, $C, $js, $smax) = @_; my $r = scalar @$js; my @q = map { Q(0) } 0..(3*$h + $r + 2);
  my @jq = map { Q($_) } @$js;
  my %wt;   # weight of leaf-exponent tuples with sum t: list of [ sorted exps, weight ]
  for my $k (keys %$C) { my $cf = $C->{$k}; next if $cf == 0; my @mu = split /,/, $k; my $dm = cdeg($k);
    my $B = 3*$h - 3 + $r - $dm; next if $B < 0;
    for my $P (@{allparts($B)}) { my $m = scalar @$P; next if defined $smax && $m > $smax; my $ex = 0; $ex += $_ for @$P; my $t = $B - $ex;
      my @D = map { $_ + 1 } @$P; my %mult; $mult{$_}++ for @D; my $den = Q(1); $den *= fact($_) for values %mult;
      my $list = $wt{"$t"} //= leafweights(\@jq, $t);
      my $acc = Q(0);
      for my $it (@$list) { my ($A, $w) = @$it; my $v = I($h, \@mu, [@$A, @D]); next if $v == 0; $acc += $w * $v; }
      next if $acc == 0;
      $q[$m] += $cf * $acc / $den; } }
  return \@q; }
# leaf exponent multisets A (size r, sum t) with weight m_A(j) = Sum_{ordered a, sorted(a) = A} prod_i j_i^{a_i}
#   = coefficient of prod_v y_v^{M_v} in prod_i (Sum_v y_v j_i^v)   (M_v = multiplicity of the value v in A)
sub parts_le { my ($t, $r, $max) = @_; $max //= $t; return ([]) if $t == 0; return () if $r == 0; my @o; for my $k (reverse 1..($max < $t ? $max : $t)) { push @o, [$k, @$_] for parts_le($t - $k, $r - 1, $k); } return @o; }
sub mA { my ($jq, $A) = @_; my (@vals, @M); for my $y (@$A) { if (@vals && $vals[-1] == $y) { $M[-1]++ } else { push @vals, $y; push @M, 1 } }
  my $nv = scalar @vals; my %P = (join(',', (0) x $nv) => Q(1));
  for my $x (@$jq) { my %N; my @pw = map { $vals[$_] ? $x ** $vals[$_] : Q(1) } 0..$nv-1;
    for my $key (keys %P) { my @e = split /,/, $key;
      for my $k (0..$nv-1) { next if $e[$k] >= $M[$k]; my @f = @e; $f[$k]++; my $nk = join(',', @f); my $c = $P{$key} * $pw[$k];
        if (exists $N{$nk}) { $N{$nk} += $c } else { $N{$nk} = $c } } }
    %P = %N; }
  return $P{join(',', @M)} // Q(0); }
sub leafweights { my ($jq, $t) = @_; my $r = scalar @$jq; my @out;
  for my $P (parts_le($t, $r)) { my @A = ((0) x ($r - scalar(@$P)), reverse @$P); my $w = mA($jq, \@A); push @out, [\@A, $w] if $w != 0; }
  return \@out; }
sub xi_value { my ($q, $s) = @_; $s //= Q(-1); my $v = Q(0); my $p = Q(1); for my $m (0..$#$q) { $v += $q->[$m] * $p; $p *= $s; } return $v; }
1;
