#!/usr/bin/perl
# ref_hirota.pl (referee 2): Route A a3_hirota_check.pl sections 1-2 verbatim + a replaced section 3 (Xi values from the R2 formulas or the engine).
# a3_hirota_check.pl -- independent validation of Theorem D (KP at the dilaton point) at higher weight, WITH ALL TERMS:
#  1. generate the KP Hirota equations from the DJKM generating identity
#        sum_{j>=0} p_j(-2y) p_{j+1}(Dt) exp(sum_i y_i D_i) tau.tau = 0,   Dt = (D_1, D_2/2, D_3/3, ...),
#     (coefficient of each y-monomial), up to weight W;
#  2. convert P(D) tau.tau / tau^2 = P(d_u) exp( sum_{|b| even >= 2} 2 F_b u^b / b! )|_{u=0} into polynomials in t-derivatives of F = log tau;
#  3. evaluate every equation at the dilaton point, coefficient of x^{2G+r} p_R, with ALL genus distributions (genus 0..G factors),
#     using the dilaton-point values Xi_g (Laguerre formula, Schur side), Xi_0 = Gamma(J)/Gamma(J-r+3), string equation for unit parts:
#        factor d^b F (t-multi-index b with a unit indices and non-unit indices K) on the block R_j at genus g contributes
#        sigma^{g+|b|-1} (prod_{i in b} i) (prod_{k in K} Pi_k) Xi_g(R_j + K + 1^a),   with sum_j (2g_j - 2 + |b_j|) = 2G.
#  The sum must vanish identically (outside the range all terms are generically non-zero).
# Usage: perl a3_hirota_check.pl G tau W "R1;R2;..."
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($G, $tau, $W, $Rlist) = @ARGV; $tau = R($tau); my $sig = $tau*(1+$tau);
my @Rs = map { [ grep { $_ ne '' } split /,/, $_ ] } split /;/, ($Rlist // '');
@Rs = ([]) unless @Rs;
# ---------- 1. Hirota polynomials ----------
# polynomials: hash "y1..yW|D1..DW" exponent strings -> BigRat ; truncate total y-weight <= W-1, D-weight <= W
my $nv = $W;
sub pkey { my ($ye, $de) = @_; return join(',', @$ye) . '|' . join(',', @$de); }
sub pmul2 { my ($A, $B) = @_; my %C;
  for my $k1 (keys %$A) { my ($y1, $d1) = split /\|/, $k1; my @y1 = split /,/, $y1; my @d1 = split /,/, $d1;
    for my $k2 (keys %$B) { my ($y2, $d2) = split /\|/, $k2; my @y2 = split /,/, $y2; my @d2 = split /,/, $d2;
      my @y = map { $y1[$_] + $y2[$_] } 0..$nv-1; my @d = map { $d1[$_] + $d2[$_] } 0..$nv-1;
      my $wy = 0; $wy += ($_+1)*$y[$_] for 0..$nv-1; my $wd = 0; $wd += ($_+1)*$d[$_] for 0..$nv-1; next if $wy > $W-1 || $wd > $W;
      my $k = pkey(\@y, \@d); my $v = $A->{$k1} * $B->{$k2}; $C{$k} = exists $C{$k} ? $C{$k} + $v : $v; } }
  for my $k (keys %C) { delete $C{$k} if $C{$k}->is_zero; } return \%C; }
sub padd2 { my ($A, $B) = @_; my %C = %$A; for my $k (keys %$B) { $C{$k} = exists $C{$k} ? $C{$k} + $B->{$k} : $B->{$k}; delete $C{$k} if $C{$k}->is_zero; } return \%C; }
my @zero = (0) x $nv;
sub mono { my ($which, $i, $c) = @_; my @y = @zero; my @d = @zero; if ($which eq 'y') { $y[$i-1] = 1; } else { $d[$i-1] = 1; } return { pkey(\@y, \@d) => $c }; }
# elementary Schur polynomials p_j(t) for t given as a linear polynomial t_i = coef * var_i
sub schur_p { my ($which, $scale) = @_;   # returns array of polynomials p_0..p_W ; t_i = scale(i) * var_i
  # exp(sum_i t_i k^i): coefficient of k^j
  my @P = ({ pkey(\@zero, \@zero) => R(1) });
  for my $j (1..$W) { my $acc = {}; for my $i (1..$j) { my $ti = mono($which, $i, $scale->($i)); my $t = pmul2($ti, $P[$j-$i]); for my $k (keys %$t) { $t->{$k} *= $i; } $acc = padd2($acc, $t); }
    for my $k (keys %$acc) { $acc->{$k} /= $j; } $P[$j] = $acc; }   # j p_j = sum_i i t_i p_{j-i}
  return @P; }
my @Py = schur_p('y', sub { R(-2) });
my @PD = schur_p('D', sub { R(1)/$_[0] });
# exp(sum y_i D_i)
my $yd = {}; for my $i (1..$W) { my @y = @zero; my @d = @zero; $y[$i-1] = 1; $d[$i-1] = 1; $yd->{pkey(\@y,\@d)} = R(1); }
my $E = { pkey(\@zero,\@zero) => R(1) }; my $pw = { pkey(\@zero,\@zero) => R(1) };
for my $n (1..$W) { $pw = pmul2($pw, $yd); my %t = map { $_ => $pw->{$_} / fact($n) } keys %$pw; $E = padd2($E, \%t); }
my $Hgen = {}; for my $j (0..$W-1) { $Hgen = padd2($Hgen, pmul2(pmul2($Py[$j], $PD[$j+1]), $E)); }
my %H;   # y-exponent string -> { D-exponent string -> coef }
for my $k (keys %$Hgen) { my ($ys, $ds) = split /\|/, $k; $H{$ys}{$ds} = $Hgen->{$k}; }
# ---------- 2. Hirota polynomial -> polynomial in F-derivatives ----------
# F-derivative = sorted multi-index list (t-indices). [u^m] exp(S): sum over multiset partitions of m into blocks of even size >= 2,
# value prod (2 F_b / b!) / (symmetry of equal blocks) ; then times m!.
sub mpartitions { my @m = @_; return ([]) unless @m; my $first = shift @m; my @out;
  # choose the block containing the first element: subsets of the rest (as index sets), dedupe via sorted keys
  my %seen; my $n = scalar @m;
  for my $mask (0..(1 << $n) - 1) { my @blk = ($first); my @rest; for my $i (0..$n-1) { if (($mask >> $i) & 1) { push @blk, $m[$i]; } else { push @rest, $m[$i]; } }
    next if (scalar(@blk) % 2) == 1; my $bk = join(',', sort { $a <=> $b } @blk); my $rk = join(',', sort { $a <=> $b } @rest);
    next if $seen{"$bk|$rk"}++;
    my $cnt = 0; for my $mask2 (0..(1 << $n) - 1) { my @b2 = ($first); my @r2; for my $i (0..$n-1) { if (($mask2 >> $i) & 1) { push @b2, $m[$i]; } else { push @r2, $m[$i]; } }
      $cnt++ if join(',', sort { $a <=> $b } @b2) eq $bk; }   # number of index-subsets giving this block (labelled multiplicity)
    for my $sub (mpartitions(@rest)) { push @out, [ $cnt, [ [ split /,/, $bk ], @$sub ] ]; } }
  # each output: [multiplicity, [blocks]]; multiplicity counts labelled choices for the first block; recursion counts the rest
  my @flat; for my $o (@out) { push @flat, $o; } return @flat; }
# labelled set partitions of the multiset (positions distinct): value = prod_blocks 2 F_b ; (b!/b! cancels for labelled positions: d_u^m exp(S) at 0
# equals sum over set partitions of the m positions into even blocks of prod 2 F_b)  [since S has 2F_b u^b/b! and d_u^m picks b! labelled]
sub setparts_even { my @pos = @_; return ([]) unless @pos; my $f = shift @pos; my @out; my $n = scalar @pos;
  for my $mask (0..(1 << $n) - 1) { my @blk = ($f); my @rest; for my $i (0..$n-1) { if (($mask >> $i) & 1) { push @blk, $pos[$i]; } else { push @rest, $pos[$i]; } }
    next if (scalar(@blk) % 2) == 1; for my $sub (setparts_even(@rest)) { push @out, [ \@blk, @$sub ]; } } return @out; }
my %Fpoly;   # equation id -> list of [coef, [multi-index blocks]]
for my $ys (sort keys %H) { my @terms;
  for my $ds (keys %{$H{$ys}}) { my $c = $H{$ys}{$ds}; my @de = split /,/, $ds; my @m; for my $i (0..$nv-1) { push @m, ($i+1) x $de[$i]; }
    next unless @m; next if (scalar(@m) % 2) == 1;   # odd Hirota monomials vanish on tau.tau
    for my $sp (setparts_even(@m)) { push @terms, [ $c * R(2) ** scalar(@$sp), [ map { [ sort { $a <=> $b } @$_ ] } @$sp ] ]; } }
  # merge identical terms
  my %agg; for my $t (@terms) { my $k = join(';', sort map { join(',', @$_) } @{$t->[1]}); $agg{$k} = ($agg{$k} // R(0)) + $t->[0]; }
  my @lst; for my $k (keys %agg) { next if $agg{$k}->is_zero; push @lst, [ $agg{$k}, [ map { [ split /,/, $_ ] } split /;/, $k ] ]; }
  $Fpoly{$ys} = \@lst if @lst; }
printf "Hirota equations up to weight %d: %d non-trivial (y-monomials)\n", $W, scalar keys %Fpoly;
# ---------- 3. dilaton-point values (REFEREE 2 replacement) ----------
# The Hirota polynomials (sections 1-2) are Route A's code verbatim (a3_hirota_check.pl). Here the dilaton-point values Xi_g come either
# from the R2 formulas (MODE=formula: (6.1) for one leaf, Theorem 1 at n = 2, 3 via the referee implementations RefXY/RefXY3, integer tau)
# or from the TR-free engine tables (MODE=engine), selected by the environment variable MODE. Unit parts by the string equation,
# genus 0 by Gamma(J)/Gamma(J-r+3), genus 1 with only unit parts (b-1)!(-(1+1/sigma)/24), genus >= 2 with only unit parts 0 (as in Route A).
use lib "$FindBin::Bin/..";
use RefR2; use RefXY; use RefXY3;
my $MODE = $ENV{MODE} // 'formula';
my $EF = "$FindBin::Bin/../../../../E_formulas"; my %TAB;
sub tab { my ($g, $r) = @_; my $k = "$g,$r"; return $TAB{$k} if exists $TAB{$k};
  for my $f ("$EF/engine_h${g}_r${r}.txt", "$EF/referee/ref_engine_h${g}_r${r}.txt", "$FindBin::Bin/../eng/ref_engine_h${g}_r${r}.txt") {
    if (-e $f) { my $t = RefR2::load_table($f, $g, $r); if (%$t) { $TAB{$k} = $t; return $t; } } } die "no engine table genus $g r=$r"; }
my %memo; my %c61;
sub XiPure { my ($g, @R) = @_; @R = sort { $b <=> $a } @R; my $key = "$g|@R"; return $memo{$key} if exists $memo{$key}; my $r = scalar @R; my $xi;
  if ($MODE eq 'engine') { $xi = RefR2::table_value(tab($g, $r), $tau, @R); }
  elsif ($r == 1) { $xi = RefR2::xi61($g, $tau, $R[0], \%c61); }
  elsif ($r == 2) { $xi = RefXY::xi2($tau, $g, @R); }
  elsif ($r == 3) { $xi = RefXY3::xi3($tau, $g, reverse @R); }
  else { die "n = $r > 3 not implemented (lower W)"; }
  $memo{$key} = $xi; return $xi; }
sub Xi { my ($g, @L) = @_; my $J = 0; $J += $_ for @L; my $r = scalar @L; die "empty" unless $r;
  if ($g == 0) { my $v = R(1); if ($r >= 3) { $v *= ($J - $_) for 1..($r-3); } elsif ($r == 2) { $v = R(1)/$J; } else { $v = R(1)/($J*($J+1)); } return $v; }
  my @Tp = grep { $_ >= 2 } @L; my $b = $r - scalar @Tp; my $Jp = 0; $Jp += $_ for @Tp;
  return XiPure($g, @Tp) if $b == 0;
  if (@Tp) { my $f = R(1); $f *= ($Jp + $_) for 0..$b-1; return $f * XiPure($g, @Tp); }
  return fact($b-1) * (-(1 + 1/$sig)/24) if $g == 1; return R(0); }
sub Pi { my $j = shift; return R(1) if $j <= 1; my $p = R(1); for my $b (1..$j-1) { $p *= ($j*$tau + $b); } return $p / fact($j-1); }
sub factor_val { my ($g, $Rj, $b) = @_; my @K = grep { $_ >= 2 } @$b; my $a = scalar(@$b) - scalar(@K);
  my $v = $sig ** ($g + scalar(@$b) - 1); $v *= $_ for @$b; $v *= Pi($_) for @K;
  return $v * Xi($g, @$Rj, @K, (1) x $a); }
my ($NZ, $BAD) = (0, 0);
for my $R (@Rs) { my @R = @$R; my $r = scalar @R; my ($nz, $bad, $nonzero_terms) = (0, 0, 0);
  for my $ys (sort keys %Fpoly) { my $tot = R(0); my $nterm = 0; my $nzt = 0;
    for my $t (@{$Fpoly{$ys}}) { my ($c, $blocks) = @$t; my $k = scalar @$blocks;
      my $need = 2*$G; my @base = map { scalar(@$_) - 2 } @$blocks; my $sb = 0; $sb += $_ for @base; my $rem = $need - $sb; next if $rem < 0 || $rem % 2;
      my $gsum = $rem / 2;
      my @gas; my $rec; $rec = sub { my ($i, $left, $acc) = @_; if ($i == $k) { push @gas, [@$acc] if $left == 0; return; } for my $gg (0..$left) { $rec->($i+1, $left-$gg, [@$acc, $gg]); } }; $rec->(0, $gsum, []);
      for my $ga (@gas) {
        my @assign = ([]); for my $part (@R) { my @na; for my $as (@assign) { for my $j (0..$k-1) { push @na, [@$as, $j]; } } @assign = @na; }
        for my $as (@assign) { my @Rj = map { [] } 1..$k; for my $q (0..$r-1) { push @{$Rj[$as->[$q]]}, $R[$q]; }
          my $v = $c; for my $j (0..$k-1) { $v *= factor_val($ga->[$j], $Rj[$j], $blocks->[$j]); last if $v->is_zero; } $tot += $v; $nterm++; $nzt++ unless $v->is_zero; } } }
    $nz++ if $nterm; $bad++ unless $tot->is_zero; $nonzero_terms += $nzt;
    print "  R=(@R) eq y^[$ys]: total = $tot\n" unless $tot->is_zero; }
  $NZ += $nz; $BAD += $bad;
  printf "MODE=%s G=%d tau=%s R=(%s): %d Hirota equations (all genera 0..G), %d non-zero individual terms, %d non-zero totals\n", $MODE, $G, $tau, join(',', @R), $nz, $nonzero_terms, $bad; }
printf "SUMMARY MODE=%s G=%d tau=%s W=%d: %d equations, %d failures\n", $MODE, $G, $tau, $W, $NZ, $BAD;
