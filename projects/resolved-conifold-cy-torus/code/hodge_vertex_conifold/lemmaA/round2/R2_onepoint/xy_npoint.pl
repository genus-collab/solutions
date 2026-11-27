#!/usr/bin/perl
# xy_npoint.pl tau g n jmax [engine_file] : n-leaf values from the ABDKS x-y swap formula (LogTR with affine y, arXiv:2312.16950,
# Prop. 3.1) on C_tau (integer tau >= 1), after integration by parts:
#   [x^mu] F_{g,n} = (-1)^{n-1} (prod 1/mu_i) [hb^{2g-2+n}] [v_1^{-1}] ... [v_n^{-1}] ( prod_i Q_{mu_i}(v_i) * Sum_{sigma n-cycle} prod_i 1/(v_i - v_sigma(i) + k_{i,sigma(i)} hb) ),
#   k_{ij} = (mu_i + mu_j)/2,  Q_mu(v) = prod_{M2(mu)} (v + 1 + m' hb/e2) / prod_{M1(mu)} (v + m hb/tau)  (|M1| = tau mu, |M2| = e2 mu),
# Laurent expansions at infinity in the region |v_1| >> ... >> |v_n| (v_1 extracted first).  n = 1: the 1-cycle factor is 1/(mu hb).
# Dictionary P2: Xi_g(mu) = sigma^{g-1} [x^mu] F_{g,n} / prod Pi_{mu_i}.  Compared with Route E's engine_h<g>_r<n>.txt.
# Also prints the homogeneity check: the coefficient of hb^p in the extracted residue is non-zero only for p <= |mu|.
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/lib";
use Math::BigRat; use MP;
my ($tau, $g, $n, $jmax, $datafile) = @ARGV; $tau //= 1; $g //= 2; $n //= 2; $jmax //= 4;
$datafile //= "$FindBin::Bin/../../E_formulas/engine_h${g}_r${n}.txt";
$tau = Q($tau); my $e2 = $tau + 1; my $sig = $tau * $e2;
my $H = 2*$g - 2 + $n;     # target hb-power
my $HT = $H + 1;           # truncation (one more, to see the homogeneity bound)
sub fact { my $k = shift; my $r = Q(1); $r *= $_ for 2..$k; return $r; }
sub Pi { my $j = shift; my $p = Q(1); $p *= ($j * $tau + $_) for 1..$j-1; return $p / fact($j-1); }
sub qser { my ($mu, $T) = @_; my $n1 = $tau * $mu; my $n2 = $e2 * $mu;
  my @P = map { [ map { Q(0) } 0..$HT ] } 0..$T; $P[0][0] = Q(1);
  my $mulf = sub { my ($F) = @_; my @C = map { [ map { Q(0) } 0..$HT ] } 0..$T;
    for my $i (0..$T) { for my $k (0..$T-$i) { for my $p (0..$HT) { next if $P[$i][$p]->is_zero; for my $q (0..$HT-$p) { next if $F->[$k][$q]->is_zero; $C[$i+$k][$p+$q] += $P[$i][$p] * $F->[$k][$q]; } } } }
    @P = @C; };
  for my $kk (1..$n2) { my $m = Q($kk) - ($n2 + 1)/2; my @F = map { [ map { Q(0) } 0..$HT ] } 0..$T; $F[0][0] = Q(1); if ($T >= 1) { $F[1][0] = Q(1); $F[1][1] = $m / $e2; } $mulf->(\@F); }
  for my $kk (1..$n1) { my $m = Q($kk) - ($n1 + 1)/2; my @F = map { [ map { Q(0) } 0..$HT ] } 0..$T; my $cf = Q(1);
    for my $d (0..$T) { $F[$d][$d] = $cf if $d <= $HT; $cf *= -$m / $tau; } $mulf->(\@F); }
  return \@P; }
# polynomials in variables v_1..v_n, hb  (index n = hb)
my $NV = $n + 1;
sub trunc_h { my ($p) = @_; my %q; for my $k (keys %$p) { my @e = split /,/, $k; next if $e[$n] > $HT; $q{$k} = $p->{$k}; } return \%q; }
sub hpow_poly { my ($arr) = @_; my %p; for my $d (0..$#$arr) { next if $arr->[$d]->is_zero; my @e = (0) x $NV; $e[$n] = $d; $p{join(',', @e)} = $arr->[$d]->copy; } return \%p; }
# n-cycles of {0..n-1}: permutations with a single cycle
sub ncycles { my ($k) = @_; return ([0]) if $k == 1; my @out;
  my $rec; $rec = sub { my ($path, $rest) = @_; if (!@$rest) { my @s; for my $i (0..$#$path) { $s[$path->[$i]] = $path->[($i+1) % @$path]; } push @out, \@s; return; }
    for my $i (0..$#$rest) { my @r2 = @$rest; my ($x) = splice(@r2, $i, 1); $rec->([@$path, $x], \@r2); } };
  $rec->([0], [1..$k-1]); return @out; }
my @cycles = ncycles($n);
my %C; if (-e $datafile) { open my $fh, '<', $datafile or die; while (my $line = <$fh>) { chomp $line; my @f = split /\s*\|\s*/, $line; next unless @f >= 4; my ($h, $r, $k, $rest) = @f; next unless $h == $g && $r == $n;
    for my $term (split /\s*;\s*/, $rest) { next unless $term =~ /^([\d,]+):(\S+)$/; $C{$k}{$1} = Q($2); } } close $fh; }
sub xi_engine { my @j = @_; return undef unless %C; my $v = Q(0); my $vv = -Q(1) / $sig;
  for my $k (keys %C) { my $s = Q(0); for my $key (keys %{$C{$k}}) { my @a = split /,/, $key; my $t = $C{$k}{$key}->copy; $t *= Q($j[$_])**$a[$_] for 0..$#a; $s += $t; } $v += $s * $vv**$k; } return $v; }
sub residue { my (@mu) = @_; my $sum = 0; $sum += $_ for @mu; my $T = 2*$sum + 3*$n + 3;
  my @q = map { qser($_, $T) } @mu;
  my @tot = map { Q(0) } 0..$HT;
  if ($n == 1) { # [v^{-1}] Q_mu(v) / (mu hb): coefficient of v^{-1} is q_{mu+1}; divide by mu hb
    my $mu = $mu[0]; for my $d (1..$HT) { $tot[$d-1] += $q[0][$mu+1][$d] / $mu; } return \@tot; }
  for my $sg (@cycles) {
    # factors: for each i, 1/(v_i - v_s(i) + k hb); orient with smaller index first: if i < s: (v_i - v_s + k hb); else -(v_s - v_i - k hb)
    my @fac; my $sign = 1;
    for my $i (0..$n-1) { my $s = $sg->[$i]; my $k = Q($mu[$i] + $mu[$s]) / 2;
      if ($i < $s) { push @fac, [$i, $s, $k]; } else { push @fac, [$s, $i, -$k]; $sign = -$sign; } }   # factor 1/(v_a - v_b + c hb)
    my $P = { join(',', (0) x $NV) => Q($sign) };
    my @remaining = @fac;
    for my $step (0..$n-1) {
      # factors with a == step: 1/(v_step - (v_b - c hb)) ; v_step is the largest remaining variable
      my @here = grep { $_->[0] == $step } @remaining; @remaining = grep { $_->[0] != $step } @remaining;
      die "unexpected factor" if grep { $_->[1] == $step } @remaining;
      my $m = scalar @here;
      my @A = map { my $f = $_; my $p = mp_var($f->[1], $NV); $p = mp_add($p, mp_scal(-$f->[2], mp_var($n, $NV))); $p } @here;   # A = v_b - c hb
      my $emax = 0; for my $key (keys %$P) { my @e = split /,/, $key; $emax = $e[$step] if $e[$step] > $emax; }
      my $smax = $emax + $mu[$step] + 1;
      # complete homogeneous h_s(A_1..A_m) for s <= smax
      my @hs;
      if ($m == 0) { @hs = ((mp_const(1,$NV)), map { {} } 1..$smax); }
      elsif ($m == 1) { @hs = (mp_const(1,$NV)); for my $s (1..$smax) { push @hs, trunc_h(mp_mul($hs[-1], $A[0])); } }
      else { my @p0 = (mp_const(1,$NV)); my @p1 = (mp_const(1,$NV)); for my $s (1..$smax) { push @p0, trunc_h(mp_mul($p0[-1], $A[0])); push @p1, trunc_h(mp_mul($p1[-1], $A[1])); }
        for my $s (0..$smax) { my $h = {}; for my $a (0..$s) { $h = mp_add($h, trunc_h(mp_mul($p0[$a], $p1[$s-$a]))); } push @hs, $h; } }
      my $NP = {};
      for my $key (keys %$P) { my @e = split /,/, $key; my $ee = $e[$step]; $e[$step] = 0; my $mono = { join(',', @e) => $P->{$key} };
        # [v^{-1}] v^ee * v^mu * sum_a q_a v^{-a} * sum_s h_s v^{-s-m}:  ee + mu - a - s - m = -1
        for my $a (0..($ee + $mu[$step] + 1 - $m)) { my $s = $ee + $mu[$step] + 1 - $m - $a; next if $s < 0 || $s > $smax; next if $a > $T;
          my $qa = hpow_poly($q[$step][$a]); next unless %$qa; next unless %{$hs[$s]};
          $NP = mp_add($NP, trunc_h(mp_mul(mp_mul($mono, $qa), $hs[$s]))); } }
      $P = $NP; }
    for my $key (keys %$P) { my @e = split /,/, $key; die "leftover variable" if grep { $e[$_] } 0..$n-1; $tot[$e[$n]] += $P->{$key}; } }
  return \@tot; }
my ($nb, $nt, $nh) = (0, 0, 0);
my @tuples; my $rec; $rec = sub { my ($pre, $max) = @_; if (@$pre == $n) { push @tuples, [@$pre]; return; } for my $x (reverse 1..$max) { $rec->([@$pre, $x], $x); } }; $rec->([], $jmax);
for my $mu (@tuples) { my @mu = @$mu; my $sum = 0; $sum += $_ for @mu;
  my $res = residue(@mu);
  # homogeneity: coefficient of hb^p must vanish for p > |mu|
  my $hviol = 0; for my $p (0..$HT) { $hviol++ if $p > $sum && !$res->[$p]->is_zero; } $nh += $hviol;
  my $pm = Q(1); $pm *= $_ for @mu;
  my $F = ($n % 2 ? 1 : -1) * $res->[$H] / $pm;
  my $pp = Q(1); $pp *= Pi($_) for @mu;
  my $xi = $sig**($g-1) * $F / $pp;
  my $xe = xi_engine(@mu);
  my $tag = defined $xe ? ($xi == $xe ? "OK" : "MISMATCH") : "(no data)";
  $nt++; $nb++ if $tag eq "MISMATCH";
  printf "tau=%s g=%d mu=(%s): xy %s  engine %s  %s%s%s\n", $tau, $g, join(',', @mu), $xi, (defined $xe ? $xe : '-'), $tag,
    ($sum <= 2*$g - 3 + $n ? "  [range]" : ""), ($hviol ? "  HOMOGENEITY-VIOLATION" : ""); }
print "total $nt, mismatches $nb, homogeneity violations $nh\n";
