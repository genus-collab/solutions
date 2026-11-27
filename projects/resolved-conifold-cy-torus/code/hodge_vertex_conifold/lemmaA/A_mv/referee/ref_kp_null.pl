#!/usr/bin/perl
# ref_kp_null.pl -- referee's independent implementation of Theorem D (linearised dispersionless Hirota at the dilaton point) and of the
# rank computations behind Theorem F.  Written from the referee's own derivation; shares no code with the route's a3_kp_nullity*.pl.
# Relation (R, al, be), al, be >= 2, e := 2G - D(R) - al - be >= 1, R labelled (positions 0..n-1):
#   sum_{a<=al, b<=be} sum_{K subset R} Pi_a Pi_b Xi_G(K + {a,b}) E[al-a][be-b][R\K]  =  Pi_c Xi_G(R + {1,c}),   c = al+be-1,
#   E = exp(Om),  Om[a][b][S] = sigma Pi_a Pi_b Xi_0(S + {a,b})  (a, b >= 1),  Xi_0(L) = Gamma(J)/Gamma(J-r+3) (unit parts included),
#   unit parts of genus-G arguments by the string equation, Xi_G(1^u) = 0.
#   E is computed with the w-degree Euler recursion  n E[m][n][S] = sum b1 Om[a1][b1][K] E[m-a1][n-b1][S\K]  (labelled Leibniz rule).
# Modes:
#   selftest G p seed      : genus-0 identity [z^-al w^-be] E[R] = Om[1][al+be-1][R] (all al, be >= 1) for all R up to D(R) <= 2G-3
#   full G p seed S        : unknowns c_{T,k} (sigma^G Xi_G(T) = sum_{k=1..G} c_{T,k} sigma^k) for ALL T in the range; all relations
#   twoleaf G p seed S     : unknowns = one-leaf c_{t,k}; two-leaf values from the R = empty relations; polynomiality of every two-leaf
#                            T imposed by Lagrange interpolation in sigma on the space sigma*Q[sigma]_{<=G-1} (first G samples as nodes)
use strict; use warnings;
my ($mode, $G, $p, $seed, $S) = @ARGV; srand($seed); $S //= 0;
sub md { my $x = $_[0] % $p; return $x < 0 ? $x + $p : $x; }
sub mi { my $a = md($_[0]); die "not invertible" unless $a; my ($t, $nt, $r, $nr) = (0, 1, $p, $a); while ($nr) { my $q = int($r / $nr); ($t, $nt) = ($nt, $t - $q*$nt); ($r, $nr) = ($nr, $r - $q*$nr); } return md($t); }
sub mpartsD { my ($D, $mx) = @_; return ([]) if $D == 0; my @o; for my $k (reverse 1..($mx < $D ? $mx : $D)) { push @o, [$k, @$_] for mpartsD($D - $k, $k); } return @o; }
sub tkey { return join('.', sort { $b <=> $a } @_); }
my $Dmax = 2*$G - 3;
# multisets T in the range (parts >= 2), indexed
my @T; my %Tix; for my $D (1..$Dmax) { for my $q (mpartsD($D, $D)) { my @t = map { $_ + 1 } @$q; $Tix{tkey(@t)} = scalar @T; push @T, \@t; } }
my ($fac, $ifac) = ([1], [1]); for my $i (1..200) { $fac->[$i] = ($fac->[$i-1] * $i) % $p; $ifac->[$i] = mi($fac->[$i]); }
sub X0 { my @L = @_; my $J = 0; $J += $_ for @L; my $r = scalar @L; die "X0 empty" unless $r;   # Gamma(J)/Gamma(J-r+3)
  if ($r == 1) { return mi($J * ($J + 1)); } if ($r == 2) { return mi($J); } my $v = 1; for my $i (1..$r-3) { $v = ($v * md($J - $i)) % $p; } return $v; }
# per-sample context
sub ctx { my ($tau) = @_; my %c = (tau => $tau, sig => ($tau * md($tau + 1)) % $p); my @Pi = (1, 1);
  for my $j (2..4*$G + 4) { my $v = 1; for my $b (1..$j-1) { $v = ($v * md($j * $tau + $b)) % $p; } $Pi[$j] = ($v * $ifac->[$j-1]) % $p; } $c{Pi} = \@Pi; return \%c; }
# E table for labelled R: returns closure E(m, n, mask)
sub Etab { my ($cx, $R, $amax, $bmax) = @_; my $n = scalar @$R; my $full = (1 << $n) - 1; my ($sig, $Pi) = ($cx->{sig}, $cx->{Pi});
  my @parts = map { my $m = $_; [ map { $R->[$_] } grep { ($m >> $_) & 1 } 0..$n-1 ] } 0..$full;
  my %Om; for my $mask (0..$full) { for my $a (1..$amax) { for my $b (1..$bmax) { $Om{"$a,$b,$mask"} = ($sig * (($Pi->[$a] * $Pi->[$b]) % $p) % $p) * X0(@{$parts[$mask]}, $a, $b) % $p; } } }
  my %E; my $sub; $sub = sub { my ($m, $nn, $mask) = @_; my $k = "$m,$nn,$mask"; return $E{$k} if exists $E{$k}; my $res = 0;
    if ($nn == 0) { $res = ($m == 0 && $mask == 0) ? 1 : 0; }
    elsif ($m <= 0) { $res = 0; }
    else { for my $a1 (1..$m) { for my $b1 (1..$nn) { my $K = $mask; while (1) { my $e = $sub->($m - $a1, $nn - $b1, $mask & ~$K);
              $res = ($res + (($b1 * $Om{"$a1,$b1,$K"}) % $p) * $e) % $p if $e; last if $K == 0; $K = ($K - 1) & $mask; } } }
      $res = ($res * mi($nn)) % $p; }
    return $E{$k} = $res; };
  return ($sub, \%Om, \@parts); }
if ($mode eq 'selftest') { my ($n, $bad) = (0, 0);
  for my $it (1..2) { my $cx = ctx(3 + int(rand($p - 9)));
    my @Rs = ([]); for my $D (1..$Dmax) { push @Rs, map { [ map { $_ + 1 } @$_ ] } mpartsD($D, $D); }
    for my $R (@Rs) { my $DR = 0; $DR += $_ - 1 for @$R; next if scalar(@$R) > 5; my $w = 2*$G + 1 - $DR; next if $w < 2;
      my ($E, $Om, $parts) = Etab($cx, $R, $w, $w); my $full = (1 << scalar(@$R)) - 1;
      for my $al (1..$w) { for my $be (1..$w) { next if $al + $be > $w + 1; my $c = $al + $be - 1;
          my $rhs = ($cx->{sig} * $cx->{Pi}[$c] % $p) * X0(@$R, 1, $c) % $p; my $lhs = $E->($al, $be, $full); $n++; $bad++ if $lhs != $rhs; } } } }
  print "selftest G=$G: genus-0 dispersionless Hirota at the dilaton point, $n relations, $bad mismatches\n"; exit; }
# ---- linear forms: hash col => coef ----
sub XGform { my ($cx, $ncol, @L) = @_;   # Xi_G(L) as linear form; full mode: columns T*G + k-1 ; one-leaf mode handled separately
  my @Tp = grep { $_ >= 2 } @L; my $u = scalar(@L) - scalar(@Tp); return {} unless @Tp; my $J = 0; $J += $_ for @Tp; my $f = 1; $f = ($f * md($J + $_)) % $p for 0..$u-1;
  my $i = $Tix{tkey(@Tp)}; die "out of range: @L" unless defined $i; my $isg = mi($cx->{sig}); my $sp = 1; for (1..$G) { $sp = ($sp * $isg) % $p; }   # sigma^-G
  my %v; my $pw = 1; for my $k (1..$G) { $pw = ($pw * $cx->{sig}) % $p; $v{$i*$G + $k - 1} = ((($f * $sp) % $p) * $pw) % $p; } return \%v; }
sub addto { my ($acc, $v, $c) = @_; for my $k (keys %$v) { $acc->{$k} = (($acc->{$k} // 0) + $c * $v->{$k}) % $p; } }
my @rows;
sub rank_of { my ($rows, $ncol) = @_; my %piv; my $rank = 0;
  for my $row (@$rows) { last if $rank == $ncol; my %v = map { $_ => $row->{$_} } grep { $row->{$_} } keys %$row; next unless %v;
    while (%v) { my ($c0) = sort { $a <=> $b } keys %v; if (exists $piv{$c0}) { my $f = $v{$c0}; my $pr = $piv{$c0};
        for my $c2 (keys %$pr) { my $nv = md(($v{$c2} // 0) - $f * $pr->{$c2}); if ($nv) { $v{$c2} = $nv; } else { delete $v{$c2}; } } next; }
      my $inv = mi($v{$c0}); my %nr = map { $_ => ($v{$_} * $inv) % $p } keys %v; $piv{$c0} = \%nr; $rank++; last; } }
  return ($rank, \%piv); }
if ($mode eq 'full') { $S ||= 3*$G + 6; my $ncol = scalar(@T) * $G;
  my @Rs = ([]); for my $D (1..2*$G - 5) { push @Rs, map { [ map { $_ + 1 } @$_ ] } mpartsD($D, $D); }
  for my $s (1..$S) { my $cx = ctx(3 + int(rand($p - 9)));
    for my $R (@Rs) { my $DR = 0; $DR += $_ - 1 for @$R; my $wmax = 2*$G - 1 - $DR; next if $wmax < 4; my $n = scalar @$R; my $full = (1 << $n) - 1;
      my ($E, $Om, $parts) = Etab($cx, $R, $wmax - 2, $wmax - 2);
      for my $al (2..$wmax - 2) { for my $be ($al..$wmax - $al) { my %row;
          for my $a (1..$al) { for my $b (1..$be) { my $K = $full; while (1) { my $e = $E->($al - $a, $be - $b, $full & ~$K);
                if ($e) { my $v = XGform($cx, $ncol, @{$parts->[$K]}, $a, $b); addto(\%row, $v, ((($cx->{Pi}[$a] * $cx->{Pi}[$b]) % $p) * $e) % $p); }
                last if $K == 0; $K = ($K - 1) & $full; } } }
          my $c = $al + $be - 1; addto(\%row, XGform($cx, $ncol, @$R, 1, $c), md(-$cx->{Pi}[$c])); push @rows, \%row; } } } }
  my ($rank) = rank_of(\@rows, $ncol);
  printf "full G=%d: %d multisets in range, %d unknowns, %d relations (%d samples, p=%d): rank %d, NULLITY %d\n", $G, scalar @T, $ncol, scalar @rows, $S, $p, $rank, $ncol - $rank; exit; }
if ($mode eq 'twoleaf') { $S ||= 2*$G + 6; my $nc = ($Dmax) * $G;   # one-leaf t = 2..2G-2 -> index (t-2)*G + k-1
  my @pairs; for my $sum (4..2*$G - 1) { for my $b (2..int($sum/2)) { push @pairs, [$sum - $b, $b]; } }   # (al >= be)
  my @samp;   # per sample: sigma, and the forms of sigma^G Xi_G(al,be)
  for my $s (1..$S) { my $cx = ctx(3 + int(rand($p - 9))); my $sig = $cx->{sig}; my $Pi = $cx->{Pi}; my $w = 2*$G - 1;
    my ($E) = Etab($cx, [], $w, $w);
    my $one = sub { my $t = shift; my %v; my $pw = 1; for my $k (1..$G) { $pw = ($pw * $sig) % $p; $v{($t-2)*$G + $k - 1} = $pw; } return \%v; };   # sigma^G Xi_G(t)
    my %two;   # key "a,b" (a>=b>=2) -> form of sigma^G Xi_G(a,b)
    my $X = sub { my ($a, $b) = @_; ($a, $b) = ($b, $a) if $b > $a; return {} if $a == 1; if ($b == 1) { my $v = $one->($a); return { map { $_ => ($v->{$_} * $a) % $p } keys %$v }; }
      return $two{"$a,$b"} // die "two-leaf ($a,$b) not yet known"; };
    for my $pr (@pairs) { my ($al, $be) = @$pr; my %rhs; my $c = $al + $be - 1; addto(\%rhs, $one->($c), ($Pi->[$c] * $c) % $p);
      for my $a (1..$al) { for my $b (1..$be) { next if $a == $al && $b == $be; my $e = $E->($al - $a, $be - $b, 0); next unless $e;
          addto(\%rhs, $X->($a, $b), md(-(($Pi->[$a] * $Pi->[$b]) % $p) * $e)); } }
      my $inv = mi(($Pi->[$al] * $Pi->[$be]) % $p); $two{"$al,$be"} = { map { $_ => ($rhs{$_} * $inv) % $p } keys %rhs }; }
    push @samp, [ $sig, \%two ]; }
  # Lagrange basis on nodes sigma_1..sigma_G for the space sigma*Q[sigma]_{<=G-1}: l_s(x) = (x/sig_s) prod_{u != s} (x - sig_u)/(sig_s - sig_u)
  my @nodes = map { $samp[$_][0] } 0..$G-1;
  for my $t ($G..$S-1) { my $x = $samp[$t][0]; my @l;
    for my $s (0..$G-1) { my $v = ($x * mi($nodes[$s])) % $p; for my $u (0..$G-1) { next if $u == $s; $v = ($v * md($x - $nodes[$u])) % $p; $v = ($v * mi(md($nodes[$s] - $nodes[$u]))) % $p; } push @l, $v; }
    for my $pr (@pairs) { my $k = "$pr->[0],$pr->[1]"; my %row; addto(\%row, $samp[$t][1]{$k}, 1); addto(\%row, $samp[$_][1]{$k}, md(-$l[$_])) for 0..$G-1; push @rows, \%row; } }
  my ($rank, $piv) = rank_of(\@rows, $nc);
  printf "twoleaf G=%d: %d one-leaf unknowns, %d two-leaf multisets, %d equations (%d samples, p=%d): rank %d, NULLITY %d\n", $G, $nc, scalar @pairs, scalar @rows, $S, $p, $rank, $nc - $rank; exit; }
die "unknown mode $mode\n";
