#!/usr/bin/perl
# r2_lemma5.pl TAU G N KMAX : check of Lemma 5 (dilaton resummation) with the Route E referee's TR-free engine RefHodge,
# dilaton/v leaves kept as explicit marked points:
#   Sum_K ((e sigma)^K/K!) <Omega prod psi^{k_i} prod_{j<=K} v(psi_{n+j})>_{g,n+K} = (1+e sigma)^(2-2g-n) H_{g,n}(e sigma/(1+e sigma); k),
#   v(psi) = -psi/(1-psi),  H(s;k) = Sum_m (s^m/m!) <Omega prod psi^{k_i} prod T(psi)>,  T = -psi^2/(1-psi),
# compared as power series in e through order e^KMAX, for every k with |k| <= 3g-3+n.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub R { Math::BigRat->new(@_) }
my ($tau, $g, $n, $KMAX) = @ARGV; $tau = R($tau); my $sig = $tau * (1 + $tau);
my $dimmax = 3*$g - 3 + $n + $KMAX;
my $om = RefHodge::omega($g, [R(1), $tau, -1 - $tau], ($g == 1 ? 1 : 3*$g - 3));
sub fact { my $m = shift; my $r = R(1); $r *= $_ for 2..$m; $r }
sub binI { my ($a, $b) = @_; return R(0) if $b < 0 || $b > $a; fact($a) / (fact($b) * fact($a - $b)) }
# <Omega prod psi^{k_i} prod_{j} psi^{d_j}> summed over d_j >= dmin, for m extra points: sum over compositions
sub extra { my ($k, $m, $dmin) = @_; my $K = 0; $K += $_ for @$k; my $dim = 3*$g - 3 + $n + $m; my $tot = R(0);
  my @mus = keys %$om;
  my $rec; $rec = sub { my ($ds, $left) = @_;
    if (@$ds == $m) { my $s = R(0); for my $mu (@mus) { my @mm = split /,/, $mu; my $dm = 0; my @CH = @RefHodge::CHL; $dm += $CH[$_] * $mm[$_] for 0..$#mm;
          next unless $dm + $K + (eval { my $t = 0; $t += $_ for @$ds; $t }) == $dim; my $v = RefHodge::I($g, \@mm, [@$k, @$ds]); $s += $om->{$mu} * $v unless $v->is_zero; }
      $tot += $s; return; }
    for my $d ($dmin..$left) { $rec->([@$ds, $d], $left - $d); } };
  $rec->([], $dim - $K); return $tot; }
my ($cmp, $bad) = (0, 0);
my @ks; my $tr; $tr = sub { my ($pre, $left) = @_; if (@$pre == $n) { push @ks, [@$pre]; return; } for my $x (0..$left) { $tr->([@$pre, $x], $left - $x); } }; $tr->([], 3*$g - 3 + $n);
for my $k (@ks) {
  # LHS coefficients of e^K
  my @L = map { R(0) } 0..$KMAX; for my $KK (0..$KMAX) { my $v = extra($k, $KK, 1); $v *= R(-1) ** $KK; $L[$KK] = $v * $sig ** $KK / fact($KK); }
  # H coefficients h_m = (1/m!) <... T^m> = ((-1)^m/m!) Sum_{d>=2}
  my $Mmax = 3*$g - 3 + $n; my @h = map { extra($k, $_, 2) * R(-1) ** $_ / fact($_) } 0..$Mmax;
  # RHS: Sum_m h_m (e sigma)^m (1+e sigma)^(2-2g-n-m), expanded in e
  my @Rt = map { R(0) } 0..$KMAX; my $N0 = 2 - 2*$g - $n;
  for my $m (0..$Mmax) { next if $h[$m]->is_zero; my $p = $N0 - $m;   # (1+x)^p = Sum_i C(p,i) x^i, x = e sigma
    for my $i (0..$KMAX - $m) { my $c = R(1); for my $t (0..$i-1) { $c *= R($p - $t) / ($t + 1); } $Rt[$m + $i] += $h[$m] * $c * $sig ** ($m + $i); } }
  for my $KK (0..$KMAX) { $cmp++; $bad++ if $L[$KK] != $Rt[$KK]; } }
print "Lemma 5 at tau = $tau (sigma = $sig), (g,n) = ($g,$n), e-orders 0..$KMAX, all k: $cmp comparisons, $bad mismatches\n";
